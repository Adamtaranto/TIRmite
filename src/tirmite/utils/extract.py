"""
Unified sequence extraction from indexed FASTA or BLAST databases.

TIRmite can retrieve hit sequences from either a ``pyfaidx``-indexed FASTA
genome or a BLAST database.  These are user-selectable alternatives, so the
same coordinates must yield the same sequence from either source.  This module
provides the single primitive both backends go through.

Notes
-----
**Coordinate contract.** All coordinates in this module are 1-based, inclusive,
plus-strand genomic coordinates, matching ``blastdbcmd -range`` and the
``hitStart``/``hitEnd`` columns produced by
:func:`tirmite.core.parsers.import_nhmmer` and
:func:`tirmite.core.parsers.import_blast`.  Conversion to ``pyfaidx``'s 0-based
half-open slicing happens inside :class:`FastaSource` and nowhere else.

**Return contract.** :func:`fetch_sequence` returns uppercase, plus-strand
sequence, or ``None`` if the region could not be retrieved in full.
Reverse-complementing is the caller's responsibility and must be applied
*after* any case-marking of flanking sequence, otherwise the case markers land
on the wrong end.

**Backend differences neutralised here.** ``blastdbcmd`` and ``pyfaidx``
disagree in several ways that would otherwise produce different sequences for
the same query:

- ``blastdbcmd`` silently returns the **entire sequence** when the requested
  range starts beyond the end of the contig; ``pyfaidx`` returns an empty
  string.  Both are rejected here by clamping first and verifying the returned
  length.
- ``blastdbcmd`` errors when a range coordinate is ``<= 0``; ``pyfaidx`` treats
  negative indices as offsets from the end.  Both are prevented by clamping.
- ``blastdbcmd`` discards soft-masking (output is always uppercase) while
  ``pyfaidx`` preserves it.  Both are uppercased here.
- A missing sequence ID raises ``KeyError`` from ``pyfaidx`` but is a non-zero
  exit from ``blastdbcmd``.  Both become ``None``.
"""

import logging
from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, NamedTuple, Optional, Tuple, Union

# Shared timeout for all blastdbcmd invocations, in seconds.  Previously these
# ranged from 30s to no timeout at all depending on the call site.
BLASTDBCMD_TIMEOUT = 60


class FastaSource:
    """
    Sequence source backed by a ``pyfaidx``-indexed FASTA genome.

    Parameters
    ----------
    genome : pyfaidx.Fasta
        Indexed genome object, as returned by
        :func:`tirmite.utils.utils.indexGenome`.
    descriptions : dict, optional
        Mapping of sequence ID to FASTA header description, used to decorate
        output records.  Defaults to an empty mapping.

    Notes
    -----
    Sequence IDs are matched exactly and case-sensitively against the keys of
    the ``pyfaidx`` index, which are the first whitespace-delimited token of
    each FASTA header.
    """

    def __init__(self, genome: Any, descriptions: Optional[Dict[str, str]] = None):
        self.genome = genome
        self.descriptions = descriptions or {}

    def contig_length(self, seqid: str) -> Optional[int]:
        """
        Return the length of a sequence in the index.

        Parameters
        ----------
        seqid : str
            Sequence identifier.

        Returns
        -------
        int or None
            Length in bases, or None if ``seqid`` is not in the index.
        """
        try:
            return len(self.genome[seqid])
        except (KeyError, TypeError):
            return None

    def fetch_raw(self, seqid: str, start: int, end: int) -> Optional[str]:
        """
        Fetch a pre-clamped 1-based inclusive region.

        Parameters
        ----------
        seqid : str
            Sequence identifier.
        start, end : int
            1-based inclusive coordinates, already clamped to the contig.

        Returns
        -------
        str or None
            Raw sequence as stored, or None if the read failed.

        Notes
        -----
        Callers should use :func:`fetch_sequence` rather than calling this
        directly; it performs no clamping or validation.
        """
        try:
            return str(self.genome[seqid][start - 1 : end])
        except (KeyError, TypeError) as e:
            logging.error(
                f'Failed to read {seqid}:{start}-{end} from genome index: {e}'
            )
            return None

    def sequence_description(self, seqid: str) -> str:
        """
        Return the FASTA header description, or '' if there is none.

        Parameters
        ----------
        seqid : str
            Sequence identifier.

        Returns
        -------
        str
            Text following the sequence ID in the header, or ''.

        Notes
        -----
        Falls back to the header stored in the index when no description
        mapping was supplied, so this source reports the same descriptions as
        :class:`BlastDBSource` reads from a BLAST database built from the same
        FASTA.
        """
        if seqid in self.descriptions:
            return self.descriptions[seqid]

        # pyfaidx records expose the full header as long_name; the description
        # is whatever follows the first whitespace.
        long_name = getattr(self.genome[seqid], 'long_name', None) if seqid else None
        if not long_name:
            return ''
        parts = long_name.split(None, 1)
        return parts[1].strip() if len(parts) > 1 else ''

    def describe(self) -> str:
        """
        Describe this source for log messages.

        Returns
        -------
        str
            Human-readable identifier for this source.
        """
        return f'genome index {getattr(self.genome, "filename", "<unknown>")}'


class BlastDBSource:
    """
    Sequence source backed by a BLAST database, read via ``blastdbcmd``.

    Parameters
    ----------
    path : str or Path
        Path to the BLAST database (without file extension).  The database must
        have been created with ``-parse_seqids`` for sequence ID lookup to work.
    descriptions : dict, optional
        Mapping of sequence ID to description, used to decorate output records.
        Defaults to an empty mapping.

    Notes
    -----
    ``blastdbcmd -entry`` matches against database *accessions*, which are not
    always identical to the ``pyfaidx`` key for the same header: a header of
    ``>lcl|contig_1`` indexes as ``lcl|contig_1`` under pyfaidx but as
    ``contig_1`` in a BLAST database.  Lookup is also case-insensitive here but
    case-sensitive under pyfaidx.  Use :func:`check_ids` to detect IDs that will
    not resolve before extraction begins.

    Contig lengths are cached, since every lookup costs a subprocess spawn.
    """

    def __init__(
        self, path: Union[str, Path], descriptions: Optional[Dict[str, str]] = None
    ):
        self.path = str(path)
        self.descriptions = descriptions or {}
        self._length_cache: Dict[str, Optional[int]] = {}
        self._description_cache: Dict[str, str] = {}

    def _run(self, args: list) -> Optional[str]:
        """
        Run blastdbcmd with the given extra arguments.

        Parameters
        ----------
        args : list
            Arguments appended after ``blastdbcmd -db <path>``.

        Returns
        -------
        str or None
            Captured stdout, or None if the call failed or timed out.
        """
        cmd = ['blastdbcmd', '-db', self.path] + args
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=BLASTDBCMD_TIMEOUT,
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            logging.debug(f'blastdbcmd failed for {" ".join(args)}: {e.stderr.strip()}')
            return None
        except subprocess.TimeoutExpired:
            logging.error(
                f'blastdbcmd timed out after {BLASTDBCMD_TIMEOUT}s for {" ".join(args)}'
            )
            return None
        except FileNotFoundError:
            logging.error('blastdbcmd command not found. Is BLAST+ installed?')
            return None

    def _load_metadata(self, seqid: str) -> None:
        """
        Fetch and cache length and description in one call.

        Parameters
        ----------
        seqid : str
            Sequence identifier.

        Returns
        -------
        None
            Populates the length and description caches.

        Notes
        -----
        Each blastdbcmd invocation costs a subprocess spawn, so length and
        title are retrieved together.
        """
        # Note the literal tab: %t may itself contain spaces.
        out = self._run(['-entry', seqid, '-outfmt', '%l\t%t'])

        length: Optional[int] = None
        description = ''

        if out:
            # A database may hold several records under one accession; take the
            # first line only rather than silently concatenating.
            first = out.strip('\n').split('\n')[0]
            parts = first.split('\t', 1)
            try:
                length = int(parts[0])
            except (ValueError, IndexError):
                logging.error(f'Could not parse length of {seqid} from blastdbcmd')
                length = None
            if len(parts) > 1:
                description = parts[1].strip()

        self._length_cache[seqid] = length
        self._description_cache[seqid] = description

    def contig_length(self, seqid: str) -> Optional[int]:
        """
        Return the length of a sequence in the database.

        Parameters
        ----------
        seqid : str
            Sequence identifier.

        Returns
        -------
        int or None
            Length in bases, or None if ``seqid`` is not in the database.
        """
        if seqid not in self._length_cache:
            self._load_metadata(seqid)
        return self._length_cache[seqid]

    def sequence_description(self, seqid: str) -> str:
        """
        Return the description recorded for a sequence, or '' if there is none.

        Parameters
        ----------
        seqid : str
            Sequence identifier.

        Returns
        -------
        str
            The BLAST database title, which corresponds to the text following
            the sequence ID in the original FASTA header.
        """
        if seqid in self.descriptions:
            return self.descriptions[seqid]
        if seqid not in self._description_cache:
            self._load_metadata(seqid)
        return self._description_cache.get(seqid, '')

    def fetch_raw(self, seqid: str, start: int, end: int) -> Optional[str]:
        """
        Fetch a pre-clamped 1-based inclusive region on the plus strand.

        Parameters
        ----------
        seqid : str
            Sequence identifier.
        start, end : int
            1-based inclusive coordinates, already clamped to the contig.

        Returns
        -------
        str or None
            Raw sequence, or None if the call failed or returned bad FASTA.

        Notes
        -----
        Callers should use :func:`fetch_sequence` rather than calling this
        directly; it performs no clamping or validation.
        """
        out = self._run(['-entry', seqid, '-range', f'{start}-{end}'])
        if out is None:
            return None

        lines = out.strip().split('\n')
        if len(lines) < 2 or not lines[0].startswith('>'):
            logging.error(
                f'Invalid FASTA output from blastdbcmd for {seqid}:{start}-{end}'
            )
            return None

        # Filter by prefix rather than position so that multi-record output is
        # not silently spliced into the sequence.
        return ''.join(line for line in lines[1:] if not line.startswith('>'))

    def describe(self) -> str:
        """
        Describe this source for log messages.

        Returns
        -------
        str
            Human-readable identifier for this source.
        """
        return f'BLAST database {self.path}'


# Either backend satisfies the same duck-typed interface.
SequenceSource = Union[FastaSource, BlastDBSource]


def make_source(
    genome: Any = None,
    blastdb: Optional[Union[str, Path]] = None,
    descriptions: Optional[Dict[str, str]] = None,
) -> SequenceSource:
    """
    Build the appropriate sequence source from CLI arguments.

    Parameters
    ----------
    genome : pyfaidx.Fasta, optional
        Indexed genome. Used when ``blastdb`` is not given.
    blastdb : str or Path, optional
        Path to a BLAST database. Takes precedence over ``genome``.
    descriptions : dict, optional
        Mapping of sequence ID to FASTA header description.

    Returns
    -------
    FastaSource or BlastDBSource
        A source wrapping whichever backend was supplied.

    Raises
    ------
    ValueError
        If neither ``genome`` nor ``blastdb`` is provided.
    """
    if blastdb is not None:
        return BlastDBSource(blastdb, descriptions)
    if genome is not None:
        return FastaSource(genome, descriptions)
    raise ValueError('One of genome or blastdb must be provided')


def annotate(
    source: SequenceSource,
    seqid: str,
    coord_info: str,
    descriptions: Optional[Dict[str, str]] = None,
) -> str:
    """
    Append a sequence's description to a coordinate string for a FASTA header.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to look the description up in.
    seqid : str
        Sequence identifier.
    coord_info : str
        Coordinate portion of the description, built by the caller.
    descriptions : dict, optional
        Explicit description mapping, which takes precedence over the source.

    Returns
    -------
    str
        ``coord_info`` with the description appended, or ``coord_info``
        unchanged when there is no description.

    Notes
    -----
    Both backends can supply descriptions, so headers match regardless of which
    is in use. An empty description appends nothing, rather than a trailing
    space.
    """
    if descriptions is not None and seqid in descriptions:
        desc = descriptions[seqid]
    else:
        desc = source.sequence_description(seqid)

    desc = (desc or '').strip()
    return f'{coord_info} {desc}' if desc else coord_info


def contig_exists(source: SequenceSource, seqid: str) -> bool:
    """
    Report whether a sequence can be resolved in a source.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to check.
    seqid : str
        Sequence identifier.

    Returns
    -------
    bool
        True if ``seqid`` can be resolved.
    """
    return source.contig_length(seqid) is not None


def clamp_region(
    source: SequenceSource, seqid: str, start: int, end: int
) -> Optional[Tuple[int, int]]:
    """
    Clamp a 1-based inclusive region to the bounds of a contig.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to resolve the contig length against.
    seqid : str
        Sequence identifier.
    start, end : int
        1-based inclusive region, which may extend beyond either contig end.

    Returns
    -------
    tuple of (int, int) or None
        Clamped ``(start, end)``, or None if the contig is unknown or the
        region lies entirely outside it.

    Notes
    -----
    Callers that report coordinates in record IDs or descriptions should clamp
    first and report the clamped values, so that both backends label the same
    sequence identically.
    """
    length = source.contig_length(seqid)
    if length is None:
        logging.warning(f'Sequence {seqid} not found in {source.describe()}')
        return None

    start = int(start)
    end = int(end)

    if start > end:
        logging.debug(f'Empty region requested for {seqid}: {start}-{end}')
        return None

    if start > length or end < 1:
        logging.warning(
            f'Region {seqid}:{start}-{end} lies entirely outside the contig '
            f'(length {length}), skipping'
        )
        return None

    return max(1, start), min(length, end)


def fetch_sequence(
    source: SequenceSource, seqid: str, start: int, end: int
) -> Optional[str]:
    """
    Fetch a genomic region as uppercase plus-strand sequence.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to read from.
    seqid : str
        Sequence identifier.
    start, end : int
        1-based inclusive coordinates on the plus strand.  Coordinates
        extending past either contig end are clamped.

    Returns
    -------
    str or None
        Uppercase plus-strand sequence, or None if the region could not be
        retrieved in full.

    Notes
    -----
    Reverse-complementing is the caller's responsibility.  When flanking
    sequence is marked by lowercasing, apply the case marking to this
    plus-strand result *first* and reverse-complement afterwards; marking a
    sequence that has already been reverse-complemented puts the markers on the
    wrong end.
    """
    region = clamp_region(source, seqid, start, end)
    if region is None:
        return None

    clamped_start, clamped_end = region
    expected = clamped_end - clamped_start + 1

    seq = source.fetch_raw(seqid, clamped_start, clamped_end)
    if seq is None:
        logging.warning(
            f'Failed to extract {seqid}:{clamped_start}-{clamped_end} '
            f'from {source.describe()}'
        )
        return None

    # blastdbcmd returns the whole sequence, with a zero exit status, for some
    # out-of-range requests.  Clamping should prevent that, but verify rather
    # than trust: a length mismatch means we do not have the region we asked
    # for and returning it would silently corrupt downstream output.
    if len(seq) != expected:
        logging.warning(
            f'Extraction of {seqid}:{clamped_start}-{clamped_end} returned '
            f'{len(seq)}bp but {expected}bp was requested; skipping'
        )
        return None

    return seq.upper()


class PaddedRegion(NamedTuple):
    """
    A fixed-width region, padded where it extended beyond the contig.

    Attributes
    ----------
    seq : str
        Sequence of exactly the requested width.
    left_pad : int
        Number of pad characters prepended because the region started before
        position 1 of the contig.
    right_pad : int
        Number of pad characters appended because the region ran past the end
        of the contig.
    start, end : int
        The clamped 1-based inclusive coordinates of the *real* sequence, i.e.
        excluding the padding. Report these rather than the requested
        coordinates so record headers describe genuine positions.
    """

    seq: str
    left_pad: int
    right_pad: int
    start: int
    end: int

    @property
    def is_padded(self) -> bool:
        """
        Report whether any part of this region was synthesised.

        Returns
        -------
        bool
            True if either pad count is non-zero.
        """
        return bool(self.left_pad or self.right_pad)


def fetch_region_padded(
    source: SequenceSource,
    seqid: str,
    start: int,
    end: int,
    pad_char: str = 'N',
) -> Optional[PaddedRegion]:
    """
    Fetch a region of exactly the requested width, padding past contig ends.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to read from.
    seqid : str
        Sequence identifier.
    start, end : int
        1-based inclusive coordinates on the plus strand. Either may fall
        outside the contig.
    pad_char : str, default 'N'
        Character used for positions that do not exist in the contig.

    Returns
    -------
    PaddedRegion or None
        The padded region, or None if the contig is unknown or the requested
        window does not overlap it at all.

    Notes
    -----
    Callers that need a fixed width — flanks compared across records, TSDs of a
    declared length — should use this rather than :func:`fetch_sequence`, which
    silently returns a short sequence at a contig boundary.

    A window with **no** overlap returns None rather than an all-pad string: a
    region entirely off the contig carries no information, and emitting one
    would assert the existence of sequence that was never observed. Only
    partial overlap is padded.

    Padding is a statement that the genome does not extend far enough, not that
    the bases are unknown. Downstream consumers must not treat padded positions
    as evidence; :attr:`PaddedRegion.left_pad` and ``right_pad`` let them
    exclude those positions from comparisons and flag the record.
    """
    start = int(start)
    end = int(end)

    if start > end:
        logging.debug(f'Empty region requested for {seqid}: {start}-{end}')
        return None

    # clamp_region rejects windows with no overlap, which is what we want here.
    region = clamp_region(source, seqid, start, end)
    if region is None:
        return None

    clamped_start, clamped_end = region

    seq = fetch_sequence(source, seqid, clamped_start, clamped_end)
    if seq is None:
        return None

    left_pad = clamped_start - start
    right_pad = end - clamped_end

    if left_pad or right_pad:
        logging.warning(
            f'Region {seqid}:{start}-{end} extends beyond the contig; '
            f'padding with {left_pad}bp + {right_pad}bp of "{pad_char}" '
            f'around the available {clamped_start}-{clamped_end}'
        )

    return PaddedRegion(
        seq=(pad_char * left_pad) + seq + (pad_char * right_pad),
        left_pad=left_pad,
        right_pad=right_pad,
        start=clamped_start,
        end=clamped_end,
    )


def check_ids(source: SequenceSource, seqids: Any) -> list:
    """
    Report sequence IDs that cannot be resolved in ``source``.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to check against.
    seqids : iterable of str
        Sequence identifiers to verify.  Duplicates are checked once.

    Returns
    -------
    list of str
        Sorted IDs that could not be resolved.

    Notes
    -----
    Intended as a pre-flight check before bulk extraction.  The common cause of
    failure on the BLAST side is a database built without ``-parse_seqids``, or
    an accession that differs from the FASTA header token (``>lcl|contig_1``
    resolves as ``contig_1``).
    """
    missing = {sid for sid in set(seqids) if not contig_exists(source, sid)}

    if missing:
        logging.warning(
            f'{len(missing)} sequence ID(s) could not be resolved in '
            f'{source.describe()}: {", ".join(sorted(missing)[:5])}'
            + (' ...' if len(missing) > 5 else '')
        )
        if isinstance(source, BlastDBSource):
            logging.warning(
                'Was the BLAST database created with -parse_seqids? '
                'Accessions may also differ from FASTA header tokens.'
            )

    return sorted(missing)


def blastdbcmd_available() -> bool:
    """
    Report whether blastdbcmd is available.

    Returns
    -------
    bool
        True if ``blastdbcmd`` is on PATH.
    """
    return shutil.which('blastdbcmd') is not None
