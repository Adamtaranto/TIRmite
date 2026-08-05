#!/usr/bin/env python3
"""
TIRmite-validate: Validate reconstructed target sites.

This module validates reconstructed target sites by searching for
matching intact target sites in a genome database:
1. Run blastn of reconstructed target sites against a validation database
2. Filter for hits spanning the target site junction
3. Extract hit regions and re-align with mafft
4. Validate TSD/DR length from alignment gaps
"""

import argparse
import csv
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple, cast

from Bio import SeqIO  # type: ignore[import-not-found]
from Bio.Seq import Seq  # type: ignore[import-not-found]
from Bio.SeqRecord import SeqRecord  # type: ignore[import-not-found]

from tirmite._version import __version__  # type: ignore[import-not-found]
from tirmite.runners.mafft import align_in_memory, mafft_available
from tirmite.utils.extract import (
    BlastDBSource,
    blastdbcmd_available,
    fetch_sequence,
)
from tirmite.utils.logs import init_logging

logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Raised when validation cannot be performed as requested.

    Used for conditions that must not be mistaken for "nothing to validate":
    a missing or unparseable BLAST file, a missing external tool, or a run in
    which no site could be checked at all.
    """


def add_validate_parser(subparsers: Any) -> argparse.ArgumentParser:
    """
    Add validate subcommand parser.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        Subparser object to add validate command to.

    Returns
    -------
    argparse.ArgumentParser
        The configured validate subcommand parser.
    """
    parser = cast(
        argparse.ArgumentParser,
        subparsers.add_parser(
            'validate',
            help='Validate reconstructed target sites',
            description=(
                'Validate reconstructed target sites by searching for matching '
                'intact (empty) sites in a genome database.'
            ),
        ),
    )
    _configure_validate_parser(parser)
    return parser


def create_validate_parser() -> argparse.ArgumentParser:
    """
    Create standalone argument parser for validate command.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser for validate workflow options.
    """
    parser = argparse.ArgumentParser(
        prog='tirmite-validate',
        description='Validate reconstructed target sites against a genome database',
    )
    _configure_validate_parser(parser)
    return parser


def _configure_validate_parser(parser: argparse.ArgumentParser) -> None:
    """
    Configure parser with validate command arguments.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to configure.

    Returns
    -------
    None
        Modifies parser in place.
    """
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__),
    )

    parser.add_argument(
        '--loglevel',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help="Set logging level. Default: 'INFO'",
    )

    parser.add_argument(
        '--target-sites',
        type=str,
        required=True,
        dest='target_sites',
        help=(
            'Path to reconstructed target site FASTA file from tirmite pair '
            'with target site reconstruction enabled.'
        ),
    )

    parser.add_argument(
        '--blastdb',
        type=str,
        required=True,
        help='Path to BLAST database for validation searches.',
    )

    parser.add_argument(
        '--blast-results',
        type=str,
        default=None,
        dest='blast_results',
        help=(
            'Path to pre-computed BLAST results in format 6. If provided, '
            'skip running blastn and use these results instead. '
            'Expected columns: qseqid sseqid pident length mismatch gapopen '
            'qstart qend sstart send evalue bitscore qlen slen sstrand.'
        ),
    )

    parser.add_argument(
        '--min-coverage',
        type=float,
        default=0.95,
        dest='min_coverage',
        help=(
            'Minimum fraction of the query that must be covered by a hit '
            'to pass filtering. Default: 0.95'
        ),
    )

    parser.add_argument(
        '--evalue',
        type=float,
        default=1e-5,
        help='E-value threshold for blastn searches. Default: 1e-5',
    )

    parser.add_argument(
        '--tsd-length',
        type=int,
        default=0,
        dest='tsd_length',
        help=(
            'Default TSD/DR length to validate. '
            'Used when --tsd-length-map is not provided.'
        ),
    )

    parser.add_argument(
        '--tsd-length-map',
        type=str,
        default=None,
        dest='tsd_length_map',
        help=(
            'Path to tab-delimited file mapping model pairs to TSD lengths. '
            'Format: left_model<TAB>right_model<TAB>tsd_length.'
        ),
    )

    parser.add_argument(
        '--tsd-in-model',
        action='store_true',
        default=False,
        dest='tsd_in_model',
        help='TSD/DR feature is part of the termini model.',
    )

    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='Directory for output files. Default: current directory.',
    )

    parser.add_argument(
        '--prefix',
        type=str,
        default=None,
        help='Prefix for output file names.',
    )

    parser.add_argument(
        '--logfile',
        action='store_true',
        default=False,
        help='Write log messages to file in output directory.',
    )


def _run_validation_blastn(
    query_fasta: str,
    blastdb: str,
    output_file: str,
    evalue: float = 1e-5,
) -> None:
    """
    Run the blastn search specific to target-site validation.

    Parameters
    ----------
    query_fasta : str
        Path to query FASTA file.
    blastdb : str
        Path to BLAST database.
    output_file : str
        Path to output file for BLAST results.
    evalue : float, default 1e-5
        E-value threshold.

    Raises
    ------
    RuntimeError
        If blastn is not found or fails.

    Notes
    -----
    Deliberately does not delegate to
    :func:`tirmite.runners.runBlastn.run_blastn`, despite the overlap. That
    function always applies ``-word_size 4 -perc_identity 60``, which suit
    finding diverged transposon termini but are wrong here: validation looks
    for a near-identical empty target site, and a 60% identity floor combined
    with a 4-base word size would both admit spurious hits and run far slower
    than needed.

    This search also requires the extended 15-column ``-outfmt`` (adding
    ``qlen slen sstrand``), which :func:`parse_blast_results` depends on, and
    ``-dust no`` so that low-complexity target sites are not masked away.

    Renamed from ``run_blastn`` to remove the name collision with the runners
    module; it is private because nothing outside this workflow should use it.
    """
    if not shutil.which('blastn'):
        raise RuntimeError('blastn not found in PATH. Please install NCBI BLAST+.')

    cmd = [
        'blastn',
        '-query',
        query_fasta,
        '-db',
        blastdb,
        '-outfmt',
        '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand',
        '-evalue',
        str(evalue),
        '-out',
        output_file,
        '-dust',
        'no',
    ]

    logger.info(f'Running blastn: {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f'blastn failed: {result.stderr}')
        raise RuntimeError(f'blastn failed with return code {result.returncode}')


def parse_blast_results(
    blast_file: str,
) -> List[Dict[str, Any]]:
    """
    Parse BLAST format 6 results with extended columns.

    Parameters
    ----------
    blast_file : str
        Path to BLAST format 6 output file.

    Returns
    -------
    list of dict
        List of hit dictionaries with keys: qseqid, sseqid, pident, length,
        mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore,
        qlen, slen, sstrand.
    """
    columns = [
        'qseqid',
        'sseqid',
        'pident',
        'length',
        'mismatch',
        'gapopen',
        'qstart',
        'qend',
        'sstart',
        'send',
        'evalue',
        'bitscore',
        'qlen',
        'slen',
        'sstrand',
    ]

    hits: List[Dict[str, Any]] = []

    if not os.path.exists(blast_file):
        # An unreadable results file must not look like "no empty sites found".
        raise ValidationError(f'BLAST results file not found: {blast_file}')

    skipped_short = 0
    skipped_malformed = 0
    data_rows = 0

    with open(blast_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            # Comment and blank lines are legitimate in BLAST output.
            if not row or row[0].startswith('#'):
                continue
            data_rows += 1

            if len(row) < len(columns):
                skipped_short += 1
                continue

            hit: Dict[str, Any] = {}
            try:
                for i, col in enumerate(columns):
                    val = row[i]
                    if col in (
                        'qstart',
                        'qend',
                        'sstart',
                        'send',
                        'length',
                        'mismatch',
                        'gapopen',
                        'qlen',
                        'slen',
                    ):
                        hit[col] = int(val)
                    elif col in ('pident', 'evalue', 'bitscore'):
                        hit[col] = float(val)
                    else:
                        hit[col] = val
            except ValueError:
                # A non-numeric value in a numeric column: report rather than
                # crashing with a traceback from deep inside the parser.
                skipped_malformed += 1
                continue

            hits.append(hit)

    # A file whose every row was discarded is a format mismatch, not an empty
    # result. The overwhelmingly common cause is a standard 12-column
    # `-outfmt 6` file: TIRmite needs the 15-column form that adds
    # qlen, slen and sstrand. Silently returning [] made that indistinguishable
    # from "the genome contains no empty sites", and the run exited 0.
    if data_rows and not hits:
        raise ValidationError(
            f'None of the {data_rows} data row(s) in {blast_file} could be parsed. '
            f'{skipped_short} row(s) had fewer than {len(columns)} columns and '
            f'{skipped_malformed} had unparseable values. TIRmite requires the '
            "extended format: -outfmt '6 qseqid sseqid pident length mismatch "
            "gapopen qstart qend sstart send evalue bitscore qlen slen sstrand'."
        )

    if skipped_short or skipped_malformed:
        logger.warning(
            f'Skipped {skipped_short} short and {skipped_malformed} malformed '
            f'row(s) in {blast_file}'
        )

    logger.info(f'Parsed {len(hits)} hits from {blast_file}')
    return hits


def filter_junction_spanning(
    hits: List[Dict[str, Any]],
    min_coverage: float = 0.95,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Filter BLAST hits for junction-spanning hits with sufficient query coverage.

    Parameters
    ----------
    hits : list of dict
        Parsed BLAST hits.
    min_coverage : float, default 0.95
        Minimum fraction of query covered by the hit.

    Returns
    -------
    dict
        Dictionary mapping query ID to list of passing hits.
    """
    filtered: Dict[str, List[Dict[str, Any]]] = {}

    for hit in hits:
        qlen = hit['qlen']
        qstart = hit['qstart']
        qend = hit['qend']
        midpoint = qlen / 2.0

        # Check if hit spans the midpoint (junction)
        hit_qstart = min(qstart, qend)
        hit_qend = max(qstart, qend)
        if hit_qstart > midpoint or hit_qend < midpoint:
            continue

        # Check coverage
        alignment_len = hit_qend - hit_qstart + 1
        coverage = alignment_len / qlen
        if coverage < min_coverage:
            continue

        qid = hit['qseqid']
        if qid not in filtered:
            filtered[qid] = []
        filtered[qid].append(hit)

    total_passing = sum(len(v) for v in filtered.values())
    logger.info(
        f'Filtered to {total_passing} junction-spanning hits '
        f'across {len(filtered)} queries'
    )
    return filtered


def extract_hit_sequence(
    source: BlastDBSource,
    subject_id: str,
    start: int,
    end: int,
    strand: str,
) -> Optional[str]:
    """
    Extract a sequence region from a BLAST database.

    Parameters
    ----------
    source : BlastDBSource
        Open sequence source for the BLAST database. Takes a source rather
        than a path so its caches survive across hits; see Notes.
    subject_id : str
        Subject sequence identifier.
    start : int
        1-based inclusive start coordinate, on the plus strand.
    end : int
        1-based inclusive end coordinate, on the plus strand.
    strand : str
        Strand: 'plus' or 'minus'. The result is reverse complemented
        for 'minus'.

    Returns
    -------
    str or None
        Extracted sequence, or None if extraction failed.

    Notes
    -----
    Delegates to :func:`tirmite.utils.extract.fetch_sequence`, which applies the
    same clamping and validation used for FASTA-indexed genomes.

    This used to take a database *path* and construct a fresh
    ``BlastDBSource`` on every call. That class caches contig lengths
    precisely because each lookup costs a ``blastdbcmd`` subprocess, and
    rebuilding it per hit threw the cache away every time -- so each hit paid
    two subprocess spawns (one to resolve the contig length for clamping, one
    to fetch the sequence) instead of amortising the first across every hit on
    the same contig. ``blastdbcmd_available()`` was re-checked per hit too;
    that is now pre-flighted once by :func:`_check_required_tools`.
    """
    seq = fetch_sequence(source, subject_id, start, end)
    if seq is None:
        return None

    if strand == 'minus':
        seq = str(Seq(seq).reverse_complement())

    return seq if seq else None


def parse_target_site_metadata(description: str) -> Dict[str, str]:
    """
    Parse the ``key=value`` metadata written into a target-site FASTA header.

    Parameters
    ----------
    description : str
        The FASTA description line produced by ``tirmite pair``.

    Returns
    -------
    dict of str to str
        Every ``key=value`` token found. Keys with no ``=`` are ignored.

    Notes
    -----
    ``writeTargetSites`` records ``flank_len``, ``tsd_len``, ``tsd_in_model``,
    the model names, the contig and the element orientation. Reading them makes
    validation self-consistent with the run that produced the sites, instead of
    relying on the user to re-supply matching ``--tsd-length`` and
    ``--tsd-in-model`` values on the command line.
    """
    metadata: Dict[str, str] = {}
    for token in description.split():
        if '=' in token:
            key, _, value = token.partition('=')
            metadata[key] = value
    return metadata


def compute_junction_position(
    flank_len: Optional[int], tsd_len: int, tsd_in_model: bool, query_len: int
) -> int:
    """
    Locate the element/genome junction within a reconstructed target site.

    Parameters
    ----------
    flank_len : int or None
        Flank width used when the site was reconstructed. When None the
        junction falls back to the query midpoint.
    tsd_len : int
        Declared TSD length.
    tsd_in_model : bool
        Whether the TSD lies inside the terminus model.
    query_len : int
        Length of the reconstructed target site.

    Returns
    -------
    int
        1-based position in the unaligned query around which to look for gaps.

    Notes
    -----
    The junction is not the midpoint of the query, which is what this code
    assumed. The two reconstruction modes build the site differently:

    - ``tsd_in_model=False``: ``left_flank + right_flank[tsd_len:]``, so the
      site is ``2 * flank_len - tsd_len`` long and the boundary sits at
      ``flank_len``. The midpoint is ``flank_len - tsd_len / 2``, i.e. off by
      half the TSD length -- 10 bp out for a 20 bp direct repeat.
    - ``tsd_in_model=True``: ``left_flank + tsd + right_flank``, so the TSD is
      centred at ``flank_len + tsd_len / 2``, which does coincide with the
      midpoint.

    Unequal flanks (a truncated flank at a contig edge) shift the junction
    further still, which is why the recorded ``flank_len`` is preferred over
    any calculation from ``query_len``.
    """
    if flank_len is None:
        return query_len // 2

    if tsd_in_model:
        # Centre of the TSD block that sits between the two flanks.
        return min(query_len, flank_len + max(tsd_len, 1) // 2)

    # Boundary between the left flank and the trimmed right flank.
    return min(query_len, flank_len)


def check_tsd_gaps(
    query_aligned: str,
    target_aligned: str,
    tsd_length: int,
    query_len: int,
    junction_pos: Optional[int] = None,
    tsd_in_model: bool = False,
) -> Tuple[Optional[int], str]:
    """
    Check for gaps around the alignment midpoint indicating TSD length errors.

    Parameters
    ----------
    query_aligned : str
        Aligned query sequence (with gaps).
    target_aligned : str
        Aligned target sequence (with gaps).
    tsd_length : int
        User-specified TSD length.
    query_len : int
        Original unaligned query length.
    junction_pos : int, optional
        1-based position of the element/genome junction in the unaligned
        query. Defaults to the query midpoint, which is only correct when the
        TSD lies inside the terminus model.
    tsd_in_model : bool, default False
        Whether the TSD lies inside the terminus model. This inverts the sign
        of the reported error; see Notes.

    Returns
    -------
    predicted_error : int or None
        Predicted error in TSD length. Positive means the declared TSD is too
        long, negative too short, and ``0`` that it agrees with the data.
        **None means the comparison was inconclusive** and carries no evidence
        either way.
    message : str
        Human-readable message about the validation result.

    Notes
    -----
    The None return is load-bearing. This previously returned ``0`` both when
    the alignment confirmed the declared length and when gaps on *both*
    sequences made the comparison meaningless, so an inconclusive result was
    averaged in as agreement and reported as "TSD length appears consistent".

    That is the same failure mode fixed elsewhere in 1.5.0, where
    :func:`tirmite.core.tsd.compare_tsds` was given an ``Optional[int]`` return
    for exactly this reason: an unverifiable TSD must never be indistinguishable
    from a verified one.
    """
    # Locate the junction in aligned coordinates.
    query_pos = 0
    midpoint_aligned = 0
    target_midpoint = query_len // 2 if junction_pos is None else junction_pos

    for i, c in enumerate(query_aligned):
        if c != '-':
            query_pos += 1
        if query_pos >= target_midpoint:
            midpoint_aligned = i
            break

    # Check a window around the midpoint for gaps
    # Padding ensures we look beyond the TSD region itself
    WINDOW_PADDING = 5
    MIN_WINDOW_SIZE = 10
    window = max(tsd_length + WINDOW_PADDING, MIN_WINDOW_SIZE)
    start = max(0, midpoint_aligned - window)
    end = min(len(query_aligned), midpoint_aligned + window + 1)

    query_gaps = query_aligned[start:end].count('-')
    target_gaps = target_aligned[start:end].count('-')

    # The two reconstruction modes respond to an over-long TSD in opposite
    # directions, so the sign of the measured difference has to be flipped for
    # one of them:
    #
    #   tsd_in_model=False: the site is left_flank + right_flank[tsd_len:], so
    #     over-declaring TRIMS real flank. The query is shorter than the true
    #     empty site, MAFFT gaps the QUERY, and gaps in the query mean "too
    #     long".
    #   tsd_in_model=True: the site is left_flank + tsd + right_flank, so
    #     over-declaring INSERTS extra element bases. The query is longer,
    #     MAFFT gaps the TARGET, and gaps in the query now mean "too short".
    #
    # Without this, an in-model TSD declared too long was reported as too
    # short, pushing the user to lengthen it further.
    sign = -1 if tsd_in_model else 1

    if query_gaps > 0 and target_gaps == 0:
        error = sign * query_gaps
        direction = 'too long' if error > 0 else 'too short'
        return error, (
            f'Query has {query_gaps} gap(s) near the junction: '
            f'TSD length may be {query_gaps}bp {direction}'
        )
    elif target_gaps > 0 and query_gaps == 0:
        error = -sign * target_gaps
        direction = 'too long' if error > 0 else 'too short'
        return error, (
            f'Target has {target_gaps} gap(s) near the junction: '
            f'TSD length may be {target_gaps}bp {direction}'
        )
    elif query_gaps == 0 and target_gaps == 0:
        return 0, 'TSD length appears consistent with validation data'
    else:
        # Gaps on both sides carry no directional information. Returning None
        # rather than 0 keeps this out of the average entirely.
        return None, (
            f'Both query ({query_gaps}) and target ({target_gaps}) have gaps '
            f'near the junction; validation inconclusive'
        )


def _check_required_tools(args: argparse.Namespace) -> None:
    """
    Verify the external tools this run will need are installed.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments, used to decide whether blastn is needed.

    Raises
    ------
    ValidationError
        If a required executable is not on PATH.

    Notes
    -----
    Checked up front rather than per hit. Both downstream helpers report a
    missing tool by returning None -- ``align_in_memory`` per alignment and
    ``extract_hit_sequence`` per hit -- and neither was aggregated, so a run
    without MAFFT or without blastdbcmd wrote a full summary reporting
    ``predicted_tsd_error 0.0`` for every site and exited 0.
    """
    missing = []

    if not blastdbcmd_available():
        missing.append('blastdbcmd (NCBI BLAST+) - needed to extract hit regions')

    if not mafft_available():
        missing.append('mafft - needed to align hits against the target site')

    # blastn is only needed when we are running the search ourselves.
    if not args.blast_results and not shutil.which('blastn'):
        missing.append('blastn (NCBI BLAST+) - needed to search for empty sites')

    if missing:
        raise ValidationError(
            'Required tool(s) not found on PATH:\n  ' + '\n  '.join(missing)
        )


def main(args: Optional[argparse.Namespace] = None) -> int:
    """
    Main entry point for tirmite-validate.

    Parameters
    ----------
    args : argparse.Namespace, optional
        Parsed command-line arguments. If None, parses from sys.argv.

    Returns
    -------
    int
        Exit code (0 for success, 1 for error).
    """
    if args is None:
        parser = create_validate_parser()
        args = parser.parse_args()

    assert args is not None

    try:
        # Set up output directory
        if args.outdir:
            outdir = Path(args.outdir)
        else:
            outdir = Path.cwd()
        outdir.mkdir(parents=True, exist_ok=True)

        # Set up logging
        logfile_path = None
        if args.logfile:
            prefix_str = f'{args.prefix}_' if args.prefix else ''
            logfile_path = outdir / f'{prefix_str}tirmite_validate.log'

        init_logging(loglevel=args.loglevel, logfile=logfile_path)

        logger.info(f'TIRmite-validate version {__version__}')
        logger.info(f'Output directory: {outdir}')

        # Pre-flight the external tools. Without these checks a missing tool
        # produced a complete report full of 0.0 errors and exit 0: MAFFT
        # failures return None per alignment, and blastdbcmd failures return
        # None per hit, and neither was ever aggregated into a run-level
        # failure.
        _check_required_tools(args)

        # Validate inputs
        if not Path(args.target_sites).exists():
            logger.error(f'Target sites file not found: {args.target_sites}')
            sys.exit(1)

        # Load query sequences
        logger.info(f'Loading target sites from {args.target_sites}')
        queries = list(SeqIO.parse(args.target_sites, 'fasta'))
        if not queries:
            logger.error('No sequences found in target sites file')
            sys.exit(1)
        logger.info(f'Loaded {len(queries)} target site queries')

        # Load TSD length map if provided
        tsd_length_map: Dict[str, int] = {}
        if args.tsd_length_map:
            from tirmite.core.tsd import load_tsd_length_map

            tsd_length_map = load_tsd_length_map(args.tsd_length_map)

        # Run or load BLAST results
        with tempfile.TemporaryDirectory() as tmpdir:
            if args.blast_results:
                logger.info(
                    f'Loading pre-computed BLAST results from {args.blast_results}'
                )
                blast_file = args.blast_results
            else:
                blast_file = os.path.join(tmpdir, 'blast_results.txt')
                _run_validation_blastn(
                    query_fasta=args.target_sites,
                    blastdb=args.blastdb,
                    output_file=blast_file,
                    evalue=args.evalue,
                )

            # Parse and filter results
            all_hits = parse_blast_results(blast_file)
            filtered_hits = filter_junction_spanning(
                all_hits, min_coverage=args.min_coverage
            )

            # One sequence source for the whole run. BlastDBSource caches
            # contig lengths and descriptions, and each cache miss costs a
            # blastdbcmd subprocess, so this must outlive the per-hit loop.
            blast_source = BlastDBSource(args.blastdb)

            # Process each query
            prefix_str = f'{args.prefix}_' if args.prefix else ''
            summary_rows: List[Dict[str, Any]] = []
            all_alignments: Dict[str, List[SeqRecord]] = {}
            # How many queries yielded at least one usable comparison. A run
            # where this stays zero validated nothing and must not exit 0.
            validated_queries = 0

            for query in queries:
                qid = query.id
                query_len = len(query.seq)
                logger.info(f'Processing query: {qid} (len={query_len})')

                # `tirmite pair` records everything needed to validate this
                # site in the FASTA header. Prefer it over the CLI flags, which
                # the user has to re-supply consistently by hand.
                metadata = parse_target_site_metadata(query.description)
                left_model = metadata.get('left_model', '')
                right_model = metadata.get('right_model', '')

                # Resolve TSD length: header first, then the map, then the flag.
                query_tsd_length = args.tsd_length
                source = '--tsd-length'

                if tsd_length_map and left_model and right_model:
                    # Try both key orders. `pair` looks the map up with
                    # coordinate-ordered models but writes ROLE-ordered ones to
                    # the header, so for a reverse-inserted element the two
                    # differ and the direct lookup misses.
                    for key in (
                        f'{left_model}\t{right_model}',
                        f'{right_model}\t{left_model}',
                    ):
                        if key in tsd_length_map:
                            query_tsd_length = tsd_length_map[key]
                            source = '--tsd-length-map'
                            break
                    else:
                        logger.warning(
                            f'{qid}: no entry for {left_model}/{right_model} in '
                            'the TSD length map (tried both orders)'
                        )

                header_tsd_len = metadata.get('tsd_len')
                if header_tsd_len is not None:
                    try:
                        query_tsd_length = int(header_tsd_len)
                        source = 'FASTA header'
                    except ValueError:
                        logger.warning(
                            f'{qid}: unparseable tsd_len={header_tsd_len!r} in header'
                        )

                # tsd_in_model changes the SIGN of the reported error, so
                # getting it from the header rather than the flag matters.
                header_in_model = metadata.get('tsd_in_model')
                if header_in_model is not None:
                    query_tsd_in_model = header_in_model.lower() == 'true'
                    if query_tsd_in_model != args.tsd_in_model:
                        logger.warning(
                            f'{qid}: --tsd-in-model={args.tsd_in_model} disagrees '
                            f'with the header (tsd_in_model={query_tsd_in_model}); '
                            'using the header, which records how the site was '
                            'actually reconstructed'
                        )
                else:
                    query_tsd_in_model = args.tsd_in_model

                header_flank_len = metadata.get('flank_len')
                query_flank_len: Optional[int] = None
                if header_flank_len is not None:
                    try:
                        query_flank_len = int(header_flank_len)
                    except ValueError:
                        logger.warning(
                            f'{qid}: unparseable flank_len={header_flank_len!r}'
                        )

                junction_pos = compute_junction_position(
                    query_flank_len,
                    query_tsd_length,
                    query_tsd_in_model,
                    query_len,
                )
                logger.debug(
                    f'{qid}: tsd_length={query_tsd_length} (from {source}), '
                    f'tsd_in_model={query_tsd_in_model}, junction={junction_pos}'
                )

                query_hits = filtered_hits.get(qid, [])
                logger.info(f'  {len(query_hits)} hits passing filters')

                if not query_hits:
                    summary_rows.append(
                        {
                            'query_id': qid,
                            'left_model': left_model,
                            'right_model': right_model,
                            'query_length': query_len,
                            'tsd_length': query_tsd_length,
                            'num_filtered_hits': 0,
                            'num_compared_hits': 0,
                            'predicted_tsd_error': 'NA',
                            'validation_message': 'No valid empty site hits found',
                        }
                    )
                    continue

                # Extract hit sequences and re-align
                alignment_seqs = [query]
                tsd_errors: List[int] = []
                tsd_messages: List[str] = []
                inconclusive = 0

                for hit in query_hits:
                    sstart = hit['sstart']
                    send = hit['send']

                    # Handle strand: when sstart > send, coordinates indicate
                    # minus strand regardless of sstrand field
                    sstrand = hit.get('sstrand', 'plus')
                    if sstart > send:
                        strand = 'minus'
                        sstart, send = send, sstart
                    else:
                        strand = sstrand

                    seq_str = extract_hit_sequence(
                        blast_source, hit['sseqid'], sstart, send, strand
                    )
                    if seq_str is None:
                        continue

                    hit_id = f'{hit["sseqid"]}_{sstart}_{send}_{strand}'
                    hit_record = SeqRecord(
                        Seq(seq_str),
                        id=hit_id,
                        name=hit_id,
                        description=f'evalue={hit["evalue"]} pident={hit["pident"]}',
                    )
                    alignment_seqs.append(hit_record)

                if len(alignment_seqs) > 1:
                    # Run MAFFT alignment
                    aligned = align_in_memory(alignment_seqs, tmpdir)

                    if aligned:
                        all_alignments[qid] = aligned
                        query_aligned = str(aligned[0].seq)

                        for aln_rec in aligned[1:]:
                            target_aligned = str(aln_rec.seq)
                            error, msg = check_tsd_gaps(
                                query_aligned,
                                target_aligned,
                                query_tsd_length,
                                query_len,
                                junction_pos=junction_pos,
                                tsd_in_model=query_tsd_in_model,
                            )
                            # None means inconclusive; count it separately so
                            # it cannot be averaged in as agreement.
                            if error is None:
                                inconclusive += 1
                            else:
                                tsd_errors.append(error)
                            tsd_messages.append(msg)

                # Determine consensus TSD error
                if tsd_errors:
                    avg_error = sum(tsd_errors) / len(tsd_errors)
                    if abs(avg_error) < 0.5:
                        consensus_msg = 'TSD length appears consistent'
                    elif avg_error > 0:
                        consensus_msg = (
                            f'TSD may be ~{abs(round(avg_error))}bp too long'
                        )
                    else:
                        consensus_msg = (
                            f'TSD may be ~{abs(round(avg_error))}bp too short'
                        )

                    if inconclusive:
                        consensus_msg += (
                            f' (from {len(tsd_errors)} comparison(s); '
                            f'{inconclusive} inconclusive)'
                        )

                    if abs(avg_error) >= 1:
                        logger.warning(
                            f'TSD length validation warning for {qid}: {consensus_msg}'
                        )
                    reported_error: Any = round(avg_error, 1)
                    validated_queries += 1
                elif inconclusive:
                    # Every comparison was inconclusive. Reporting 0.0 here is
                    # what made an unverifiable result look like a confirmed one.
                    reported_error = 'NA'
                    consensus_msg = (
                        f'All {inconclusive} comparison(s) inconclusive: '
                        'gaps on both sequences near the junction'
                    )
                else:
                    reported_error = 'NA'
                    consensus_msg = 'No alignments available for validation'

                summary_rows.append(
                    {
                        'query_id': qid,
                        'left_model': left_model,
                        'right_model': right_model,
                        'query_length': query_len,
                        'tsd_length': query_tsd_length,
                        # Hits that survived BLAST filtering, and hits that were
                        # actually extracted and compared. These differ whenever
                        # extraction fails, and only the second says how much
                        # evidence the verdict rests on.
                        'num_filtered_hits': len(query_hits),
                        'num_compared_hits': len(tsd_errors) + inconclusive,
                        'predicted_tsd_error': reported_error,
                        'validation_message': consensus_msg,
                    }
                )

            # Write summary report
            summary_file = outdir / f'{prefix_str}validation_summary.tsv'
            with open(summary_file, 'w', newline='') as f:
                if summary_rows:
                    writer = csv.DictWriter(
                        f, fieldnames=summary_rows[0].keys(), delimiter='\t'
                    )
                    writer.writeheader()
                    writer.writerows(summary_rows)
            logger.info(f'Validation summary written to {summary_file}')

            # Write aligned FASTA files
            for qid, aligned_records in all_alignments.items():
                safe_qid = qid.replace('/', '_').replace('\\', '_')
                aln_file = outdir / f'{prefix_str}{safe_qid}_alignment.fasta'
                with open(aln_file, 'w') as handle:
                    SeqIO.write(aligned_records, handle, 'fasta')
                logger.info(f'Alignment written to {aln_file}')

        # A run that could not validate a single site is a failure, however
        # many rows the summary has. Previously this reported success and
        # exited 0 whether it validated 500 sites or none of them.
        if validated_queries == 0:
            logger.error(
                f'No target site could be validated ({len(summary_rows)} '
                'queries processed). See the summary for per-query reasons.'
            )
            return 1

        logger.info(
            f'TIRmite-validate completed: {validated_queries} of '
            f'{len(summary_rows)} target site(s) validated'
        )

    except ValidationError as e:
        logger.error(f'Validation cannot proceed: {e}')
        return 1
    except KeyboardInterrupt:
        logger.info('Analysis interrupted by user')
        sys.exit(130)
    except Exception as e:
        logger.error(f'Unexpected error: {e}')
        logger.exception('Full traceback:')
        sys.exit(1)

    return 0


if __name__ == '__main__':
    main()
