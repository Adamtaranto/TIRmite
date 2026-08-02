"""External flank extraction around element termini.

Computes the genomic window immediately outside each terminus and extracts it,
optionally padding past a contig boundary so the window is always the
requested width.
"""

import logging
import os
from typing import Any, Dict, List, NamedTuple, Optional, Set, Tuple

from Bio import Seq, SeqIO  # type: ignore[import-not-found]
from Bio.SeqRecord import SeqRecord  # type: ignore[import-not-found]
import pandas as pd  # type: ignore[import-untyped]

from tirmite.core.termini import _model_deficit, _pair_roles, flipTIRs, resolve_terminus
from tirmite.utils.extract import (
    SequenceSource,
    annotate,
    clamp_region,
    fetch_region_padded,
    fetch_sequence,
    make_source,
)
from tirmite.utils.utils import cleanID


def compute_flank_coordinates(
    hit_start: int,
    hit_end: int,
    strand: str,
    is_left_terminus: bool,
    hmm_start: int,
    hmm_end: int,
    model_len: int,
    flank_len: int,
) -> Tuple[int, int, int]:
    """
    Compute genomic coordinates for the external flanking region of a terminus hit.

    The "external" end of a terminus hit is the side that faces away from the TE
    body. For the left terminus this is the side at lower genomic coordinates; for
    the right terminus it is the side at higher genomic coordinates.

    When a hit does not cover position 1 of the model (hmmStart > 1) the external
    boundary must be shifted by the number of uncovered model positions so that the
    reported flank begins at the correct genomic position.

    Parameters
    ----------
    hit_start : int
        1-based start coordinate of the hit in genomic coordinates (always < hit_end).
    hit_end : int
        1-based end coordinate of the hit in genomic coordinates.
    strand : str
        Strand of the hit: '+' or '-'.
    is_left_terminus : bool
        True if the hit represents the left (5') terminus of the element.
    hmm_start : int
        Alignment start position on the HMM model (1-based).
        For + strand hits this aligns to hit_start; for - strand hits it aligns to hit_end.
    hmm_end : int
        Alignment end position on the HMM model (1-based).
        For + strand hits this aligns to hit_end; for - strand hits it aligns to hit_start.
    model_len : int
        Total length of the HMM model in positions.
    flank_len : int
        Number of bases to extract in the flanking region.

    Returns
    -------
    flank_start : int
        1-based start coordinate of the flank region.
    flank_end : int
        1-based end coordinate of the flank region (inclusive).
    offset : int
        Number of model positions between the hit alignment and the external
        end of the model (0 means the alignment reaches the model end).

    Notes
    -----
    Coordinate system:
      - For + strand: hmmStart aligns to hit_start, hmmEnd aligns to hit_end.
      - For - strand: hmmStart aligns to hit_end (higher coord), hmmEnd aligns to hit_start.

    Left terminus external boundary:
      - + strand: external_pos = hit_start - (hmm_start - 1)
      - - strand: external_pos = hit_start - (model_len - hmm_end)

    Right terminus external boundary:
      - + strand: external_pos = hit_end + (model_len - hmm_end)
      - - strand: external_pos = hit_end + (hmm_start - 1)
    """
    if is_left_terminus:
        # External end faces LEFT (lower genomic coordinates)
        if strand == '+':
            offset = _model_deficit(
                hmm_start - 1, 'hmmStart', hmm_start, hmm_end, model_len
            )
            external_pos = hit_start - offset
        else:  # '-'
            # hmmStart aligns to hit_end (higher coord); hmmEnd aligns to hit_start
            offset = _model_deficit(
                model_len - hmm_end, 'hmmEnd', hmm_start, hmm_end, model_len
            )
            external_pos = hit_start - offset
        flank_start = external_pos - flank_len
        flank_end = external_pos - 1
    else:
        # External end faces RIGHT (higher genomic coordinates)
        if strand == '+':
            offset = _model_deficit(
                model_len - hmm_end, 'hmmEnd', hmm_start, hmm_end, model_len
            )
            external_pos = hit_end + offset
        else:  # '-'
            # hmmStart aligns to hit_end; external end is at hit_end side
            offset = _model_deficit(
                hmm_start - 1, 'hmmStart', hmm_start, hmm_end, model_len
            )
            external_pos = hit_end + offset
        flank_start = external_pos + 1
        flank_end = external_pos + flank_len

    return flank_start, flank_end, offset


def compute_inner_tsd_coordinates(
    hit_start: int,
    hit_end: int,
    strand: str,
    is_left_terminus: bool,
    hmm_start: int,
    hmm_end: int,
    model_len: int,
    tsd_length: int,
) -> Tuple[int, int]:
    """
    Compute genomic coordinates of the TSD at the inner boundary of a terminus hit.

    The "inner boundary" of a terminus hit is the side that faces the element body.
    For a left terminus this is the RIGHT (higher coordinate) end; for a right
    terminus it is the LEFT (lower coordinate) end.  When the TSD is modelled as
    part of the terminus HMM it occupies the innermost ``tsd_length`` model
    positions of the terminus model.

    Parameters
    ----------
    hit_start : int
        1-based start coordinate of the hit (always ≤ hit_end).
    hit_end : int
        1-based end coordinate of the hit (always ≥ hit_start).
    strand : str
        Strand of the hit: '+' or '-'.
    is_left_terminus : bool
        True if the hit is the left (5') terminus of the element.
    hmm_start : int
        1-based alignment start on the HMM model.
        For + strand hits this aligns with ``hit_start``; for - strand
        hits it aligns with ``hit_end``.
    hmm_end : int
        1-based alignment end on the HMM model.
        For + strand hits this aligns with ``hit_end``; for - strand
        hits it aligns with ``hit_start``.
    model_len : int
        Total number of positions in the HMM model.
    tsd_length : int
        Number of bases in the TSD feature.

    Returns
    -------
    tsd_start : int
        1-based genomic start coordinate of the TSD (always ≤ tsd_end).
    tsd_end : int
        1-based genomic end coordinate of the TSD.

    Notes
    -----
    Coordinate conventions (same as :func:`compute_flank_coordinates`):

    For + strand hits, model positions increase in the same direction as
    genomic coordinates (hmmStart aligns to hit_start).

    For - strand hits, model positions increase in the *opposite* direction
    (hmmStart aligns to hit_end; hmmEnd aligns to hit_start).

    Left terminus inner boundary
        - ``+`` strand: model position ``model_len`` is at the RIGHT (higher)
          end of the hit.  ``inner_pos = hit_end + (model_len - hmm_end)``.
          TSD occupies ``[inner_pos - tsd_length + 1, inner_pos]``.
        - ``-`` strand: model position 1 is at the RIGHT (higher) end of the
          hit.  ``inner_pos = hit_end + (hmm_start - 1)``.
          TSD occupies ``[inner_pos - tsd_length + 1, inner_pos]``.

    Right terminus inner boundary
        - ``+`` strand: model position 1 is at the LEFT (lower) end of the
          hit.  ``inner_pos = hit_start - (hmm_start - 1)``.
          TSD occupies ``[inner_pos, inner_pos + tsd_length - 1]``.
        - ``-`` strand: model position ``model_len`` is at the LEFT (lower)
          end of the hit.  ``inner_pos = hit_start - (model_len - hmm_end)``.
          TSD occupies ``[inner_pos, inner_pos + tsd_length - 1]``.
    """
    # Deficits are clamped to >= 0; see _model_deficit for why a negative value
    # is possible and what it means.
    end_deficit = _model_deficit(
        model_len - hmm_end, 'hmmEnd', hmm_start, hmm_end, model_len
    )
    start_deficit = _model_deficit(
        hmm_start - 1, 'hmmStart', hmm_start, hmm_end, model_len
    )

    if is_left_terminus:
        # Inner end faces RIGHT (higher genomic coords)
        if strand == '+':
            # Model pos model_len aligns to right of hit
            inner_pos = hit_end + end_deficit
        else:  # '-'
            # Model pos 1 (hmmStart) aligns to hit_end (rightmost genomic coord)
            inner_pos = hit_end + start_deficit
        tsd_start = inner_pos - tsd_length + 1
        tsd_end = inner_pos
    else:
        # Right terminus: inner end faces LEFT (lower genomic coords)
        if strand == '+':
            # Model pos 1 (hmmStart) aligns to hit_start (leftmost genomic coord)
            inner_pos = hit_start - start_deficit
        else:  # '-'
            # Model pos model_len aligns to hit_start (leftmost genomic coord)
            inner_pos = hit_start - end_deficit
        tsd_start = inner_pos
        tsd_end = inner_pos + tsd_length - 1

    return tsd_start, tsd_end


class FlankResult(NamedTuple):
    """
    External flanking sequence for one terminus hit.

    Attributes
    ----------
    seq : str
        Flank sequence on the plus strand, always ``flank_len`` bases when
        ``pad`` was requested.
    start, end : int
        1-based inclusive coordinates of the real (non-padded) sequence.
    offset : int
        Uncovered model positions between the hit alignment and the projected
        external edge of the element.
    left_pad, right_pad : int
        Bases synthesised because the flank ran past a contig boundary.
    """

    seq: str
    start: int
    end: int
    offset: int
    left_pad: int
    right_pad: int

    @property
    def is_padded(self) -> bool:
        """
        Report whether any part of this flank was synthesised.

        Returns
        -------
        bool
            True if either pad count is non-zero.
        """
        return bool(self.left_pad or self.right_pad)


def extract_terminus_flank(
    source: SequenceSource,
    hit: Any,
    is_left: bool,
    model_len: Optional[int],
    hmm_start: Optional[int],
    hmm_end: Optional[int],
    flank_len: int,
    flank_max_offset: Optional[int] = None,
    pad: bool = True,
) -> Optional[FlankResult]:
    """
    Extract the external flanking region of a single terminus hit.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to read from.
    hit : namedtuple
        Hit record with target, hitStart, hitEnd, strand, idx.
    is_left : bool
        True if the hit is a left terminus; False for a right terminus.
    model_len : int or None
        Length of the terminus model. None skips the hit with a warning.
    hmm_start, hmm_end : int or None
        Alignment coordinates on the model. None skips the hit.
    flank_len : int
        Number of flanking bases to extract.
    flank_max_offset : int, optional
        Reject hits whose model-coverage deficit exceeds this.
    pad : bool, default True
        Pad with N when the flank runs past a contig boundary, so every flank
        is ``flank_len`` bases. When False the flank is truncated instead.

    Returns
    -------
    FlankResult or None
        None when the flank cannot be extracted: missing model length or HMM
        coordinates, offset above ``flank_max_offset``, or a region that does
        not overlap the contig at all.

    Notes
    -----
    Flanks are always taken on the plus strand for genomic orientation; the
    hit's own strand only affects which side of the hit the flank is on, via
    :func:`compute_flank_coordinates`.

    This is the single implementation used by both ``writeFlanks`` and
    ``writeTargetSites``. They previously carried near-identical private copies
    that had already drifted in their logging.
    """
    if model_len is None:
        logging.warning(
            f'Model length not found for {hit.model}, skipping flank for hit {hit.idx}'
        )
        return None

    if hmm_start is None or hmm_end is None:
        logging.warning(
            f'HMM coordinates unavailable for hit {hit.idx}, skipping flank'
        )
        return None

    flank_start, flank_end, offset = compute_flank_coordinates(
        hit_start=int(hit.hitStart),
        hit_end=int(hit.hitEnd),
        strand=hit.strand,
        is_left_terminus=is_left,
        hmm_start=hmm_start,
        hmm_end=hmm_end,
        model_len=model_len,
        flank_len=flank_len,
    )

    if flank_max_offset is not None and offset > flank_max_offset:
        logging.debug(
            f'Skipping flank for hit {hit.idx}: offset {offset} > max {flank_max_offset}'
        )
        return None

    # A flank entirely before the contig start is reported distinctly from one
    # that merely overhangs, since it means the element sits at the very edge of
    # the assembly and no flanking context exists at all.
    if flank_end < 1:
        logging.warning(
            f'Flank for hit {hit.idx} on {hit.target} falls entirely before '
            f'contig start (computed coords {flank_start}–{flank_end}), skipping'
        )
        return None

    if pad:
        region = fetch_region_padded(source, hit.target, flank_start, flank_end)
        if region is None:
            logging.debug(f'Empty flank region for hit {hit.idx}, skipping')
            return None
        return FlankResult(
            seq=region.seq,
            start=region.start,
            end=region.end,
            offset=offset,
            left_pad=region.left_pad,
            right_pad=region.right_pad,
        )

    clamped = clamp_region(source, hit.target, flank_start, flank_end)
    if clamped is None:
        logging.debug(f'Empty flank region for hit {hit.idx}, skipping')
        return None
    flank_start, flank_end = clamped

    seq_str = fetch_sequence(source, hit.target, flank_start, flank_end)
    if seq_str is None:
        logging.warning(f'Failed to extract flank sequence for hit {hit.idx}, skipping')
        return None

    if len(seq_str) < flank_len:
        logging.warning(
            f'Flank for hit {hit.idx} on {hit.target} is truncated at '
            f'contig boundary: expected {flank_len}bp, '
            f'extracted {len(seq_str)}bp (coords {flank_start}–{flank_end})'
        )

    return FlankResult(
        seq=seq_str,
        start=flank_start,
        end=flank_end,
        offset=offset,
        left_pad=0,
        right_pad=0,
    )


def writeFlanks(
    outDir: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    model_lengths: Optional[Dict[str, int]] = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    config: Any = None,
    genome: Any = None,
    prefix: Optional[str] = None,
    flank_len: int = 50,
    flank_max_offset: Optional[int] = None,
    write_all: bool = False,
    write_paired: bool = False,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
    pad_flanks: bool = True,
) -> None:
    """
    Extract and write external flanking sequences for terminus hits to FASTA files.

    The "external flank" is the genomic sequence immediately outside each terminus
    hit, i.e. upstream of the left terminus and downstream of the right terminus.
    Flank coordinates are corrected for any gap between the hit alignment and the
    external end of the model (offset correction).

    Parameters
    ----------
    outDir : str, optional
        Output directory for flank FASTA files.
    hitTable : pandas.DataFrame
        DataFrame with columns model, target, hitStart, hitEnd, strand, evalue,
        hmmStart, hmmEnd. Used to look up HMM alignment coordinates by hit index.
    model_lengths : dict
        Dictionary mapping model name to model length in positions.
    paired : dict
        Dictionary of paired hits: paired[model] = [list of pair sets {id1, id2}].
    hitIndex : dict
        Hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig
        Configuration with orientation and model assignments.
    genome : pyfaidx.Fasta, optional
        Indexed genome for sequence extraction. Required if blastdb is None.
    prefix : str, optional
        Prefix for output filenames and sequence IDs.
    flank_len : int, default 50
        Number of bases to extract in each flanking region.
    flank_max_offset : int, optional
        Maximum allowed offset (uncovered model positions) between hit alignment
        and model end. Hits with offset > this value are skipped.
    write_all : bool, default False
        If True, write flanks for all hits (paired and unpaired) to output files.
        For symmetric same-strand orientations (F,F or R,R), both left and right
        flanks are written to separate files with a warning.
    write_paired : bool, default False
        If True, write outer flanks for termini assigned to pairs to separate files.
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database. Alternative to genome.
    pad_flanks : bool, default True
        Pad flanks with N where they run past a contig boundary, so every
        record is ``flank_len`` bases. Padded records carry a ``padded:L,R``
        field in their description. When False, such flanks are truncated.

    Returns
    -------
    None
        Writes FASTA files to disk.

    Notes
    -----
    Output files are named:
      {prefix}{model}_left_flank_{count}.fasta  – flanks for left terminus hits
      {prefix}{model}_right_flank_{count}.fasta – flanks for right terminus hits
      {prefix}{model}_paired_left_flank_{count}.fasta – paired left flanks
      {prefix}{model}_paired_right_flank_{count}.fasta – paired right flanks

    For asymmetric pairings the left and right model names may differ, so each
    model gets its own pair of files. For symmetric pairings the same model name
    is used for both files (distinguished by the _left_/_right_ suffix).

    For symmetric same-strand orientations (F,F or R,R) when write_all is True,
    both left and right flanks are written for all hits to separate files and a
    warning is raised advising the user to use --flanks-paired instead.

    Flanking sequences are always reported in the forward (+) genomic strand
    orientation. Coordinates are 1-based.
    """
    assert outDir is not None, 'outDir cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert model_lengths is not None, 'model_lengths cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    source = make_source(
        genome=genome, blastdb=blastdb, descriptions=genome_descriptions
    )

    # When config is None (e.g. pairing-map mode with multiple configs) unpaired
    # hits cannot be attributed to a terminus; fall back to paired-only processing.
    if config is None and write_all:
        logging.debug(
            'config is None: cannot determine terminus type for unpaired hits; '
            'processing paired hits only.'
        )
        write_all = False
        write_paired = True

    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''

    # left_flanks[model] and right_flanks[model] accumulate SeqRecords
    left_flanks: Dict[str, List[Any]] = {}
    right_flanks: Dict[str, List[Any]] = {}

    # ------------------------------------------------------------------
    # Helper: look up hmmStart / hmmEnd for a hit by its DataFrame index
    # ------------------------------------------------------------------
    def get_hmm_coords(hit_idx: int) -> Tuple[Optional[int], Optional[int]]:
        """
        Retrieve HMM alignment coordinates for a hit.

        Parameters
        ----------
        hit_idx : int
            DataFrame row index for the hit.

        Returns
        -------
        tuple of (int or None, int or None)
            (hmmStart, hmmEnd) as integers, or (None, None) if unavailable.
        """
        if hit_idx not in hitTable.index:  # type: ignore[union-attr]
            return None, None
        row = hitTable.loc[hit_idx]  # type: ignore[union-attr]
        try:
            h_start = int(row['hmmStart'])
            h_end = int(row['hmmEnd'])
        except (ValueError, TypeError):
            return None, None
        return h_start, h_end

    # ------------------------------------------------------------------
    # Helper: retrieve hit record from nested or flat hitIndex
    # ------------------------------------------------------------------
    is_nested = bool(hitIndex) and isinstance(next(iter(hitIndex.values())), dict)

    def get_hit_record(hit_id: int) -> Any:
        """
        Retrieve hit record from hitIndex (handles nested or flat structure).

        Parameters
        ----------
        hit_id : int
            Index of the hit record.

        Returns
        -------
        namedtuple
            Hit record with model, target, hitStart, hitEnd, strand, idx, evalue.
        """
        if is_nested:
            for _m, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in hitIndex')
        return hitIndex[hit_id]['rec']  # type: ignore[index]

    # ------------------------------------------------------------------
    # Helper: extract a flank SeqRecord for one hit
    # ------------------------------------------------------------------
    def build_flank_record(
        hit: Any,
        is_left: bool,
        record_id: str,
        role: Optional[str] = None,
    ) -> Optional[SeqRecord]:
        """
        Build a SeqRecord for the external flank of a single terminus hit.

        Parameters
        ----------
        hit : namedtuple
            Hit record with model, target, hitStart, hitEnd, strand, idx.
        is_left : bool
            True if the hit is a left terminus; False for right terminus.
        record_id : str
            Identifier to assign to the resulting SeqRecord.

        role : str, optional
            Terminus role, 'left' or 'right', used for the _L/_R suffix.
            Defaults to the genomic side given by ``is_left``, which is the
            same thing except for reverse-inserted asymmetric elements.

        Returns
        -------
        SeqRecord or None
            The flank SeqRecord, or None if the flank cannot be extracted
            (missing model length, missing HMM coords, offset exceeds max,
            or a region that does not overlap the contig).
        """
        hmm_start, hmm_end = get_hmm_coords(hit.idx)
        flank = extract_terminus_flank(
            source=source,
            hit=hit,
            is_left=is_left,
            model_len=model_lengths.get(hit.model) if model_lengths else None,  # type: ignore[union-attr]
            hmm_start=hmm_start,
            hmm_end=hmm_end,
            flank_len=flank_len,
            flank_max_offset=flank_max_offset,
            pad=pad_flanks,
        )
        if flank is None:
            return None

        record = SeqRecord(Seq.Seq(flank.seq))
        # The _L/_R suffix names the terminus ROLE. It defaults to the genomic
        # side, which is the same thing for symmetric elements and for forward
        # insertions; callers pass an explicit role where the two can differ.
        effective_role = role if role is not None else ('left' if is_left else 'right')
        side = 'L' if effective_role == 'left' else 'R'
        record.id = f'{record_id}_{side}'
        record.name = record.id

        coord_info = (
            f'[{hit.target}:+ {flank.start}_{flank.end}'
            f' hit:{hit.strand}:{hit.hitStart}_{hit.hitEnd}'
            f' offset:{flank.offset}'
        )
        if flank.is_padded:
            coord_info += f' padded:{flank.left_pad},{flank.right_pad}'
        coord_info += ']'

        record.description = annotate(
            source, hit.target, coord_info, genome_descriptions
        )

        return record

    # ------------------------------------------------------------------
    # Process paired hits
    # ------------------------------------------------------------------
    paired_hit_ids: Set[int] = set()

    # Separate accumulators for paired-only flanks (with element IDs)
    paired_left_flanks: Dict[str, List[Any]] = {}
    paired_right_flanks: Dict[str, List[Any]] = {}

    def _make_paired_flank_record(
        source_rec: SeqRecord, element_id: str, pair_id: str, suffix: str
    ) -> SeqRecord:
        """
        Create a paired-only flank record with element ID in the header.

        Parameters
        ----------
        source_rec : SeqRecord
            Record to copy the sequence and description from.
        element_id : str
            Element identifier for the new record ID.
        pair_id : str
            Pair identifier for the new record ID.
        suffix : str
            Trailing component of the new record ID, e.g. 'L' or 'R'.

        Returns
        -------
        SeqRecord
            A new record carrying the element ID in its header.
        """
        rec = SeqRecord(source_rec.seq)
        rec.id = f'{element_id}_{pair_id}_{suffix}'
        rec.name = rec.id
        rec.description = source_rec.description
        return rec

    for model in paired.keys():
        model_counter = 0
        for pair in paired[model]:
            model_counter += 1
            hit_ids = list(pair)
            x_id, y_id = hit_ids[0], hit_ids[1]
            x = get_hit_record(x_id)
            y = get_hit_record(y_id)
            leftHit, rightHit = flipTIRs(x, y)
            pair_id = f'{model}_{model_counter}'
            element_id = f'Element_{model_counter}'

            # Geometry is positional: the lower-coordinate hit's outer edge
            # faces lower coordinates. That stays true whichever way the
            # element is inserted, so the extracted sequence is unaffected by
            # the routing decision below.
            left_role, right_role = _pair_roles(leftHit, rightHit, config)
            left_rec = build_flank_record(
                leftHit, is_left=True, record_id=pair_id, role=left_role
            )
            right_rec = build_flank_record(
                rightHit, is_left=False, record_id=pair_id, role=right_role
            )

            # Routing is by terminus ROLE, so each model's termini always land
            # in that model's file. For a reverse-inserted asymmetric element
            # the roles are swapped relative to genomic order; grouping by
            # position would put the right model's terminus in the left
            # model's file, mixing the two models' sequences in one output.
            for rec, hit, role in (
                (left_rec, leftHit, left_role),
                (right_rec, rightHit, right_role),
            ):
                if not rec:
                    continue
                model_key = (
                    (config.left_model if role == 'left' else config.right_model)
                    if config is not None and config.is_asymmetric
                    else hit.model
                )
                suffix = 'L' if role == 'left' else 'R'
                if role == 'left':
                    left_flanks.setdefault(model_key, []).append(rec)
                    paired_left_flanks.setdefault(model_key, []).append(
                        _make_paired_flank_record(rec, element_id, pair_id, suffix)
                    )
                else:
                    right_flanks.setdefault(model_key, []).append(rec)
                    paired_right_flanks.setdefault(model_key, []).append(
                        _make_paired_flank_record(rec, element_id, pair_id, suffix)
                    )

            paired_hit_ids.add(leftHit.idx)
            paired_hit_ids.add(rightHit.idx)

    # ------------------------------------------------------------------
    # Process unpaired hits (only when write_all=True)
    # ------------------------------------------------------------------
    # Determine if this is a symmetric same-strand pairing (F,F or R,R)
    is_symmetric_same_strand = (
        config is not None
        and not config.is_asymmetric
        and config.left_strand == config.right_strand
    )

    if write_all and is_symmetric_same_strand:
        logging.warning(
            'Symmetric same-strand orientation detected (%s). '
            'Cannot determine the outer edge for unpaired hits from a single model. '
            'Unpaired hits will be skipped in --flanks output. '
            'Use --flanks-paired for external flanks of confirmed paired termini.',
            ','.join(config.orientation),
        )

    # For asymmetric pairings the hitIndex may contain hits for models that
    # belong to *other* pairs in a multi-pair pairing-map run.  These foreign
    # models should be silently skipped here because their terminus type is
    # determined correctly in their own pair's writeFlanks call.
    config_models: Optional[Set[str]] = None
    if config is not None and config.is_asymmetric:
        # is_asymmetric guarantees both left_model and right_model are non-None,
        # but we filter defensively in case of any future config construction path.
        config_models = {
            m for m in (config.left_model, config.right_model) if m is not None
        }

    if write_all:
        for model in hitIndex.keys():
            # Skip models that do not belong to the current asymmetric pair.
            if config_models is not None and model not in config_models:
                continue

            for hit_id, hit_data in hitIndex[model].items():
                if hit_data['partner'] is not None:
                    continue  # already handled above
                hit = hit_data['rec']

                if is_symmetric_same_strand:
                    # For F,F or R,R symmetric pairings, we cannot determine which
                    # side of an unpaired hit is the external (outer) flank without
                    # knowing whether it is the left or right terminus.  Writing both
                    # flanks would include an internal flank, so we skip unpaired hits
                    # entirely and advise the user to use --flanks-paired instead.
                    logging.debug(
                        f'Skipping unpaired hit {hit_id} (model={hit.model}): '
                        'cannot determine external flank in symmetric same-strand mode'
                    )
                    continue
                else:
                    terminus = resolve_terminus(hit, config)
                    if terminus is None:
                        logging.debug(
                            f'Cannot determine terminus type for unpaired hit {hit.idx} '
                            f'(model={hit.model}, strand={hit.strand}), skipping'
                        )
                        continue

                    # Geometry follows the genomic side the outer edge faces;
                    # grouping follows the terminus role. These differ for a
                    # reverse-inserted asymmetric element, where the left
                    # model's hit lies at the higher coordinate.
                    record_id = f'{model}_{hit_id}_unpaired'
                    rec = build_flank_record(
                        hit, is_left=terminus.is_lower, record_id=record_id
                    )
                    if rec:
                        if terminus.role == 'left':
                            left_flanks.setdefault(hit.model, []).append(rec)
                        else:
                            right_flanks.setdefault(hit.model, []).append(rec)

    # ------------------------------------------------------------------
    # Write output files (all flanks) - only when write_all=True
    # ------------------------------------------------------------------
    if write_all:
        for model, flanks in left_flanks.items():
            if flanks:
                outfile = os.path.join(
                    outDir,
                    prefix + model + '_left_flank_' + str(len(flanks)) + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in flanks:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')

        for model, flanks in right_flanks.items():
            if flanks:
                outfile = os.path.join(
                    outDir,
                    prefix + model + '_right_flank_' + str(len(flanks)) + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in flanks:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')

    # ------------------------------------------------------------------
    # Write paired-only flank files (with element IDs in headers)
    # Only when write_paired=True
    # ------------------------------------------------------------------
    if write_paired:
        for model, flanks in paired_left_flanks.items():
            if flanks:
                outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_left_flank_'
                    + str(len(flanks))
                    + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in flanks:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')

        for model, flanks in paired_right_flanks.items():
            if flanks:
                outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_right_flank_'
                    + str(len(flanks))
                    + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in flanks:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')
