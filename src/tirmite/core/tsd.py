"""Target site duplication (TSD) reconstruction and comparison.

Reconstructs the target site an element inserted into by joining the sequence
flanking each terminus, and compares the duplicated segments to check that the
insertion carries the TSD its family is expected to produce.
"""

import logging
import os
from typing import Any, Dict, List, Optional, Set, Tuple

from Bio import Seq  # type: ignore[import-not-found]
from Bio.SeqRecord import SeqRecord  # type: ignore[import-not-found]
import pandas as pd  # type: ignore[import-untyped]

from tirmite.core.flanks import (
    FlankResult,
    compute_inner_tsd_coordinates,
    extract_terminus_flank,
)
from tirmite.core.termini import _pair_roles, flipTIRs
from tirmite.utils.extract import (
    fetch_region_padded,
    make_source,
)
from tirmite.utils.utils import cleanID


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Compute the Hamming distance between two equal-length strings.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : str
        Second sequence.

    Returns
    -------
    int
        Number of positions at which the corresponding characters differ.

    Raises
    ------
    ValueError
        If the sequences are not of equal length.
    """
    if len(seq1) != len(seq2):
        raise ValueError(
            f'Sequences must be equal length, got {len(seq1)} and {len(seq2)}'
        )
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def load_tsd_length_map(tsd_map_file: str) -> Dict[str, int]:
    """
    Load TSD (Target Site Duplication) lengths from a tab-delimited file.

    Parameters
    ----------
    tsd_map_file : str
        Path to tab-delimited file mapping model pair keys to TSD lengths.
        Format: left_model<TAB>right_model<TAB>tsd_length.

    Returns
    -------
    dict
        Dictionary mapping 'left_model<TAB>right_model' keys to integer TSD lengths.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file format is invalid.
    """
    tsd_lengths: Dict[str, int] = {}

    try:
        with open(tsd_map_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) != 3:
                    raise ValueError(
                        f'Invalid format on line {line_num}: expected 3 tab-delimited '
                        f'columns (left_model, right_model, tsd_length), got {len(parts)}'
                    )

                left_model, right_model, tsd_len_str = (
                    parts[0].strip(),
                    parts[1].strip(),
                    parts[2].strip(),
                )
                try:
                    tsd_len = int(tsd_len_str)
                except ValueError:
                    raise ValueError(
                        f'Invalid TSD length on line {line_num}: {tsd_len_str}'
                    ) from None

                if tsd_len < 0:
                    raise ValueError(
                        f'TSD length must be non-negative on line {line_num}: {tsd_len}'
                    )

                key = f'{left_model}\t{right_model}'
                tsd_lengths[key] = tsd_len

    except FileNotFoundError:
        raise FileNotFoundError(
            f'TSD length map file not found: {tsd_map_file}'
        ) from None

    if not tsd_lengths:
        raise ValueError(f'No valid TSD lengths found in {tsd_map_file}')

    logging.info(
        f'Loaded TSD lengths for {len(tsd_lengths)} model pair(s) from {tsd_map_file}'
    )
    return tsd_lengths


def reconstruct_target_site(
    left_flank_seq: str,
    right_flank_seq: str,
    tsd_length: int = 0,
    tsd_in_model: bool = False,
) -> Tuple[str, str, str, Optional[int]]:
    """
    Reconstruct a target site by joining left and right flanking sequences.

    When a TSD (Target Site Duplication) is present, it is de-duplicated
    so only one copy appears in the reconstructed target site.

    This function handles the ``tsd_in_model=False`` case, where the TSD
    is encoded as the innermost ``tsd_length`` bases of the flanking
    sequences (last ``tsd_length`` bases of the left flank and first
    ``tsd_length`` bases of the right flank).  One copy is trimmed before
    joining.

    .. note::

       When ``tsd_in_model=True`` the TSD lies *inside* the terminus hit
       (part of the HMM model), not in the external flanking sequence.
       :func:`writeTargetSites` handles this case by calling
       :func:`compute_inner_tsd_coordinates` to extract the TSD from
       genomic coordinates within the hit before joining the clean flanks.
       Passing ``tsd_in_model=True`` to this function is therefore a
       no-op — it is treated identically to ``tsd_in_model=False``, both
       operating on the boundary of the supplied flank sequences.

    Parameters
    ----------
    left_flank_seq : str
        External flanking sequence upstream of the left terminus.
    right_flank_seq : str
        External flanking sequence downstream of the right terminus.
    tsd_length : int, default 0
        Length of the terminal duplication feature.
    tsd_in_model : bool, default False
        Retained for API compatibility.  Has no effect on the output
        because this function always treats the innermost bases of the
        flanks as the TSD.  When the TSD is inside the model, the correct
        approach is to use :func:`writeTargetSites` which extracts the TSD
        from the hit genomic coordinates via
        :func:`compute_inner_tsd_coordinates`.

    Returns
    -------
    target_site : str
        Reconstructed target site sequence.
    left_tsd : str
        TSD sequence extracted from the left side (empty string if
        ``tsd_length`` is 0).
    right_tsd : str
        TSD sequence extracted from the right side (empty string if
        ``tsd_length`` is 0).
    tsd_hamming : int or None
        Hamming distance between the left and right TSD sequences over
        informative (non-padded) positions. ``0`` when ``tsd_length`` is 0.
        ``None`` when the TSDs cannot be compared — unequal lengths, or every
        position padded at a contig boundary — which callers must render as
        "unverified" rather than as a distance of 0.

    Notes
    -----
    TSD is at the inner boundary of each flank:
        - ``left_tsd``  = last ``tsd_length`` bases of ``left_flank_seq``
        - ``right_tsd`` = first ``tsd_length`` bases of ``right_flank_seq``

    One copy (the right TSD) is trimmed before joining, yielding:
        ``target_site = left_flank_seq + right_flank_seq[tsd_length:]``
    """
    left_tsd = ''
    right_tsd = ''
    tsd_hamming = 0

    if tsd_length > 0:
        # For both tsd_in_model modes, the TSD appears at the inner boundary
        # of the flanks: the tail of the left flank and the head of the right
        # flank. The distinction (tsd_in_model vs not) affects how the user
        # interprets the duplication relative to the termini model, but the
        # trimming logic is the same: remove one copy from the right flank.
        left_tsd = (
            left_flank_seq[-tsd_length:]
            if len(left_flank_seq) >= tsd_length
            else left_flank_seq
        )
        right_tsd = (
            right_flank_seq[:tsd_length]
            if len(right_flank_seq) >= tsd_length
            else right_flank_seq
        )
        # Trim one TSD copy from the right flank to de-duplicate. The min() is
        # for clarity only - slicing past the end of a shorter flank already
        # yields '' - but it makes the intent explicit: never trim more than is
        # there.
        trimmed_right = right_flank_seq[min(tsd_length, len(right_flank_seq)) :]
        target_site = left_flank_seq + trimmed_right

        # None means the TSDs could not be compared (unequal length, or every
        # position padded). It must stay None rather than collapse to 0, which
        # would report an unverifiable TSD as a perfect duplication.
        tsd_hamming, _compared = compare_tsds(left_tsd, right_tsd)
        if tsd_hamming is None:
            logging.warning(
                f'TSD could not be verified: left={left_tsd or "-"} '
                f'right={right_tsd or "-"}'
            )
        elif tsd_hamming > 0:
            logging.warning(
                f'TSD mismatch (hamming={tsd_hamming}): '
                f'left={left_tsd} right={right_tsd}'
            )
    else:
        target_site = left_flank_seq + right_flank_seq

    return target_site, left_tsd, right_tsd, tsd_hamming


def compare_tsds(
    left_tsd: str,
    right_tsd: str,
    pad_char: str = 'N',
) -> Tuple[Optional[int], int]:
    """
    Compare two TSD sequences, ignoring positions that were padded.

    Parameters
    ----------
    left_tsd, right_tsd : str
        TSD sequences of equal length, possibly containing pad characters where
        the region ran past a contig boundary.
    pad_char : str, default 'N'
        Character marking synthesised positions.

    Returns
    -------
    hamming : int or None
        Hamming distance over informative positions, or None if the TSDs cannot
        be compared at all (different lengths, or no informative position).
    compared : int
        Number of positions actually compared.

    Notes
    -----
    A padded position is not a mismatch and is not a match — no base was
    observed there. Counting it either way misreports the evidence: scoring it
    as a match inflates confidence in a duplication that was never seen, and
    scoring it as a mismatch invents a discrepancy. Both are excluded, and a
    ``None`` result means the comparison could not be made and must not be
    rendered as a distance.
    """
    if not left_tsd or not right_tsd or len(left_tsd) != len(right_tsd):
        return None, 0

    pairs = [
        (a, b)
        for a, b in zip(left_tsd.upper(), right_tsd.upper())
        if a != pad_char and b != pad_char
    ]
    if not pairs:
        return None, 0

    return sum(1 for a, b in pairs if a != b), len(pairs)


def format_interleaved_flanks(
    left_flank_seq: str,
    right_flank_seq: str,
    tsd_length: int = 0,
) -> Tuple[str, str]:
    """
    Format left and right flanks as interleaved gap-padded strings.

    Creates two rows where TSD regions overlap visually:
    - Left flank is right-padded with gaps by the length of the right
      flank minus the TSD overlap.
    - Right flank is left-padded by the length of the left flank minus
      the TSD overlap.

    Parameters
    ----------
    left_flank_seq : str
        Left flanking sequence.
    right_flank_seq : str
        Right flanking sequence.
    tsd_length : int, default 0
        Length of the TSD overlap region.

    Returns
    -------
    left_row : str
        Left flank with gap padding on the right.
    right_row : str
        Right flank with gap padding on the left.

    Examples
    --------
    >>> format_interleaved_flanks('AAAAATSD', 'TSDCCCCC', 3)
    ('AAAAATSD-----', '-----TSDCCCCC')
    """
    overlap = min(tsd_length, len(left_flank_seq), len(right_flank_seq))
    right_pad = len(right_flank_seq) - overlap
    left_pad = len(left_flank_seq) - overlap

    left_row = left_flank_seq + '-' * right_pad
    right_row = '-' * left_pad + right_flank_seq

    return left_row, right_row


def writeTargetSites(
    outDir: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    model_lengths: Optional[Dict[str, int]] = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    config: Any = None,
    genome: Any = None,
    prefix: Optional[str] = None,
    flank_len: int = 0,
    flank_max_offset: Optional[int] = None,
    tsd_length: int = 0,
    tsd_in_model: bool = False,
    tsd_length_map: Optional[Dict[str, int]] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
    pad_flanks: bool = True,
) -> None:
    """
    Reconstruct and write target sites for paired terminus hits.

    Extracts external flanking sequences for each pair, de-duplicates the
    TSD feature, and writes reconstructed target sites, interleaved flanks,
    and a summary report.

    Parameters
    ----------
    outDir : str, optional
        Output directory for target site FASTA files.
    hitTable : pandas.DataFrame
        DataFrame with columns model, target, hitStart, hitEnd, strand, evalue,
        hmmStart, hmmEnd.
    model_lengths : dict
        Dictionary mapping model name to model length.
    paired : dict
        Paired hits: paired[model] = [list of pair sets {id1, id2}].
    hitIndex : dict
        Hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig, optional
        Configuration with orientation and model assignments. May be None
        when using pairing-map mode.
    genome : pyfaidx.Fasta, optional
        Indexed genome for sequence extraction.
    prefix : str, optional
        Prefix for output filenames.
    flank_len : int, default 0
        Number of bases to extract in each flanking region.
    flank_max_offset : int, optional
        Maximum allowed offset from model end.
    tsd_length : int, default 0
        Default TSD length (used when tsd_length_map is not provided).
    tsd_in_model : bool, default False
        Whether the TSD is inside the termini model.
    tsd_length_map : dict, optional
        Map of 'left_model\\tright_model' to TSD length.
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to descriptions.
    blastdb : str, optional
        Path to BLAST database.
    pad_flanks : bool, default True
        Pad flanks and in-model TSDs with N where they run past a contig
        boundary. Reconstructed sites built from padded sequence are marked
        with a ``padded=`` field, and a TSD that cannot be verified reports
        ``tsd_hamming=NA`` rather than 0.

    Returns
    -------
    None
        Writes FASTA files to disk.

    Notes
    -----
    Output files, written per model pair:
      {prefix}{pair_label}_target_sites_{N}.fasta – reconstructed target sites
      {prefix}{pair_label}_interleaved_flanks_{N}.fasta – interleaved flanks
    """
    assert outDir is not None, 'outDir cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert model_lengths is not None, 'model_lengths cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    source = make_source(
        genome=genome, blastdb=blastdb, descriptions=genome_descriptions
    )

    if prefix:
        prefix_str = cleanID(prefix) + '_'
    else:
        prefix_str = ''

    target_site_records: List[Any] = []
    interleaved_records: List[Any] = []

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
            (hmmStart, hmmEnd), or (None, None) if unavailable.
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
        Retrieve a hit record from a nested or flat hitIndex.

        Parameters
        ----------
        hit_id : int
            Index of the hit record.

        Returns
        -------
        namedtuple
            Hit record with model, target, hitStart, hitEnd, strand, idx, evalue.

        Raises
        ------
        KeyError
            If hit_id is not present in a nested hitIndex.
        """
        if is_nested:
            for _m, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in hitIndex')
        return hitIndex[hit_id]['rec']  # type: ignore[index]

    # ------------------------------------------------------------------
    # Helper: extract flank sequence for one hit
    # ------------------------------------------------------------------
    def extract_flank(hit: Any, is_left: bool) -> Optional[FlankResult]:
        """
        Extract the external flank for one terminus hit of a pair.

        Parameters
        ----------
        hit : namedtuple
            Hit record with target, hitStart, hitEnd, strand, idx.
        is_left : bool
            True if the hit is the left terminus of the pair.

        Returns
        -------
        FlankResult or None
            The flank, or None if it could not be extracted.
        """
        hmm_start, hmm_end = get_hmm_coords(hit.idx)
        return extract_terminus_flank(
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

    # ------------------------------------------------------------------
    # Helper: resolve TSD length for a pair of models
    # ------------------------------------------------------------------
    def get_tsd_length_for_pair(left_model: str, right_model: str) -> int:
        """
        Resolve the TSD length to use for a pair of models.

        Parameters
        ----------
        left_model, right_model : str
            Model names of the paired termini.

        Returns
        -------
        int
            TSD length from the map, falling back to the default tsd_length.
        """
        if tsd_length_map:
            key = f'{left_model}\t{right_model}'
            if key in tsd_length_map:
                return tsd_length_map[key]
            # Try symmetric key
            key_sym = f'{right_model}\t{left_model}'
            if key_sym in tsd_length_map:
                return tsd_length_map[key_sym]
            logging.warning(
                f'No TSD length found for model pair ({left_model}, {right_model}) '
                f'in TSD length map, using default tsd_length={tsd_length}'
            )
        return tsd_length

    # ------------------------------------------------------------------
    # Helper: extract TSD sequence from the inner boundary of a terminus hit
    # Used when tsd_in_model=True (TSD is part of the termini model, not in flank)
    # ------------------------------------------------------------------
    def extract_inner_tsd(hit: Any, is_left: bool, tsd_len: int) -> Optional[str]:
        """
        Extract the TSD sequence from the inner (element-facing) boundary of a hit.

        When ``tsd_in_model=True`` the TSD is encoded at the element-facing end
        of the terminus HMM and therefore lies *within* the hit coordinates
        rather than in the external flanking sequence.  This helper computes the
        genomic coordinates of the TSD and extracts the bases from the genome or
        BLAST database.

        Parameters
        ----------
        hit : namedtuple
            Hit record (model, target, hitStart, hitEnd, strand, idx).
        is_left : bool
            True if the hit is the left (upstream) terminus of the element.
        tsd_len : int
            Length of the TSD to extract.

        Returns
        -------
        PaddedRegion or None
            The TSD region (always on the forward/+ strand), or None if the
            coordinates cannot be determined or the region does not overlap the
            contig.

        Notes
        -----
        The region is padded to ``tsd_len`` when it runs past a contig
        boundary. A padded TSD is not evidence of a duplication and callers
        must exclude the padded positions from any comparison — see
        ``_compare_tsds``.
        """
        model = hit.model
        model_len = model_lengths.get(model) if model_lengths else None  # type: ignore[union-attr]
        if model_len is None:
            logging.warning(
                f'Model length not found for {model}, skipping TSD extraction'
            )
            return None

        hmm_start, hmm_end = get_hmm_coords(hit.idx)
        if hmm_start is None or hmm_end is None:
            logging.debug(
                f'HMM coordinates unavailable for hit {hit.idx}, skipping TSD'
            )
            return None

        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=int(hit.hitStart),
            hit_end=int(hit.hitEnd),
            strand=hit.strand,
            is_left_terminus=is_left,
            hmm_start=hmm_start,
            hmm_end=hmm_end,
            model_len=model_len,
            tsd_length=tsd_len,
        )

        if tsd_end < tsd_start:
            logging.debug(
                f'Invalid TSD coordinates for hit {hit.idx}: {tsd_start}-{tsd_end}'
            )
            return None

        # Pad rather than clamp. A TSD silently returned short at a contig
        # boundary used to be compared against a full-length partner, and the
        # length mismatch made the comparison be skipped entirely and reported
        # as hamming 0 - a truncated, unverifiable TSD was indistinguishable
        # from a perfect match.
        return fetch_region_padded(source, hit.target, tsd_start, tsd_end)

    # ------------------------------------------------------------------
    # Determine canonical pair key for file naming
    # ------------------------------------------------------------------
    # Use config model assignments when available so that all pairs for the
    # same model combination are grouped into one output file regardless of
    # which model happens to be at lower genomic coordinates.
    if config is not None and config.left_model is not None:
        _canonical_pair_key = f'{config.left_model}\t{config.right_model}'
    else:
        _canonical_pair_key = None  # will be derived per-pair below

    # ------------------------------------------------------------------
    # Process paired hits – group by model pair
    # ------------------------------------------------------------------
    # Records grouped by model-pair key for per-pair output files
    pair_key_records: Dict[str, List[Any]] = {}
    pair_key_interleaved: Dict[str, List[Any]] = {}

    for model in paired.keys():
        model_counter = 0
        for pair in paired[model]:
            model_counter += 1
            hit_ids = list(pair)
            x_id, y_id = hit_ids[0], hit_ids[1]
            x = get_hit_record(x_id)
            y = get_hit_record(y_id)
            leftHit, rightHit = flipTIRs(x, y)
            pair_id = f'{prefix_str}{model}_{model_counter}'

            left_flank = extract_flank(leftHit, is_left=True)
            right_flank = extract_flank(rightHit, is_left=False)

            if left_flank is None or right_flank is None:
                logging.debug(
                    f'Could not extract flanks for pair {pair_id}, skipping target site'
                )
                continue

            left_seq = left_flank.seq
            right_seq = right_flank.seq
            # Any padded base in the flanks, or in an in-model TSD below, makes
            # the reconstructed target site partly synthetic.
            padded = left_flank.is_padded or right_flank.is_padded

            pair_tsd_len = get_tsd_length_for_pair(leftHit.model, rightHit.model)

            if tsd_in_model and pair_tsd_len > 0:
                # TSD is inside the terminus model – extract from the inner boundary
                # of each hit, not from the external flank.
                # Reconstruction: left_flank + TSD + right_flank
                left_region = extract_inner_tsd(
                    leftHit, is_left=True, tsd_len=pair_tsd_len
                )
                right_region = extract_inner_tsd(
                    rightHit, is_left=False, tsd_len=pair_tsd_len
                )
                left_tsd = left_region.seq if left_region else ''
                right_tsd = right_region.seq if right_region else ''
                tsd_padded = bool(
                    (left_region and left_region.is_padded)
                    or (right_region and right_region.is_padded)
                )
                padded = padded or tsd_padded

                # Use left TSD (arbitrarily) to fill the gap; warn if mismatched
                tsd_seq = left_tsd if left_tsd else right_tsd
                target_site = left_seq + tsd_seq + right_seq

                tsd_hamming, tsd_compared = compare_tsds(left_tsd, right_tsd)
                if tsd_hamming is None:
                    logging.warning(
                        f'TSD for pair {pair_id} could not be verified '
                        f'(left={left_tsd or "-"} right={right_tsd or "-"}); '
                        'reporting as unverified rather than as a match'
                    )
                elif tsd_hamming > 0:
                    logging.warning(
                        f'TSD mismatch for pair {pair_id} '
                        f'(hamming={tsd_hamming} over {tsd_compared}bp): '
                        f'left={left_tsd} right={right_tsd}'
                    )
            else:
                # TSD is outside the model – it is the innermost n bases of each flank.
                # reconstruct_target_site() trims one copy and joins.
                target_site, left_tsd, right_tsd, tsd_hamming = reconstruct_target_site(
                    left_flank_seq=left_seq,
                    right_flank_seq=right_seq,
                    tsd_length=pair_tsd_len,
                    tsd_in_model=False,
                )

            # Report the header in terms of terminus ROLE, so left_model and
            # left_flank_hit always describe the same terminus. For a
            # reverse-inserted asymmetric element that is the hit at the higher
            # coordinate, which element_orientation records explicitly.
            left_role, _right_role = _pair_roles(leftHit, rightHit, config)
            reversed_insertion = left_role == 'right'
            role_left_hit = rightHit if reversed_insertion else leftHit
            role_right_hit = leftHit if reversed_insertion else rightHit

            # Build metadata for FASTA header. A TSD that could not be compared
            # is reported as 'NA', never as 0 - a padded or absent TSD must not
            # read as a confirmed perfect duplication.
            meta_parts = [
                f'flank_len={flank_len}',
                f'tsd_len={pair_tsd_len}',
                f'tsd_in_model={tsd_in_model}',
                f'left_model={role_left_hit.model}',
                f'right_model={role_right_hit.model}',
                f'contig={leftHit.target}',
                f'element_orientation={"reverse" if reversed_insertion else "forward"}',
                f'left_flank_hit={role_left_hit.strand}:{role_left_hit.hitStart}_{role_left_hit.hitEnd}',
                f'right_flank_hit={role_right_hit.strand}:{role_right_hit.hitStart}_{role_right_hit.hitEnd}',
                f'tsd_hamming={"NA" if tsd_hamming is None else tsd_hamming}',
            ]
            if padded:
                meta_parts.append(
                    f'padded=left:{left_flank.left_pad},{left_flank.right_pad};'
                    f'right:{right_flank.left_pad},{right_flank.right_pad}'
                )
            if left_tsd:
                meta_parts.append(f'left_tsd={left_tsd}')
            if right_tsd:
                meta_parts.append(f'right_tsd={right_tsd}')

            description = ' '.join(meta_parts)

            ts_record = SeqRecord(
                Seq.Seq(target_site),
                id=pair_id,
                name=pair_id,
                description=description,
            )
            target_site_records.append(ts_record)

            # Group by model pair for per-pair output — use canonical key
            pair_key = (
                _canonical_pair_key
                if _canonical_pair_key is not None
                else f'{leftHit.model}\t{rightHit.model}'
            )
            pair_key_records.setdefault(pair_key, []).append(ts_record)

            # Build interleaved flanks
            left_row, right_row = format_interleaved_flanks(
                left_flank_seq=left_seq,
                right_flank_seq=right_seq,
                tsd_length=pair_tsd_len,
            )

            il_left = SeqRecord(
                Seq.Seq(left_row),
                id=f'{pair_id}_left',
                name=f'{pair_id}_left',
                description=description,
            )
            il_right = SeqRecord(
                Seq.Seq(right_row),
                id=f'{pair_id}_right',
                name=f'{pair_id}_right',
                description=description,
            )
            interleaved_records.append(il_left)
            interleaved_records.append(il_right)

            # Group by model pair for per-pair output — use same canonical key
            pair_key_interleaved.setdefault(pair_key, []).extend([il_left, il_right])

    # ------------------------------------------------------------------
    # Helper: write single-line non-wrapped FASTA
    # ------------------------------------------------------------------
    def _write_single_line_fasta(records: List[Any], filepath: str) -> None:
        """
        Write SeqRecords as single-line non-wrapped FASTA.

        Parameters
        ----------
        records : list
            SeqRecords to write.
        filepath : str
            Destination path.

        Returns
        -------
        None
            Writes to disk.
        """
        with open(filepath, 'w') as handle:
            for rec in records:
                handle.write(f'>{rec.id} {rec.description}\n')
                handle.write(f'{str(rec.seq)}\n')

    # ------------------------------------------------------------------
    # Write output files – one per model pair
    # ------------------------------------------------------------------
    if target_site_records:
        # Write per-pair target site files
        for pair_key, records in pair_key_records.items():
            left_model, right_model = pair_key.split('\t')
            pair_label = (
                f'{left_model}_{right_model}'
                if left_model != right_model
                else left_model
            )
            ts_outfile = os.path.join(
                outDir, f'{prefix_str}{pair_label}_target_sites_{len(records)}.fasta'
            )
            _write_single_line_fasta(records, ts_outfile)
            logging.info(
                f'Wrote {len(records)} reconstructed target sites to {ts_outfile}'
            )
    else:
        logging.warning('No target sites could be reconstructed')

    if interleaved_records:
        # Write per-pair interleaved flank files
        for pair_key, records in pair_key_interleaved.items():
            left_model, right_model = pair_key.split('\t')
            pair_label = (
                f'{left_model}_{right_model}'
                if left_model != right_model
                else left_model
            )
            il_outfile = os.path.join(
                outDir,
                f'{prefix_str}{pair_label}_interleaved_flanks_{len(records)}.fasta',
            )
            _write_single_line_fasta(records, il_outfile)
            logging.info(f'Wrote interleaved flanks to {il_outfile}')
