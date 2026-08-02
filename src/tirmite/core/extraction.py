"""Extraction of TIR and element sequences from a genome.

Pulls the sequence for individual terminus hits, for paired termini, and for
whole elements spanning a pair, and writes each to FASTA.

Notes
-----
All coordinates handled here are 1-based inclusive and plus-strand, matching
the contract documented in :mod:`tirmite.utils.extract`. Sequence is fetched
through a ``SequenceSource``, so the same code serves an indexed FASTA and a
BLAST database.
"""

from collections import namedtuple
import logging
import os
from typing import Any, Dict, List, Optional, Set, Tuple

from Bio import Seq, SeqIO  # type: ignore[import-not-found]
from Bio.SeqRecord import SeqRecord  # type: ignore[import-not-found]
import pandas as pd  # type: ignore[import-untyped]

from tirmite.core.termini import _model_deficit, _pair_roles, flipTIRs
from tirmite.utils.extract import (
    SequenceSource,
    annotate,
    clamp_region,
    fetch_region_padded,
    fetch_sequence,
    make_source,
)
from tirmite.utils.utils import cleanID


def fetch_padded_hit(
    source: SequenceSource,
    seqid: str,
    start: int,
    end: int,
    strand: str = '+',
    padlen: Optional[int] = None,
    pad: bool = True,
) -> Optional[str]:
    """
    Fetch a hit region, optionally with flanking sequence marked in lowercase.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to read from.
    seqid : str
        Sequence identifier.
    start, end : int
        1-based inclusive hit coordinates on the plus strand.
    strand : str, default '+'
        Strand of the hit. The result is reverse-complemented when '-'.
    padlen : int, optional
        Number of flanking bases to include either side of the hit. Flanking
        bases are lowercased; the hit itself stays uppercase.
    pad : bool, default True
        When the padded window runs past a contig boundary, pad it with N so
        the record is always ``(end - start + 1) + 2 * padlen`` bases. Padding
        is lowercased along with the rest of the flank, keeping the "not part
        of the hit" convention intact.

    Returns
    -------
    str or None
        Extracted sequence, or None if extraction failed.

    Notes
    -----
    Case marking is applied to the plus-strand sequence and the reverse
    complement is taken afterwards, so the lowercase markers stay attached to
    the flanks they came from. Marking after reverse-complementing would swap
    the 5' and 3' markers.

    Only the pad region can extend past a contig; the hit itself came from the
    contig and is always in bounds.
    """
    if not padlen:
        seq = fetch_sequence(source, seqid, start, end)
        if seq is None:
            return None
        return seq if strand != '-' else str(Seq.Seq(seq).reverse_complement())

    # An inverted window has no hit region to mark. Padding would still produce
    # a valid-looking span, so reject it here rather than emit a nonsense
    # arrangement of lowercase flanks around nothing.
    if start > end:
        logging.debug(f'Inverted hit window for {seqid}: {start}-{end}, skipping')
        return None

    pad_start = start - padlen
    pad_end = end + padlen

    if pad:
        region = fetch_region_padded(source, seqid, pad_start, pad_end)
        if region is None:
            return None
        full_seq = region.seq
    else:
        # Clamp the padded window so the case-marking offsets below are
        # computed against the region actually returned.
        clamped = clamp_region(source, seqid, pad_start, pad_end)
        if clamped is None:
            return None
        pad_start, pad_end = clamped
        maybe_seq = fetch_sequence(source, seqid, pad_start, pad_end)
        if maybe_seq is None:
            return None
        full_seq = maybe_seq

    # Offsets of the hit within the padded window, as 0-based half-open slice
    # bounds. Both coordinate systems here are 1-based inclusive.
    lead = max(0, start - pad_start)
    tail = lead + (min(end, pad_end) - max(start, pad_start) + 1)

    marked = full_seq[:lead].lower() + full_seq[lead:tail] + full_seq[tail:].lower()

    if strand == '-':
        return str(Seq.Seq(marked).reverse_complement())
    return marked


def _model_extension(
    row: Any,
    model: str,
    model_lengths: Optional[Dict[str, int]],
) -> Optional[Tuple[int, int]]:
    """
    Compute bases to add either side of a hit to span the full model.

    Parameters
    ----------
    row : pandas.Series
        Hit row carrying hmmStart, hmmEnd and strand.
    model : str
        Model name, used to look up the length and for the warning.
    model_lengths : dict or None
        Mapping of model name to model length.

    Returns
    -------
    tuple of (int, int) or None
        ``(lower, upper)`` bases to extend at the lower and higher genomic
        coordinate ends, or None if the model length is unavailable.

    Notes
    -----
    The deficits are the same quantities used for flank projection: how many
    model positions the alignment failed to cover at each end. On the plus
    strand the ``hmmStart`` deficit sits at the lower-coordinate end; on the
    minus strand the mapping is reversed.
    """
    model_len = model_lengths.get(model) if model_lengths else None
    if model_len is None:
        logging.warning(
            f'Model length not found for {model}; cannot extend hits to model '
            'length, extracting the aligned region only'
        )
        return None

    try:
        hmm_start = int(row['hmmStart'])
        hmm_end = int(row['hmmEnd'])
    except (ValueError, TypeError, KeyError):
        logging.debug(f'HMM coordinates unavailable for a {model} hit, not extending')
        return None

    start_deficit = _model_deficit(
        hmm_start - 1, 'hmmStart', hmm_start, hmm_end, model_len
    )
    end_deficit = _model_deficit(
        model_len - hmm_end, 'hmmEnd', hmm_start, hmm_end, model_len
    )

    if row['strand'] == '-':
        # hmmStart aligns to the higher genomic coordinate on the minus strand.
        return end_deficit, start_deficit
    return start_deficit, end_deficit


def extractTIRs(
    model: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    maxeval: float = 0.001,
    source: Optional[SequenceSource] = None,
    padlen: Optional[int] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    model_lengths: Optional[Dict[str, int]] = None,
    extend_to_model: bool = False,
    pad: bool = True,
) -> Tuple[List[SeqRecord], int]:
    """
    Extract TIR sequences for hits of a specific model.

    Parameters
    ----------
    model : str
        Name of HMM model to extract hits for.
    hitTable : pandas.DataFrame
        DataFrame containing all hits with columns: model, target, hitStart, hitEnd,
        strand, evalue, hmmStart, hmmEnd.
    maxeval : float, default 0.001
        Maximum e-value threshold for extracting hits.
    source : FastaSource or BlastDBSource
        Sequence source, from :func:`tirmite.utils.extract.make_source`.
    padlen : int, optional
        Number of flanking bases to extract on each side of hit (shown in lowercase).
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    model_lengths : dict, optional
        Mapping of model name to model length. Required when
        ``extend_to_model`` is set.
    extend_to_model : bool, default False
        Extend each hit outward by the number of model positions it failed to
        cover, so partial hits are emitted at full model length.
    pad : bool, default True
        Pad with N where the extracted window runs past a contig boundary.

    Returns
    -------
    seqList : list of Bio.SeqRecord.SeqRecord
        List of SeqRecord objects containing extracted TIR sequences.
    hitcount : int
        Number of hits extracted that passed e-value threshold.

    Notes
    -----
    Reverse complement is applied for hits on the negative strand, and such
    records get an '_rc' ID suffix. Padded flanking sequence is shown in
    lowercase letters. Coordinates are 1-based inclusive throughout; see
    :mod:`tirmite.utils.extract` for the extraction contract.

    With ``extend_to_model``, the extension is added on both sides using the
    same model-coverage deficits that :func:`compute_flank_coordinates` uses for
    the external edge and :func:`compute_inner_tsd_coordinates` for the inner
    edge. The extended bases were not matched by the model, so they are a
    projection of where the terminus should end, not evidence that it does.
    """
    assert model is not None, 'model cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert source is not None, 'source cannot be None'
    hitcount = 0
    seqList = []

    if genome_descriptions is None:
        genome_descriptions = getattr(source, 'descriptions', None)

    eligible_hits = hitTable[
        (hitTable['model'] == model) & (hitTable['evalue'].astype(float) <= maxeval)
    ]
    total_eligible = len(eligible_hits)
    logging.info(
        f'Extracting {total_eligible} hits for model "{model}" '
        f'from {source.describe()}...'
    )
    _log_step = max(1, min(100, total_eligible // 10)) if total_eligible > 0 else 1

    for index, row in eligible_hits.iterrows():
        hitcount += 1
        if hitcount % _log_step == 0 or hitcount == total_eligible:
            logging.info(
                f'  Extracted {hitcount}/{total_eligible} hits for model "{model}"'
            )

        hit_start = int(row['hitStart'])
        hit_end = int(row['hitEnd'])
        extension = ''

        if extend_to_model:
            span = _model_extension(row, model, model_lengths)
            if span is not None:
                lower_ext, upper_ext = span
                hit_start -= lower_ext
                hit_end += upper_ext
                if lower_ext or upper_ext:
                    extension = f' extended:{lower_ext},{upper_ext}'

        hit_seq_str = fetch_padded_hit(
            source,
            row['target'],
            hit_start,
            hit_end,
            row['strand'],
            padlen,
            pad=pad,
        )

        if hit_seq_str is None:
            logging.warning(f'Failed to extract sequence for {model}_{index}, skipping')
            continue

        # Create SeqRecord
        hitrecord = SeqRecord(Seq.Seq(hit_seq_str))
        hitrecord.id = model + '_' + str(index)

        if row['strand'] == '-':
            # Sequence is already reverse-complemented by fetch_padded_hit;
            # only the ID suffix is applied here.
            hitrecord.id = hitrecord.id + '_rc'

        hitrecord.name = hitrecord.id

        # Build description with genome description. Report the extracted span
        # so the header always describes the sequence actually emitted.
        coord_info = '_'.join(
            [
                '[' + str(row['target']) + ':' + str(row['strand']),
                str(hit_start),
                str(hit_end) + ' modelAlignment:' + str(row['hmmStart']),
                str(row['hmmEnd']) + ' E-value:' + str(row['evalue']) + extension + ']',
            ]
        )

        # Add genome description if available
        hitrecord.description = annotate(
            source, row['target'], coord_info, genome_descriptions
        )

        seqList.append(hitrecord)

    return seqList, hitcount


# Fix: Do not load fasta into genome!
def writeTIRs(
    outDir: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    maxeval: float = 0.001,
    genome: Any = None,
    prefix: Optional[str] = None,
    padlen: Optional[int] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
    model_lengths: Optional[Dict[str, int]] = None,
    extend_hits_to_model: bool = False,
    pad: bool = True,
) -> None:
    """
    Write extracted TIR sequences to FASTA files organized by model.

    Parameters
    ----------
    outDir : str, optional
        Output directory for FASTA files. If None, uses current directory.
    hitTable : pandas.DataFrame
        DataFrame containing all hits with columns: model, target, hitStart, hitEnd,
        strand, evalue, hmmStart, hmmEnd.
    maxeval : float, default 0.001
        Maximum e-value threshold for extracting hits.
    genome : pyfaidx.Fasta, optional
        Indexed genome object for sequence extraction. Required if blastdb is None.
    prefix : str, optional
        Prefix to add to output filenames and sequence IDs.
    padlen : int, optional
        Number of flanking bases to extract on each side of hit (shown in lowercase).
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database for sequence extraction. Alternative to genome.
    model_lengths : dict, optional
        Mapping of model name to model length; required for
        ``extend_hits_to_model``.
    extend_hits_to_model : bool, default False
        Emit partial hits at full model length by extending them outward by the
        uncovered model positions.
    pad : bool, default True
        Pad with N where an extracted window runs past a contig boundary.

    Returns
    -------
    None
        Writes FASTA files to disk but returns nothing.

    Notes
    -----
    Creates one FASTA file per model with filename format:
    {prefix}{model}_hits_{count}.fasta

    Either genome or blastdb must be provided for sequence extraction.

    Set ``extend_hits_to_model`` to emit partial hits at full model length;
    ``model_lengths`` is required for that and the option is ignored with a
    warning without it.
    """
    assert hitTable is not None, 'hitTable cannot be None'

    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''
    if outDir:
        outDir = os.path.abspath(outDir)
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
    else:
        outDir = os.getcwd()

    source = make_source(genome=genome, blastdb=blastdb)

    for model in hitTable['model'].unique():
        # List of TIR seqrecords, and count of hits
        seqList, hitcount = extractTIRs(
            model=model,
            hitTable=hitTable,
            maxeval=maxeval,
            source=source,
            padlen=padlen,
            genome_descriptions=genome_descriptions,
            model_lengths=model_lengths,
            extend_to_model=extend_hits_to_model,
            pad=pad,
        )
        outfile = os.path.join(
            outDir, prefix + model + '_hits_' + str(hitcount) + '.fasta'
        )
        # Write extracted hits to model outfile
        with open(outfile, 'w') as handle:
            for seq in seqList:
                seq.id = prefix + str(seq.id)
                SeqIO.write(seq, handle, 'fasta')


def fetchElements(
    paired: Optional[Dict[str, List[Set[int]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    genome: Any = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
) -> Dict[str, List[Any]]:
    """
    Extract full-length transposon element sequences from paired TIR hits.

    Parameters
    ----------
    paired : dict
        Dictionary of paired hits: paired[model] = [list of pair sets {id1, id2}].
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.
        Supports both nested (model-keyed) and flat structures.
    genome : pyfaidx.Fasta, optional
        Indexed genome object for sequence extraction. Required if blastdb is None.
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database for sequence extraction. Alternative to genome.

    Returns
    -------
    dict
        Dictionary of element records keyed by model: TIRelements[model] = [element_records].
        Each element record is a namedtuple with fields: model, chromosome, start, end,
        strand, type, id, leftHit, rightHit, seq, evalue.

    Notes
    -----
    Extracts sequence from leftHit.hitStart to rightHit.hitEnd.
    Handles both symmetric (same model) and asymmetric (different models) pairing.
    Element IDs have format: Element_{counter}.

    Either genome or blastdb must be provided for sequence extraction.
    """
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    source = make_source(
        genome=genome, blastdb=blastdb, descriptions=genome_descriptions
    )
    # Check if hitIndex is nested or flat
    is_nested = isinstance(next(iter(hitIndex.values())), dict)

    def get_hit_record(hit_id: int) -> Any:
        """
        Retrieve hit record from either nested or flat hitIndex structure.

        Parameters
        ----------
        hit_id : int
            Index of the hit record to retrieve.

        Returns
        -------
        namedtuple
            Hit record namedtuple with fields: model, target, hitStart, hitEnd,
            strand, idx, evalue.

        Raises
        ------
        KeyError
            If hit_id is not found in any model within the hitIndex.

        Notes
        -----
        Handles both nested (model-keyed) and flat hitIndex structures for
        backward compatibility.
        """
        if is_nested:
            # Search through all models to find the hit
            for _model_name, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in any model')
        else:
            # Flat structure - direct access
            return hitIndex[hit_id]['rec']  # type: ignore[index]

    TIRelements: Dict[str, List[Any]] = {}
    gffTup = namedtuple(
        'gffTup',  # type: ignore[name-match]
        [
            'model',
            'chromosome',
            'start',
            'end',
            'strand',
            'type',
            'id',
            'leftHit',
            'rightHit',
            'seq',
            'evalue',
        ],
    )

    # Only process models that actually have pairs
    for model in paired.keys():
        if len(paired[model]) > 0:
            TIRelements[model] = []
            model_counter = 0
            total_pairs = len(paired[model])
            logging.info(
                f'Extracting sequences for {total_pairs} paired elements (model "{model}")...'
            )
            _log_step = max(1, min(100, total_pairs // 10)) if total_pairs > 0 else 1

            for pair in paired[model]:
                model_counter += 1
                if model_counter % _log_step == 0 or model_counter == total_pairs:
                    logging.info(
                        f'  Extracted {model_counter}/{total_pairs} elements for model "{model}"'
                    )
                # Convert set to list for indexing
                hit_ids = list(pair)
                x_id, y_id = hit_ids[0], hit_ids[1]

                # Get hit records using helper function
                x = get_hit_record(x_id)
                y = get_hit_record(y_id)

                leftHit, rightHit = flipTIRs(x, y)

                # Create element ID with counter only (avoids model name duplication)
                eleID = f'Element_{model_counter}'

                # Extract element sequence. Always on the + strand so that the
                # sequence matches the genomic coordinates reported in the GFF,
                # consistent with the flank and target-site extractors.
                ele_seq_str = fetch_sequence(
                    source,
                    leftHit.target,
                    int(leftHit.hitStart),
                    int(rightHit.hitEnd),
                )
                if ele_seq_str is None:
                    logging.warning(f'Failed to extract element {eleID}, skipping')
                    continue

                eleSeq = SeqRecord(Seq.Seq(ele_seq_str))
                eleSeq.id = eleID
                eleSeq.name = eleID

                # Build description with genome description
                coord_info = (
                    '_'.join(
                        [
                            '[' + leftHit.target + ':' + str(leftHit.hitStart),
                            str(rightHit.hitEnd),
                        ]
                    )
                    + ' len='
                    + str(rightHit.hitEnd - leftHit.hitStart)
                    + ']'
                )

                # Add genome description if available
                eleSeq.description = annotate(
                    source, leftHit.target, coord_info, genome_descriptions
                )

                TIRelement = gffTup(
                    model,  # This is the pairing model (left model for asymmetric)
                    leftHit.target,
                    leftHit.hitStart,
                    rightHit.hitEnd,
                    leftHit.strand,
                    'Element',
                    eleID,
                    leftHit,
                    rightHit,
                    eleSeq,
                    'NA',
                )
                TIRelements[model].append(TIRelement)

    return TIRelements


def writeElements(
    outDir: str,
    eleDict: Optional[Dict[str, List[Any]]] = None,
    prefix: Optional[str] = None,
) -> None:
    """
    Write extracted element sequences to FASTA files organized by model.

    Parameters
    ----------
    outDir : str
        Output directory for element FASTA files.
    eleDict : dict, optional
        Dictionary of element records keyed by model: eleDict[model] = [element_records].
    prefix : str, optional
        Prefix to add to output filenames and sequence IDs.

    Returns
    -------
    None
        Writes FASTA files to disk but returns nothing.

    Notes
    -----
    Only creates files for models that have at least one element.
    Output filename format: {prefix}{model}_elements_{count}.fasta
    """
    assert eleDict is not None, 'eleDict cannot be None'
    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''

    for model in eleDict.keys():
        if len(eleDict[model]) > 0:  # Only write files for models with actual elements
            count = len(eleDict[model])
            outfile = os.path.join(
                outDir,
                prefix + model + '_elements_' + str(count) + '.fasta',
            )
            with open(outfile, 'w') as handle:
                for element in eleDict[model]:
                    element.seq.id = prefix + str(element.seq.id)
                    SeqIO.write(element.seq, handle, 'fasta')


def writePairedTIRs(
    outDir: Optional[str] = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    genome: Any = None,
    prefix: Optional[str] = None,
    padlen: Optional[int] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
    config: Any = None,
) -> None:
    """
    Extract and write left and right TIR sequences from paired hits to FASTA.

    Parameters
    ----------
    outDir : str, optional
        Output directory for TIR FASTA files.
    paired : dict
        Dictionary of paired hits: paired[model] = [list of pair sets {id1, id2}].
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.
        Supports both nested (model-keyed) and flat structures.
    genome : pyfaidx.Fasta, optional
        Indexed genome object for sequence extraction. Required if blastdb is None.
    prefix : str, optional
        Prefix to add to output filenames and sequence IDs.
    padlen : int, optional
        Number of flanking bases to extract on each side (shown in lowercase).
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database for sequence extraction. Alternative to genome.

    config : PairingConfig, optional
        Pairing configuration. When it describes an asymmetric element, the
        _L/_R labels and output files follow the terminus role rather than
        genomic order, so each model's termini stay together.

    Returns
    -------
    None
        Writes FASTA files to disk but returns nothing.

    Notes
    -----
    Right TIRs are reverse complemented. Outputs are split into two files per model:
    {model}_{counter}_L (left TIR) → {prefix}{model}_paired_left_term_hits_{count}.fasta
    {model}_{counter}_R (right TIR) → {prefix}{model}_paired_right_term_hits_{count}.fasta

    Either genome or blastdb must be provided for sequence extraction.
    """
    assert outDir is not None, 'outDir cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    source = make_source(
        genome=genome, blastdb=blastdb, descriptions=genome_descriptions
    )
    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''

    # Check if hitIndex is nested (new format) or flat (old format)
    is_nested = isinstance(next(iter(hitIndex.values())), dict)

    def get_hit_record(hit_id: int) -> Any:
        """
        Retrieve hit record from either nested or flat hitIndex structure.

        Parameters
        ----------
        hit_id : int
            Index of the hit record to retrieve.

        Returns
        -------
        namedtuple
            Hit record namedtuple with fields: model, target, hitStart, hitEnd,
            strand, idx, evalue.

        Raises
        ------
        KeyError
            If hit_id is not found in any model within the hitIndex.

        Notes
        -----
        Handles both nested (model-keyed) and flat hitIndex structures for
        backward compatibility.
        """
        if is_nested:
            # Search through all models to find the hit
            for _model_name, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in any model')
        else:
            # Flat structure - direct access
            return hitIndex[hit_id]['rec']  # type: ignore[index]

    # Only process models that actually have pairs
    for model in paired.keys():
        if len(paired[model]) > 0:  # Only write files for models with actual pairs
            model_counter = 0
            left_seqList: List[Any] = []  # Left terminus hit sequences
            right_seqList: List[Any] = []  # Right terminus hit sequences
            total_pairs = len(paired[model])
            logging.info(
                f'Extracting TIR sequences for {total_pairs} pairs (model "{model}")...'
            )
            _log_step = max(1, min(100, total_pairs // 10)) if total_pairs > 0 else 1

            for pair in paired[model]:
                model_counter += 1
                if model_counter % _log_step == 0 or model_counter == total_pairs:
                    logging.info(
                        f'  Extracted {model_counter}/{total_pairs} TIR pairs for model "{model}"'
                    )
                # Convert set to list for indexing
                hit_ids = list(pair)
                x_id, y_id = hit_ids[0], hit_ids[1]

                # Get hit records using helper function
                x = get_hit_record(x_id)
                y = get_hit_record(y_id)

                leftHit, rightHit = flipTIRs(x, y)
                eleID = f'{model}_{model_counter}'

                # Extract TIR sequences. Both are taken on the + strand; the
                # right TIR is reverse complemented below so that the pair reads
                # 5'->3' inward.
                left_seq_str = fetch_padded_hit(
                    source,
                    leftHit.target,
                    int(leftHit.hitStart),
                    int(leftHit.hitEnd),
                    '+',
                    padlen,
                )
                right_seq_str = fetch_padded_hit(
                    source,
                    rightHit.target,
                    int(rightHit.hitStart),
                    int(rightHit.hitEnd),
                    '+',
                    padlen,
                )

                if left_seq_str is None or right_seq_str is None:
                    logging.warning(f'Failed to extract TIRs for {eleID}, skipping')
                    continue

                # Create SeqRecords for FASTA output only
                eleSeqLeft = SeqRecord(Seq.Seq(left_seq_str))
                eleSeqRight = SeqRecord(Seq.Seq(right_seq_str))
                eleSeqRight = eleSeqRight.reverse_complement()

                # Sequence orientation is positional (the higher-coordinate
                # terminus is reverse complemented so both read inward), but
                # the _L/_R label and the output file follow terminus role, so
                # each model's termini stay together for downstream alignment.
                left_role, right_role = _pair_roles(leftHit, rightHit, config)
                left_label = 'L' if left_role == 'left' else 'R'
                right_label = 'L' if right_role == 'left' else 'R'

                eleSeqLeft.id = eleID + '_' + left_label
                eleSeqLeft.name = eleSeqLeft.id

                # Build left description with genome description
                left_coord = (
                    '_'.join(
                        [
                            '[' + leftHit.target + ':' + str(leftHit.hitStart),
                            str(leftHit.hitEnd),
                        ]
                    )
                    + ']'
                )

                eleSeqLeft.description = annotate(
                    source, leftHit.target, left_coord, genome_descriptions
                )

                eleSeqRight.id = eleID + '_' + right_label
                eleSeqRight.name = eleSeqRight.id

                # Build right description with genome description
                right_coord = (
                    '_'.join(
                        [
                            '[' + leftHit.target + ':' + str(rightHit.hitEnd),
                            str(rightHit.hitStart),
                        ]
                    )
                    + ']'
                )

                eleSeqRight.description = annotate(
                    source, leftHit.target, right_coord, genome_descriptions
                )

                # File routing follows the same role labels.
                for rec, label in (
                    (eleSeqLeft, left_label),
                    (eleSeqRight, right_label),
                ):
                    if label == 'L':
                        left_seqList.append(rec)
                    else:
                        right_seqList.append(rec)

            # Write separate FASTA files for left and right terminus hits
            if left_seqList:
                left_outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_left_term_hits_'
                    + str(len(left_seqList))
                    + '.fasta',
                )
                with open(left_outfile, 'w') as handle:
                    for seq in left_seqList:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')

            if right_seqList:
                right_outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_right_term_hits_'
                    + str(len(right_seqList))
                    + '.fasta',
                )
                with open(right_outfile, 'w') as handle:
                    for seq in right_seqList:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')
