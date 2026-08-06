"""
Stacked alignments of every hit to a terminus model.

Each panel answers one question: across the whole genome, what does this
terminus actually look like, and where do the hits agree? Rows are the hits
themselves, in model orientation, so a minus-strand hit is reverse-complemented
before it is shown.

Two kinds of missing sequence are tracked separately, because they mean
different things:

``m`` (model pad)
    The hit did not match the model at this position, but the genome does
    extend there. Drawn grey: sequence exists, the alignment simply did not
    claim it.
``g`` (gap)
    The contig ended, so this position does not exist at all. Drawn blank, and
    the sequence carries a gap character.

Conflating the two would report a partial match where the truth is a truncated
assembly.

MAFFT is used when it is on ``PATH``. When it is not, or when it fails, rows
are placed by their model coordinates instead. That fallback is exact wherever
the alignment coordinates exist, which is every input format except BED.
"""

import logging
import os
from typing import Any, Dict, List, Optional, Sequence, Tuple

from tirmite.report.model import MsaPanel, MsaRow

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['ANCHOR', 'LEFT', 'MAFFT', 'build_msa_panels', 'pad_runs']

# Aligner labels as they appear in the report.
MAFFT = 'mafft'
ANCHOR = 'anchor'
LEFT = 'left'

GAP = '-'
# Kinds of padded column, matching the module docstring.
PAD_MODEL = 'm'
PAD_GAP = 'g'

_COMPLEMENT = str.maketrans(
    'ACGTUNRYSWKMBVDHacgtunryswkmbvdh-.',
    'TGCAANYRSWMKVBHDtgcaanyrswmkvbhd-.',
)


def _reverse_complement(seq: str) -> str:
    """
    Reverse-complement a nucleotide string.

    Parameters
    ----------
    seq : str
        Nucleotide sequence, which may contain ambiguity codes and gaps.

    Returns
    -------
    str
        The reverse complement.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def pad_runs(seq: str, model_pad: Sequence[Tuple[int, int]]) -> List[List[Any]]:
    """
    Compress padded positions into runs for the report payload.

    Parameters
    ----------
    seq : str
        The row's aligned sequence.
    model_pad : sequence of (int, int)
        Half-open ``(start, stop)`` column ranges where the model was not
        matched although the genome continues.

    Returns
    -------
    list of list
        Runs as ``[start_column, length, kind]``, where kind is 'm' for a model
        pad and 'g' for a position that does not exist.

    Notes
    -----
    Run-length encoding matters: a panel can be thousands of rows by hundreds
    of columns, and emitting one entry per padded cell would dominate the size
    of the whole report.
    """
    runs: List[List[Any]] = []

    for start, stop in model_pad:
        if stop > start:
            runs.append([start, stop - start, PAD_MODEL])

    # Gap characters are found in the sequence itself rather than passed in:
    # MAFFT introduces gaps of its own, and they read the same way.
    run_start = None
    for i, char in enumerate(seq):
        if char == GAP:
            if run_start is None:
                run_start = i
        elif run_start is not None:
            runs.append([run_start, i - run_start, PAD_GAP])
            run_start = None
    if run_start is not None:
        runs.append([run_start, len(seq) - run_start, PAD_GAP])

    runs.sort(key=lambda run: (run[0], run[2]))
    return runs


def _fetch_row(
    source: Any,
    hit: Any,
    model_length: Optional[int],
    pad_model: bool = False,
) -> Optional[Tuple[str, int, int, int, int]]:
    """
    Fetch one hit's sequence, optionally extended to the width of its model.

    Parameters
    ----------
    source : SequenceSource
        Sequence source to read from.
    hit : HitRecord
        The hit to fetch.
    model_length : int or None
        Declared model length, used only for logging.
    pad_model : bool, default False
        When True, extend the window past the hit to cover the model positions
        the alignment did not claim, and show that sequence greyed. When
        False, fetch only what the hit matched and leave the rest as gaps.

    Returns
    -------
    tuple or None
        ``(sequence, left_pad, right_pad, left_model, right_model)`` in model
        orientation, where the pads count positions that do not exist and the
        model counts are unmatched model positions carrying real sequence.
        None when the sequence could not be read.

    Notes
    -----
    Either way every row spans the same model interval, so the rows line up
    without any alignment; the difference is only whether the unclaimed
    columns show sequence or gaps.

    Padding is off by default because the extra bases are not evidence for the
    model: they are whatever happens to sit beside the hit, and showing them
    in an alignment invites reading them as part of the match.
    :func:`tirmite.utils.extract.fetch_region_padded` reports exactly how much
    of a requested window fell off the contig, which is what separates a gap
    from a model pad when padding is on.
    """
    from tirmite.utils.extract import fetch_region_padded

    left_deficit = hit.trunc_left or 0
    right_deficit = hit.trunc_right or 0

    if not pad_model:
        # Only the matched region is read; the deficit becomes gap columns, so
        # no sequence is shown that the model did not claim.
        region = fetch_region_padded(
            source, hit.contig, hit.start, hit.end, pad_char=GAP
        )
        if region is None:
            logger.debug(
                f'No sequence available for {hit.contig}:{hit.start}-{hit.end} '
                f'({hit.model}, model length {model_length})'
            )
            return None
        seq = GAP * left_deficit + region.seq + GAP * right_deficit
        left_pad = left_deficit + region.left_pad
        right_pad = right_deficit + region.right_pad
        if hit.strand == '-':
            seq = _reverse_complement(seq)
            left_pad, right_pad = right_pad, left_pad
        return seq, left_pad, right_pad, 0, 0

    start = hit.start - left_deficit
    end = hit.end + right_deficit

    region = fetch_region_padded(source, hit.contig, start, end, pad_char=GAP)
    if region is None:
        logger.debug(
            f'No sequence available for {hit.contig}:{start}-{end} '
            f'({hit.model}, model length {model_length})'
        )
        return None

    seq = region.seq
    left_pad, right_pad = region.left_pad, region.right_pad

    # Whatever of the deficit was not lost to the contig end is a model pad:
    # real sequence that the alignment did not claim.
    left_model = max(0, left_deficit - left_pad)
    right_model = max(0, right_deficit - right_pad)

    if hit.strand == '-':
        # Rows are shown in model orientation, so both the sequence and the
        # per-end bookkeeping flip.
        seq = _reverse_complement(seq)
        left_pad, right_pad = right_pad, left_pad
        left_model, right_model = right_model, left_model

    return seq, left_pad, right_pad, left_model, right_model


def _anchor_rows(
    rows: List[Dict[str, Any]],
) -> Tuple[List[str], List[List[Tuple[int, int]]], int]:
    """
    Place rows into a common frame using their model coordinates.

    Parameters
    ----------
    rows : list of dict
        Row records carrying 'seq', 'left_pad' and 'right_pad'.

    Returns
    -------
    seqs : list of str
        Sequences padded to a common width.
    model_pad : list of list of (int, int)
        Per-row column ranges where the model was not matched.
    width : int
        Panel width in columns.

    Notes
    -----
    Every row was fetched over its own full model interval, so the rows already
    describe the same model positions and only need padding to a common width.
    Rows of unequal length after that mean insertions relative to the model,
    which this fallback cannot resolve -- they are right-padded, and the panel
    is labelled so the reader knows an aligner would do better.
    """
    width = max((len(row['seq']) for row in rows), default=0)
    seqs = []
    model_pad = []
    for row in rows:
        seq = row['seq']
        trailing = width - len(seq)
        seqs.append(seq + GAP * trailing)
        # The model pad sits *inside* any contig-end padding: the outermost
        # columns are positions that do not exist, and the model pad is the
        # real sequence just inside them.
        pads = []
        if row['left_model']:
            start = row['left_pad']
            pads.append((start, start + row['left_model']))
        if row['right_model']:
            stop = len(seq) - row['right_pad']
            pads.append((stop - row['right_model'], stop))
        model_pad.append(pads)
    return seqs, model_pad, width


def _mafft_rows(
    rows: List[Dict[str, Any]],
    tmpdir: str,
    model: str,
) -> Optional[Tuple[List[str], List[List[Tuple[int, int]]], int]]:
    """
    Align rows with MAFFT and map the model pads onto the alignment columns.

    Parameters
    ----------
    rows : list of dict
        Row records carrying 'seq', 'left_model' and 'right_model'.
    tmpdir : str
        Directory for MAFFT's intermediate files.
    model : str
        Model name, for logging.

    Returns
    -------
    tuple or None
        ``(sequences, model_pad, width)`` as for :func:`_anchor_rows`, or None
        if MAFFT was unavailable or failed.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from tirmite.runners.mafft import align_in_memory, mafft_available

    if not mafft_available():
        return None
    if len(rows) < 2:
        # MAFFT needs at least two sequences; a single row is already aligned.
        return None

    try:
        os.makedirs(tmpdir, exist_ok=True)
    except OSError as exc:
        # An unusable temporary directory costs the alignment, not the report.
        logger.warning(f'Cannot write MAFFT temporary files for {model}: {exc}')
        return None

    records = [
        SeqRecord(Seq(row['seq'].replace(GAP, '')), id=f'row{i}', description='')
        for i, row in enumerate(rows)
    ]
    # Rows reduced to nothing by contig truncation cannot be aligned.
    if any(len(record.seq) == 0 for record in records):
        logger.debug(f'Skipping MAFFT for {model}: a row has no sequence')
        return None

    aligned = align_in_memory(records, tmpdir)
    if aligned is None:
        return None

    by_id = {record.id: str(record.seq) for record in aligned}
    if len(by_id) != len(records):
        logger.warning(
            f'MAFFT returned {len(by_id)} of {len(records)} rows for {model}; '
            'falling back to model-coordinate anchoring'
        )
        return None

    seqs = []
    model_pad = []
    for i, row in enumerate(rows):
        seq = by_id.get(f'row{i}')
        if seq is None:
            return None
        seqs.append(seq)
        # The model pads were leading and trailing runs of the ungapped
        # sequence; after alignment they occupy however many columns those
        # residues now span.
        model_pad.append(_project_pads(seq, row['left_model'], row['right_model']))

    width = max((len(seq) for seq in seqs), default=0)
    seqs = [seq + GAP * (width - len(seq)) for seq in seqs]
    return seqs, model_pad, width


def _project_pads(
    aligned: str, left_model: int, right_model: int
) -> List[Tuple[int, int]]:
    """
    Map leading and trailing model pads onto aligned column positions.

    Parameters
    ----------
    aligned : str
        The row's aligned sequence, with gap characters.
    left_model, right_model : int
        Number of leading and trailing residues that are model pad.

    Returns
    -------
    list of (int, int)
        Half-open column ranges covering those residues.
    """
    pads: List[Tuple[int, int]] = []
    if left_model:
        seen = 0
        for col, char in enumerate(aligned):
            if char != GAP:
                seen += 1
                if seen == left_model:
                    pads.append((0, col + 1))
                    break
    if right_model:
        seen = 0
        for col in range(len(aligned) - 1, -1, -1):
            if aligned[col] != GAP:
                seen += 1
                if seen == right_model:
                    pads.append((col, len(aligned)))
                    break
    return pads


def build_msa_panels(
    hits: Sequence[Any],
    groups: Sequence[Any],
    models: Sequence[Any],
    *,
    source: Any,
    tmpdir: str,
    mode: str = 'auto',
    max_rows: int = 500,
    max_cells: int = 2_000_000,
    pad_model: bool = False,
) -> List[MsaPanel]:
    """
    Build one stacked alignment panel per terminus model.

    Parameters
    ----------
    hits : sequence of HitRecord
        All hits in the report.
    groups : sequence of PairingGroup
        Pairing groups, used to colour rows by group.
    models : sequence of ModelInfo
        Models in the report.
    source : SequenceSource
        Sequence source to read hit sequences from.
    tmpdir : str
        Directory for MAFFT's intermediate files. Should be inside the run's
        own temporary directory so it is cleaned up with everything else.
    mode : {'auto', 'mafft', 'anchor', 'off'}, default 'auto'
        'auto' uses MAFFT when available and falls back to model-coordinate
        anchoring; 'mafft' and 'anchor' force one or the other; 'off' builds
        no panels.
    max_rows : int, default 500
        Row cap per panel, applied best-e-value-first.
    max_cells : int, default 2000000
        Cell budget per panel. Rows are dropped beyond it.
    pad_model : bool, default False
        Show the sequence beside a hit where its model went unmatched, greyed.
        Off by default: those bases are not evidence for the model, and an
        alignment invites reading them as part of the match.

    Returns
    -------
    list of MsaPanel
        One panel per model with at least one usable row.
    """
    if mode == 'off' or source is None:
        return []

    group_index = {group.group_id: i for i, group in enumerate(groups)}
    model_lengths = {model.name: model.length for model in models}

    by_model: Dict[str, List[Any]] = {}
    for hit in hits:
        by_model.setdefault(hit.model, []).append(hit)

    panels = []
    for name in sorted(by_model):
        panel = _build_panel(
            name,
            by_model[name],
            group_index,
            model_lengths.get(name),
            source=source,
            tmpdir=os.path.join(tmpdir, f'msa_{_safe_name(name)}'),
            mode=mode,
            max_rows=max_rows,
            max_cells=max_cells,
            pad_model=pad_model,
        )
        if panel is not None:
            panels.append(panel)
    return panels


def _safe_name(name: str) -> str:
    """
    Reduce a model name to something usable as a directory name.

    Parameters
    ----------
    name : str
        Model name, which comes from a user-supplied file.

    Returns
    -------
    str
        The name with anything outside a conservative safe set replaced.
    """
    return ''.join(c if c.isalnum() or c in '._-' else '_' for c in name)[:64]


def _build_panel(
    model: str,
    hits: Sequence[Any],
    group_index: Dict[str, int],
    model_length: Optional[int],
    *,
    source: Any,
    tmpdir: str,
    mode: str,
    max_rows: int,
    max_cells: int,
    pad_model: bool = False,
) -> Optional[MsaPanel]:
    """
    Build the alignment panel for one terminus model.

    Parameters
    ----------
    model : str
        Model name.
    hits : sequence of HitRecord
        Every hit for this model.
    group_index : dict
        Group id to index.
    model_length : int or None
        Declared model length.
    source : SequenceSource
        Sequence source.
    tmpdir : str
        Directory for MAFFT's intermediate files.
    mode : str
        Alignment mode, as for :func:`build_msa_panels`.
    max_rows : int
        Row cap.
    max_cells : int
        Cell budget.
    pad_model : bool, default False
        Show unmatched model positions as sequence rather than gaps.

    Returns
    -------
    MsaPanel or None
        The panel, or None if no row could be built.
    """
    total = len(hits)
    # Best hits first, so the cap keeps the most informative rows.
    ordered = sorted(
        hits, key=lambda h: (h.evalue if h.evalue is not None else float('inf'), h.uid)
    )[:max_rows]

    rows: List[Dict[str, Any]] = []
    for hit in ordered:
        fetched = _fetch_row(source, hit, model_length, pad_model=pad_model)
        if fetched is None:
            continue
        seq, left_pad, right_pad, left_model, right_model = fetched
        rows.append(
            {
                'uid': hit.uid,
                'group_i': next(
                    (group_index[g] for g in hit.group_ids if g in group_index), None
                ),
                'role': hit.role,
                'seq': seq,
                'left_pad': left_pad,
                'right_pad': right_pad,
                'left_model': left_model,
                'right_model': right_model,
                'has_model_coords': hit.hmm_start is not None,
            }
        )

    if not rows:
        return None

    note = None
    aligned = None
    if mode in ('auto', 'mafft'):
        aligned = _mafft_rows(rows, tmpdir, model)
        if aligned is not None:
            aligner = MAFFT
    if aligned is None:
        if mode == 'mafft':
            note = 'MAFFT was unavailable or failed; rows are placed by model position.'
        seqs, model_pad, width = _anchor_rows(rows)
        # Without model coordinates there is nothing to anchor on and the rows
        # are merely left-aligned, which the caption has to admit.
        if not any(row['has_model_coords'] for row in rows):
            aligner = LEFT
            note = (
                'This input format carries no model coordinates, so rows are '
                'left-aligned rather than anchored to the model.'
            )
        else:
            aligner = ANCHOR
            if note is None:
                note = (
                    'MAFFT was not available; rows are placed by their model '
                    'coordinates.'
                )
    else:
        seqs, model_pad, width = aligned

    # Drop rows beyond the cell budget rather than shipping a panel the browser
    # cannot draw. Rows are already ordered best-first.
    if width and len(rows) * width > max_cells:
        keep = max(1, max_cells // width)
        if keep < len(rows):
            rows = rows[:keep]
            seqs = seqs[:keep]
            model_pad = model_pad[:keep]

    panel_rows = [
        MsaRow(
            uid=row['uid'],
            group_i=row['group_i'],
            role=row['role'],
            seq=seqs[i],
            pad=pad_runs(seqs[i], model_pad[i]),
        )
        for i, row in enumerate(rows)
    ]

    return MsaPanel(
        model=model,
        aligner=aligner,
        n_cols=width,
        n_rows_shown=len(panel_rows),
        n_rows_total=total,
        rows=panel_rows,
        note=note,
    )
