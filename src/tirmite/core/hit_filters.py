"""Anchor filtering: keep only hits that reach the outer edge of their model.

A terminus HMM models one end of a transposable element. A genuine terminus
hit should align to the *outer* edge of that model -- the end that faces away
from the element body -- because that edge corresponds to the element/genome
junction. A hit that covers only the model's interior is usually an internal
repeat or a spurious match, and ``--max-offset`` exists to discard it.

Notes
-----
This module exists because ``compute_outer_edge_offset`` and
``filter_hits_by_anchor`` were implemented twice, independently, in
:mod:`tirmite.cli.hmm_pair` and :mod:`tirmite.cli.ensemble_search`. The
``hmm_pair`` copy received two corrections that the ``ensemble_search`` copy
never did:

1. The reverse-insertion strand swap shipped in 1.5.0. Without it,
   ``tirmite search --max-offset --pairing-map`` measured the offset against
   each hit's *inner* edge whenever an element was inserted in reverse, and
   silently discarded valid hits.
2. Symmetric-model handling. A pairing map row naming the same feature on both
   sides describes a symmetric element, whose model has no fixed terminus
   role. The ``ensemble_search`` copy labelled such a model ``'left'`` and then
   immediately overwrote it with ``'right'``.

Both subcommands now call this single implementation.
"""

import logging
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

import pandas as pd  # type: ignore[import-untyped]

# A pairing map reaches this module in two shapes: ``tirmite pair`` carries an
# ordered list of (left, right) tuples, while ``tirmite search`` carries a dict
# keyed on the left feature. Both are accepted and normalised on entry.
PairingMap = Union[Mapping[str, str], Sequence[Tuple[str, str]]]


class AnchorFilterError(Exception):
    """Raised when the anchor filter cannot be applied to the given inputs."""


def _normalise_pairing_map(pairing_map: Optional[PairingMap]) -> List[Tuple[str, str]]:
    """
    Reduce either accepted pairing-map shape to a list of (left, right) pairs.

    Parameters
    ----------
    pairing_map : mapping or sequence of tuple, optional
        Either ``{left: right}`` or ``[(left, right), ...]``.

    Returns
    -------
    list of tuple of str
        Pairs in a single canonical form. Empty if ``pairing_map`` is None.
    """
    if not pairing_map:
        return []
    if isinstance(pairing_map, Mapping):
        return list(pairing_map.items())
    return [tuple(pair) for pair in pairing_map]  # type: ignore[misc]


def compute_outer_edge_offset(
    hmm_start: int,
    hmm_end: int,
    model_len: int,
    strand: str,
    terminus_type: str,
) -> int:
    """
    Compute the offset between a hit alignment boundary and the outer edge of the query model.

    The "outer edge" is the external end of the terminus model:

    - For a left terminus the outer edge faces upstream (model position 1 on + strand).
    - For a right terminus the outer edge faces downstream (model position ``model_len`` on + strand).

    Parameters
    ----------
    hmm_start : int
        1-based start position of the alignment on the query/HMM model.
    hmm_end : int
        1-based end position of the alignment on the query/HMM model.
    model_len : int
        Total length of the query/HMM model.
    strand : str
        Strand of the hit: '+' or '-'.
    terminus_type : str
        'left' or 'right' terminus. This is the genomic side the outer edge
        faces, which is not necessarily the model's pairing role -- see
        :func:`filter_hits_by_anchor` for the distinction.

    Returns
    -------
    int
        Number of unaligned model positions between the hit and the outer edge.
        Zero means the hit reaches the outer edge exactly.
    """
    if terminus_type == 'left':
        if strand == '+':
            return hmm_start - 1
        else:  # '-'
            return model_len - hmm_end
    else:  # 'right'
        if strand == '+':
            return model_len - hmm_end
        else:  # '-'
            return hmm_start - 1


def filter_hits_by_anchor(
    hit_table: pd.DataFrame,
    model_lengths: Dict[str, int],
    max_offset: int,
    orientation: str = 'F,R',
    pairing_map: Optional[PairingMap] = None,
    on_missing_length: str = 'warn',
) -> pd.DataFrame:
    """
    Filter hits to only those anchored near the outer edge of the query model.

    For asymmetric pairings (different left/right models) or symmetric pairings
    with different strands (e.g. F,R), the terminus type is determined and only
    the external edge offset is checked.

    For symmetric same-strand pairings (F,F or R,R) without a pairing map, the
    hit must be within ``max_offset`` bases of **both** ends of the query model
    (that is, both ``hmmStart - 1 <= max_offset`` and
    ``model_len - hmmEnd <= max_offset``), since neither end can be identified
    as the outer one.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with columns: model, strand, hmmStart, hmmEnd.
    model_lengths : dict
        Mapping of model name to model length in positions.
    max_offset : int
        Maximum allowed offset from the outer edge, in model positions.
    orientation : str, default 'F,R'
        Comma-separated orientation codes (F=Forward/+, R=Reverse/-).
    pairing_map : mapping or sequence of tuple, optional
        Left/right feature pairs. Accepts either ``{left: right}`` or
        ``[(left, right), ...]``. When provided, terminus role is taken from
        the model name rather than inferred from strand.
    on_missing_length : {'warn', 'raise'}, default 'warn'
        What to do when a hit's model has no entry in ``model_lengths``.
        ``'warn'`` logs and keeps the hit unchecked; ``'raise'`` collects all
        such models and raises :class:`AnchorFilterError` at the end.

    Returns
    -------
    pandas.DataFrame
        Filtered hit table.

    Raises
    ------
    AnchorFilterError
        If ``on_missing_length='raise'`` and any hit's model length is absent.

    Notes
    -----
    Role versus side. For an asymmetric pair, the model name tells you the
    terminus *role* (this is the left-terminus model), but
    :func:`compute_outer_edge_offset` needs to know which genomic *side* the
    outer edge faces. Those agree only for a forward insertion. When a hit's
    strand is not the one the orientation expects for its role, the element is
    inserted in reverse and the two sides swap. Omitting that swap measures the
    offset against the hit's inner edge and discards valid reverse-inserted
    hits -- the bug fixed for ``tirmite pair`` in 1.5.0.
    """
    if hit_table.empty:
        return hit_table

    if on_missing_length not in ('warn', 'raise'):
        raise ValueError(
            f"on_missing_length must be 'warn' or 'raise', got {on_missing_length!r}"
        )

    # Imported lazily so that this module does not execute the pairing engine
    # at import time, which keeps tirmite.core free of import-order surprises.
    from tirmite.tirmitetools import parse_orientation

    # Parse orientation with the same validator PairingConfig uses, so the two
    # can never disagree about what a given orientation string means.
    try:
        orientation_parts = parse_orientation(orientation)
    except ValueError as e:
        logging.warning(f'{e} Skipping anchor filter.')
        return hit_table

    left_strand = '+' if orientation_parts[0] == 'F' else '-'
    right_strand = '+' if orientation_parts[1] == 'F' else '-'
    strands_differ = left_strand != right_strand

    # Build model-to-terminus map from the pairing map. A row that names the
    # same feature on both sides describes a symmetric element; such models
    # have no fixed terminus role, so they must fall through to the
    # strand-based or both-ends test rather than being labelled.
    model_terminus: Dict[str, str] = {}
    symmetric_models: Set[str] = set()
    for left_feature, right_feature in _normalise_pairing_map(pairing_map):
        if left_feature == right_feature:
            symmetric_models.add(left_feature)
            continue
        model_terminus[left_feature] = 'left'
        model_terminus[right_feature] = 'right'

    # A model listed symmetrically anywhere is symmetric everywhere.
    for model_name in symmetric_models:
        model_terminus.pop(model_name, None)

    kept: List[bool] = []
    skipped_no_terminus = 0
    removed = 0
    removed_per_model: Dict[str, int] = {}
    missing_lengths: Set[str] = set()

    for _, row in hit_table.iterrows():
        model = row['model']
        strand = row['strand']

        model_len = model_lengths.get(model)
        if model_len is None:
            missing_lengths.add(model)
            if on_missing_length == 'warn':
                logging.warning(
                    f'Model length not found for {model}, '
                    'keeping hit without anchor check'
                )
            kept.append(True)
            continue

        try:
            hmm_start = int(row['hmmStart'])
            hmm_end = int(row['hmmEnd'])
        except (ValueError, TypeError):
            # A hit with unparseable model coordinates cannot be checked; keep
            # it rather than silently dropping data on a formatting problem.
            kept.append(True)
            continue

        if model in model_terminus:
            # Asymmetric: translate pairing ROLE into the genomic SIDE the
            # outer edge faces. See the Notes section above.
            role = model_terminus[model]
            expected_strand = left_strand if role == 'left' else right_strand
            forward_insertion = strand == expected_strand
            if forward_insertion:
                terminus_type: Optional[str] = role
            else:
                terminus_type = 'right' if role == 'left' else 'left'
        elif strands_differ:
            # Symmetric with different strands: strand alone distinguishes the
            # two termini.
            if strand == left_strand:
                terminus_type = 'left'
            elif strand == right_strand:
                terminus_type = 'right'
            else:
                terminus_type = None
        else:
            # Same-strand symmetric (F,F or R,R) without a pairing map: neither
            # end can be identified as outer, so require both to be anchored.
            offset_start = hmm_start - 1
            offset_end = model_len - hmm_end
            if offset_start <= max_offset and offset_end <= max_offset:
                kept.append(True)
            else:
                kept.append(False)
                removed += 1
                removed_per_model[model] = removed_per_model.get(model, 0) + 1
            continue

        if terminus_type is None:
            skipped_no_terminus += 1
            kept.append(True)
            continue

        offset = compute_outer_edge_offset(
            hmm_start, hmm_end, model_len, strand, terminus_type
        )

        if offset <= max_offset:
            kept.append(True)
        else:
            kept.append(False)
            removed += 1
            removed_per_model[model] = removed_per_model.get(model, 0) + 1

    if missing_lengths and on_missing_length == 'raise':
        missing_list = ', '.join(sorted(missing_lengths))
        raise AnchorFilterError(
            f'Anchor filter requires model lengths for {len(missing_lengths)} '
            f'model(s) that are not available: {missing_list}. '
            'Provide lengths via --fasta, --hmm, or --lengths-file.'
        )

    if skipped_no_terminus:
        logging.debug(
            f'Anchor filter: {skipped_no_terminus} hit(s) kept without checking - '
            'terminus type could not be determined. For same-strand symmetric '
            'pairing, supply --pairing-map to identify left/right models.'
        )

    result = hit_table[kept].copy()

    logging.info(
        f'Anchor filter (max_offset={max_offset}): '
        f'{len(hit_table)} -> {len(result)} hits ({removed} removed)'
    )
    if removed_per_model:
        for model_name, count in sorted(removed_per_model.items()):
            logging.info(f'  {model_name}: {count} hit(s) excluded by anchor filter')
    return result


def expand_pairing_map_to_components(
    pairing_map: Optional[PairingMap],
    cluster_map: Optional[Mapping[str, Iterable[str]]],
) -> List[Tuple[str, str]]:
    """
    Rewrite a cluster-level pairing map in terms of component model names.

    Parameters
    ----------
    pairing_map : mapping or sequence of tuple, optional
        Pairs of features, which may name clusters rather than models.
    cluster_map : mapping of str to iterable of str, optional
        Mapping of cluster name to its component model names.

    Returns
    -------
    list of tuple of str
        Pairs expressed entirely in component model names. A feature with no
        entry in ``cluster_map`` passes through unchanged, so maps that mix
        cluster names and bare model names work.

    Notes
    -----
    ``tirmite search`` applies the anchor filter *before*
    ``merge_overlapping_cluster_hits`` renames hit models to cluster names,
    because anchor offsets are measured against per-model query lengths, which
    only exist at component level. A cluster-level pairing map therefore
    matches nothing at that point in the pipeline unless it is expanded here.

    A self-paired cluster (``LEFT_TIR`` paired with ``LEFT_TIR``) expands to
    every component against every component, including each against itself.
    The resulting ``left == right`` rows are exactly what
    :func:`filter_hits_by_anchor` treats as symmetric, which is correct.

    Examples
    --------
    >>> expand_pairing_map_to_components(
    ...     {'CL': 'CR'}, {'CL': ['a1', 'a2'], 'CR': ['b1']}
    ... )
    [('a1', 'b1'), ('a2', 'b1')]
    """
    pairs = _normalise_pairing_map(pairing_map)
    if not pairs:
        return []
    if not cluster_map:
        return pairs

    def components(feature: str) -> List[str]:
        """
        Return the component models of a cluster, or the feature itself.

        Parameters
        ----------
        feature : str
            A cluster name or a bare model name.

        Returns
        -------
        list of str
            Component model names, or ``[feature]`` if it names no cluster.
        """
        members = cluster_map.get(feature)
        if not members:
            return [feature]
        return list(members)

    expanded: List[Tuple[str, str]] = []
    seen: Set[Tuple[str, str]] = set()
    for left_feature, right_feature in pairs:
        for left_component in components(left_feature):
            for right_component in components(right_feature):
                pair = (left_component, right_component)
                if pair not in seen:
                    seen.add(pair)
                    expanded.append(pair)
    return expanded
