"""
Summary statistics for a single terminus pairing.

The counting and the formatting live in separate functions on purpose. The
plain-text ``*_summary.txt`` file written by ``tirmite pair`` and the HTML
report both need the same numbers; deriving them twice is how the two outputs
end up disagreeing. :func:`pair_summary_stats` produces the numbers,
:func:`format_pair_summary` renders the text file, and the HTML report consumes
the :class:`PairSummary` directly.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Set

__all__ = [
    'PairSummary',
    'format_pair_summary',
    'pair_summary_stats',
]


@dataclass(frozen=True)
class PairSummary:
    """
    Counts describing the outcome of one terminus pairing procedure.

    Attributes
    ----------
    left_feature : str
        Name of the left model/feature in the pairing.
    right_feature : str
        Name of the right model/feature in the pairing. Equal to
        `left_feature` for symmetric pairings.
    hits_per_model : dict of str to int
        Number of hits retained for each model involved in this pairing,
        after all filters.
    total_hits : int
        Sum of `hits_per_model` values.
    total_pairs : int
        Number of terminus pairs found.
    total_elements : int
        Number of elements successfully extracted from those pairs. May be
        lower than `total_pairs` when extraction fails or is disabled.
    total_unpaired : int
        Number of retained hits that were not assigned a partner.
    filter_stats : dict
        Raw hit-filtering statistics as accumulated by the CLI. Keys used
        when rendering: 'initial_hits', 'mincov', 'coverage_excluded',
        'maxeval', 'evalue_excluded', 'max_offset', 'anchor_excluded',
        'after_filtering', 'pairing_map_models_ignored',
        'pairing_map_hits_ignored'.
    """

    left_feature: str
    right_feature: str
    hits_per_model: Dict[str, int]
    total_hits: int
    total_pairs: int
    total_elements: int
    total_unpaired: int
    filter_stats: Dict[str, Any] = field(default_factory=dict)

    @property
    def pair_label(self) -> str:
        """
        Return the filesystem/UI label for this pairing.

        Returns
        -------
        str
            ``'<left>_<right>'`` for an asymmetric pairing, or just the model
            name when the pairing is symmetric.
        """
        if self.left_feature != self.right_feature:
            return f'{self.left_feature}_{self.right_feature}'
        return self.left_feature

    @property
    def group_id(self) -> str:
        """
        Return a globally unique identifier for this pairing.

        Returns
        -------
        str
            ``'<left>__<right>'``. Unlike :attr:`pair_label` this does not
            collapse for symmetric pairings, so it stays unique when a run
            mixes symmetric and asymmetric rows from a pairing map.
        """
        return f'{self.left_feature}__{self.right_feature}'


def pair_summary_stats(
    left_feature: str,
    right_feature: str,
    pair_paired: Dict[str, List[Any]],
    pair_hitIndex: Dict[str, Dict[int, Dict[str, Any]]],
    total_pairs: int,
    total_elements: int,
    filter_stats: Dict[str, Any] | None = None,
) -> PairSummary:
    """
    Count the outcome of one terminus pairing procedure.

    Parameters
    ----------
    left_feature : str
        Name of the left model/feature.
    right_feature : str
        Name of the right model/feature.
    pair_paired : dict
        Paired hit indices for this pairing, keyed by model:
        ``paired[model] = [{idxA, idxB}, ...]``.
    pair_hitIndex : dict
        Hit index for this pairing: ``hitIndex[model][idx] = {...}``.
    total_pairs : int
        Number of pairs found.
    total_elements : int
        Number of elements extracted.
    filter_stats : dict, optional
        Hit-filtering statistics to carry through to the report.

    Returns
    -------
    PairSummary
        The counted summary.
    """
    # Both sides of the pairing are counted even when one contributed no hits,
    # so a zero row is visible in the report rather than silently absent.
    models = {left_feature, right_feature}
    hits_per_model: Dict[str, int] = {
        model: len(pair_hitIndex[model]) if model in pair_hitIndex else 0
        for model in models
    }

    # A hit index may appear in more than one pair set only if pairing is
    # buggy, but counting through a set keeps this robust either way.
    paired_hit_ids: Set[int] = set()
    for pairs in pair_paired.values():
        for pair_set in pairs:
            paired_hit_ids.update(pair_set)

    total_hits = sum(hits_per_model.values())

    return PairSummary(
        left_feature=left_feature,
        right_feature=right_feature,
        hits_per_model=hits_per_model,
        total_hits=total_hits,
        total_pairs=total_pairs,
        total_elements=total_elements,
        total_unpaired=total_hits - len(paired_hit_ids),
        filter_stats=dict(filter_stats) if filter_stats else {},
    )


def format_pair_summary(summary: PairSummary) -> str:
    """
    Render a :class:`PairSummary` as the plain-text summary report.

    Parameters
    ----------
    summary : PairSummary
        The counted summary to render.

    Returns
    -------
    str
        The full text of the ``*_summary.txt`` file, including its trailing
        newline.
    """
    lines: List[str] = []
    lines.append('TIRmite Pair Summary Report\n')
    lines.append('===========================\n\n')
    lines.append(f'Model pair: {summary.left_feature} <-> {summary.right_feature}\n\n')

    filter_stats = summary.filter_stats
    if filter_stats:
        lines.append('Filtering criteria applied\n')
        lines.append('--------------------------\n')

        initial = filter_stats.get('initial_hits')
        if initial is not None:
            lines.append(f'  Initial hits imported: {initial}\n')

        # Pairing-map pre-filter: models absent from the map never reach pairing.
        ignored_models = filter_stats.get('pairing_map_models_ignored')
        if ignored_models:
            lines.append(
                f'  Pairing-map model filter: {len(ignored_models)} model(s) ignored '
                f'({", ".join(sorted(ignored_models))}), '
                f'{filter_stats.get("pairing_map_hits_ignored", 0)} hits excluded\n'
            )

        mincov = filter_stats.get('mincov')
        if mincov is not None:
            cov_excluded = filter_stats.get('coverage_excluded', 0)
            lines.append(
                f'  Coverage filter (min coverage >= {mincov}): '
                f'{cov_excluded} hit(s) excluded\n'
            )

        maxeval = filter_stats.get('maxeval')
        if maxeval is not None:
            eval_excluded = filter_stats.get('evalue_excluded', 0)
            lines.append(
                f'  E-value filter (max e-value <= {maxeval}): '
                f'{eval_excluded} hit(s) excluded\n'
            )

        max_offset = filter_stats.get('max_offset')
        if max_offset is not None:
            anchor_excluded = filter_stats.get('anchor_excluded')
            excl_str = (
                str(anchor_excluded) if anchor_excluded is not None else 'unknown'
            )
            lines.append(
                f'  Anchor offset filter (max offset <= {max_offset}): '
                f'{excl_str} hit(s) excluded\n'
            )

        after_filtering = filter_stats.get('after_filtering')
        if after_filtering is not None:
            lines.append(f'  Hits remaining after all filters: {after_filtering}\n')
        lines.append('\n')

    lines.append('Hits per model (after all filters)\n')
    lines.append('----------------------------------\n')
    for model, count in sorted(summary.hits_per_model.items()):
        lines.append(f'  {model}: {count}\n')
    lines.append(f'\nTotal hits for this pair: {summary.total_hits}\n')
    lines.append(f'Total pairs found: {summary.total_pairs}\n')
    lines.append(f'Total paired elements extracted: {summary.total_elements}\n')
    lines.append(f'Total unpaired hits: {summary.total_unpaired}\n')

    return ''.join(lines)
