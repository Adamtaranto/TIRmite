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
from statistics import median
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set

if TYPE_CHECKING:  # pragma: no cover - import cycle avoidance for annotations
    from tirmite.report.model import ReportData

__all__ = [
    'PairSummary',
    'contig_table',
    'filter_table',
    'format_pair_summary',
    'group_table',
    'model_table',
    'pair_summary_stats',
    'search_cluster_table',
    'search_contested_table',
    'search_cross_match_matrix',
    'search_cross_match_table',
    'search_model_overlap_table',
    'search_query_table',
    'search_stage_table',
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


# ---------------------------------------------------------------------------
# Report tables
#
# Each returns rows ready for the shared table template. A cell is either a
# plain value or a dict of {value, sort, colour}: `sort` carries the underlying
# number when the display text is formatted for reading and would sort wrongly
# as a string, and `colour` attaches a pairing-group swatch so identity is not
# carried by position alone.
# ---------------------------------------------------------------------------


def _bp(value: Optional[int]) -> Dict[str, Any]:
    """
    Format a base-pair count as a readable, correctly sortable cell.

    Parameters
    ----------
    value : int or None
        Length in base pairs.

    Returns
    -------
    dict
        Cell with a human-readable value and its numeric sort key.
    """
    if value is None:
        return {'value': '–', 'sort': None}
    if value >= 1_000_000:
        text = f'{value / 1_000_000:.2f} Mb'
    elif value >= 10_000:
        text = f'{value / 1000:.1f} kb'
    else:
        text = f'{value:,} bp'
    return {'value': text, 'sort': value}


def _pct(value: Optional[float]) -> Dict[str, Any]:
    """
    Format a fraction as a percentage cell.

    Parameters
    ----------
    value : float or None
        Fraction between 0 and 1.

    Returns
    -------
    dict
        Cell with a percentage string and its numeric sort key.
    """
    if value is None:
        return {'value': '–', 'sort': None}
    return {'value': f'{value * 100:.0f}%', 'sort': value}


def group_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the per-pairing-group statistics rows.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per pairing group.
    """
    lengths_by_group: Dict[int, List[int]] = {}
    for element in data.elements:
        lengths_by_group.setdefault(element.group_i, []).append(element.length)

    rows = []
    for i, group in enumerate(data.groups):
        lengths = lengths_by_group.get(i, [])
        rows.append(
            {
                'group': {
                    'value': group.label,
                    'sort': group.label,
                    # A custom property rather than the hex, so the swatch
                    # follows the surface into dark mode.
                    'colour': f'var(--group-{i})',
                },
                'orientation': group.orientation or '–',
                'pairs': {'value': f'{group.n_pairs:,}', 'sort': group.n_pairs},
                'elements': {
                    'value': f'{group.n_elements:,}',
                    'sort': group.n_elements,
                },
                'unpaired': {
                    'value': f'{group.n_unpaired:,}',
                    'sort': group.n_unpaired,
                },
                'median_length': _bp(int(median(lengths)) if lengths else None),
            }
        )
    return rows


def model_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the per-model statistics rows.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per terminus model.
    """
    return [
        {
            'model': model.name,
            'length': {
                'value': f'{model.length:,}' if model.length else '–',
                'sort': model.length,
            },
            'hits': {'value': f'{model.n_hits:,}', 'sort': model.n_hits},
            'contigs': {'value': f'{model.n_contigs:,}', 'sort': model.n_contigs},
            'paired': {'value': f'{model.n_paired:,}', 'sort': model.n_paired},
            'coverage': _pct(model.median_model_coverage),
            'full': _pct(model.frac_full_length),
            'clipped': _pct(model.frac_clipped),
        }
        for model in data.models
    ]


def contig_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the per-sequence statistics rows.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per sequence carrying at least one hit, ordered by descending
        pair count so the sequences worth looking at appear first.
    """
    ordered = sorted(data.contigs, key=lambda c: (-c.n_pairs, -c.n_hits, c.name))
    rows: List[Dict[str, Any]] = []
    for contig in ordered:
        per_mb = contig.n_elements / (contig.length / 1_000_000) if contig.length else 0
        rows.append(
            {
                'contig': contig.name
                + (' (estimated length)' if contig.length_source == 'inferred' else ''),
                'length': _bp(contig.length),
                'hits': {'value': f'{contig.n_hits:,}', 'sort': contig.n_hits},
                'pairs': {'value': f'{contig.n_pairs:,}', 'sort': contig.n_pairs},
                'element_bp': _bp(contig.element_bp),
                'density': {'value': f'{per_mb:.2f}', 'sort': per_mb},
            }
        )
    return rows


def filter_table(filter_stats: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Build the filtering statistics rows.

    Parameters
    ----------
    filter_stats : dict
        Statistics as accumulated by the CLI.

    Returns
    -------
    list of dict
        One row per filter that was actually applied. A filter that was not
        configured is omitted rather than shown as zero, which would imply it
        ran and excluded nothing.
    """
    if not filter_stats:
        return []

    rows: List[Dict[str, Any]] = []

    initial = filter_stats.get('initial_hits')
    if initial is not None:
        rows.append(
            {
                'filter': 'Hits imported',
                'setting': '–',
                'excluded': {'value': f'{initial:,}', 'sort': initial},
            }
        )

    ignored = filter_stats.get('pairing_map_models_ignored')
    if ignored:
        excluded = filter_stats.get('pairing_map_hits_ignored', 0)
        rows.append(
            {
                'filter': 'Pairing-map model filter',
                'setting': f'{len(ignored)} model(s) not in map: '
                + ', '.join(sorted(ignored)),
                'excluded': {'value': f'{excluded:,}', 'sort': excluded},
            }
        )

    for key, label, setting_key in (
        ('coverage_excluded', 'Coverage filter', 'mincov'),
        ('evalue_excluded', 'E-value filter', 'maxeval'),
        ('anchor_excluded', 'Anchor offset filter', 'max_offset'),
    ):
        setting = filter_stats.get(setting_key)
        if setting is None:
            continue
        excluded = filter_stats.get(key)
        rows.append(
            {
                'filter': label,
                'setting': f'{setting_key} = {setting}',
                'excluded': {
                    'value': f'{excluded:,}' if excluded is not None else 'unknown',
                    'sort': excluded,
                },
            }
        )

    remaining = filter_stats.get('after_filtering')
    if remaining is not None:
        rows.append(
            {
                'filter': 'Hits remaining',
                'setting': '–',
                'excluded': {'value': f'{remaining:,}', 'sort': remaining},
            }
        )

    return rows


# ---------------------------------------------------------------------------
# Search report tables
#
# These read the `stats` payload the search accumulator assembles, rather than
# the run's internals, so they stay usable from a saved report.
# ---------------------------------------------------------------------------


def search_query_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the per-query-terminus statistics rows.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per query model that produced a hit.
    """
    hits_before = data.stats.get('hits_per_model_before_merge', {})
    cluster_map = data.stats.get('cluster_map', {})
    rows = []
    for model in data.models:
        # After merging a row is a cluster, while the recorded counts are per
        # component, so a cluster's figure is the sum of what it absorbed.
        if model.name in cluster_map:
            components = cluster_map[model.name]
            contributions = [hits_before[c] for c in components if c in hits_before]
            before = sum(contributions) if contributions else None
        else:
            before = hits_before.get(model.name)
        rows.append(
            {
                'model': model.name,
                'length': {
                    'value': f'{model.length:,}' if model.length else '–',
                    'sort': model.length,
                },
                'hits': {'value': f'{model.n_hits:,}', 'sort': model.n_hits},
                'before_merge': {
                    'value': f'{before:,}' if before is not None else '–',
                    'sort': before,
                },
                'contigs': {'value': f'{model.n_contigs:,}', 'sort': model.n_contigs},
                'coverage': _pct(model.median_model_coverage),
                'full': _pct(model.frac_full_length),
                'clipped': _pct(model.frac_clipped),
            }
        )
    return rows


def search_stage_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the pipeline-stage rows, showing where hits were lost.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per stage, with the change from the previous stage.

    Notes
    -----
    The per-stage delta is the point of this table. The log records absolute
    counts at each step, which leaves the reader subtracting to find the step
    that removed everything.
    """
    stages = data.stats.get('stages', [])
    rows = []
    previous = None
    for stage in stages:
        count = stage['hits']
        if previous is None:
            change = {'value': '–', 'sort': None}
        else:
            delta = count - previous
            change = {
                'value': f'{delta:+,}' if delta else '0',
                'sort': delta,
            }
        rows.append(
            {
                'stage': stage['label'],
                'hits': {'value': f'{count:,}', 'sort': count},
                'change': change,
            }
        )
        previous = count
    return rows


def search_cross_match_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the cross-match rows: loci where two groups both hit.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per contested locus, longest overlap first.
    """
    matches = data.stats.get('cross_matches', [])
    rows = []
    for match in sorted(matches, key=lambda m: -m.get('overlap', 0)):
        rows.append(
            {
                'contig': match['contig'],
                'a': match['a_cluster'],
                'a_span': {
                    'value': f'{match["a_start"]:,}–{match["a_end"]:,}',
                    'sort': match['a_start'],
                },
                'b': match['b_cluster'],
                'b_span': {
                    'value': f'{match["b_start"]:,}–{match["b_end"]:,}',
                    'sort': match['b_start'],
                },
                'overlap': {
                    'value': f'{match["overlap"]:,} bp',
                    'sort': match['overlap'],
                },
            }
        )
    return rows


def search_cross_match_matrix(data: 'ReportData') -> Dict[str, Any]:
    """
    Summarise cross-matches as a group-by-group count matrix.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    dict
        With 'labels' (the groups involved) and 'counts', a mapping of
        ``'a\\tb'`` to the number of contested loci. Empty when nothing
        cross-matched.

    Notes
    -----
    The pair is stored unordered -- the two groups are symmetric partners in a
    conflict, not a direction -- so each pair is counted once under its
    alphabetically first member.
    """
    matches = data.stats.get('cross_matches', [])
    counts: Dict[str, int] = {}
    labels = set()
    for match in matches:
        a, b = sorted((match['a_cluster'], match['b_cluster']))
        labels.add(a)
        labels.add(b)
        key = f'{a}\t{b}'
        counts[key] = counts.get(key, 0) + 1
    return {'labels': sorted(labels), 'counts': counts}


def search_cluster_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the clustering rows: what each cluster was assembled from.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per cluster, with its components and their hit counts.
    """
    cluster_map = data.stats.get('cluster_map', {})
    hits_before = data.stats.get('hits_per_model_before_merge', {})
    by_group = {group.group_id: group for group in data.groups}

    rows = []
    for cluster, components in sorted(cluster_map.items()):
        group = by_group.get(cluster)
        before = sum(hits_before.get(component, 0) for component in components)
        after = group.n_unpaired if group else 0
        rows.append(
            {
                'cluster': {
                    'value': cluster,
                    'sort': cluster,
                    'colour': (
                        f'var(--group-{list(by_group).index(cluster)})'
                        if cluster in by_group
                        else None
                    ),
                },
                'components': ', '.join(sorted(components)),
                'n_components': {
                    'value': f'{len(components):,}',
                    'sort': len(components),
                },
                'before': {'value': f'{before:,}', 'sort': before},
                'after': {'value': f'{after:,}', 'sort': after},
                'merged': {
                    'value': f'{before - after:,}' if before >= after else '–',
                    'sort': before - after,
                },
            }
        )
    return rows


def search_contested_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the rows describing which model won each contested locus.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per (removed model, kept model) pair, from the nested-hit and
        cross-model filter steps.
    """
    rows = []
    for entry in data.stats.get('cross_model_removed', []):
        rows.append(
            {
                'step': 'Cross-model overlap',
                'removed': entry['removed'],
                'kept': entry['kept'],
                'hits': {'value': f'{entry["hits"]:,}', 'sort': entry['hits']},
            }
        )
    for removed, containers in sorted(data.stats.get('nested_removed', {}).items()):
        for container, count in sorted(containers.items()):
            rows.append(
                {
                    'step': 'Nested in partner',
                    'removed': removed,
                    'kept': container,
                    'hits': {'value': f'{count:,}', 'sort': count},
                }
            )
    for model, count in sorted(data.stats.get('excluded_not_in_map', {}).items()):
        rows.append(
            {
                'step': 'Not in pairing map',
                'removed': model,
                'kept': '—',
                'hits': {'value': f'{count:,}', 'sort': count},
            }
        )
    return rows


def search_model_overlap_table(data: 'ReportData') -> List[Dict[str, Any]]:
    """
    Build the pairwise model-overlap rows shown beside the heatmap.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        One row per pair of models sharing at least one locus, busiest first.
        Same-model pairs are labelled as redundancy rather than as a shared
        locus between two models.

    Notes
    -----
    This is the heatmap's table view: the figure shows the pattern, and the
    exact counts are readable here without hovering a cell.
    """
    cluster_map = data.stats.get('cluster_map', {})
    owner: Dict[str, str] = {}
    for cluster in sorted(cluster_map):
        for component in cluster_map[cluster]:
            owner.setdefault(component, cluster)

    rows = []
    for entry in data.stats.get('model_overlaps', []):
        a, b = entry['a'], entry['b']
        cluster_a, cluster_b = owner.get(a), owner.get(b)
        if a == b:
            relation = 'same model (redundant hits)'
        elif cluster_a is not None and cluster_a == cluster_b:
            relation = f'same cluster ({cluster_a})'
        else:
            relation = 'different clusters' if cluster_a or cluster_b else 'unclustered'
        rows.append(
            {
                'a': a,
                'b': b,
                'relation': relation,
                'hits': {'value': f'{entry["hits"]:,}', 'sort': entry['hits']},
            }
        )
    return rows
