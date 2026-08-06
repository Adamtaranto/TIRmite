"""
Static statistical figures for the report, rendered as inline SVG.

Matplotlib is imported inside the functions that need it rather than at module
scope. It is a heavy import and most TIRmite runs never build a report, so the
CLI should not pay for it.

Figures are styled for print: single-column width, 7 pt sans, hairline axes,
no grid. The house style is applied through ``plt.rc_context`` rather than by
mutating ``plt.rcParams``, because TIRmite is importable as a library and must
not reconfigure its host's plotting.
"""

import io
import logging
import re
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from tirmite.report.model import FigureSpec, ReportData

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['build_figures', 'nature_rcparams', 'namespace_svg_ids']

# Nature's single-column measure is 89 mm; the double column is 183 mm.
SINGLE_COLUMN_IN = 89 / 25.4
FIGURE_HEIGHT_IN = SINGLE_COLUMN_IN * 0.68


def nature_rcparams() -> Dict[str, Any]:
    """
    Return the report's figure style as matplotlib rc parameters.

    Returns
    -------
    dict
        Parameters suitable for ``matplotlib.pyplot.rc_context``.

    Notes
    -----
    Several sans faces are listed rather than one: continuous integration
    machines rarely have Helvetica, and a missing font makes matplotlib fall
    back with a warning on every figure.
    """
    return {
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'axes.titlesize': 7.5,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6.5,
        'axes.linewidth': 0.5,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'lines.linewidth': 1.0,
        'patch.linewidth': 0.5,
        'legend.frameon': False,
        'figure.dpi': 150,
        'savefig.transparent': True,
        'svg.fonttype': 'none',
    }


# matplotlib emits internal ids such as 'p1a2b3c4d5' for clip paths and glyph
# definitions. Several inline SVGs in one document collide on those ids and
# render blank or corrupted, so each figure's ids get a unique prefix.
_ID_ATTR = re.compile(r'\bid="([^"]+)"')
_HREF = re.compile(r'\b((?:xlink:)?href)="#([^"]+)"')
_URL_REF = re.compile(r'url\(#([^)]+)\)')
_XML_PROLOGUE = re.compile(r'^\s*(<\?xml[^>]*\?>|<!DOCTYPE[^>]*>)\s*', re.IGNORECASE)
# matplotlib writes an RDF block naming itself and its homepage. It is dead
# weight in an embedded figure, and the URL it carries would otherwise be the
# only external-looking reference in a document that has none.
_METADATA = re.compile(r'<metadata>.*?</metadata>\s*', re.DOTALL | re.IGNORECASE)


def namespace_svg_ids(svg: str, prefix: str) -> str:
    """
    Rewrite every internal id in an SVG so it cannot collide with another.

    Parameters
    ----------
    svg : str
        SVG markup.
    prefix : str
        Prefix to apply to each id.

    Returns
    -------
    str
        The markup with ids and their references rewritten.
    """
    ids = set(_ID_ATTR.findall(svg))
    if not ids:
        return svg

    def rename(name: str) -> str:
        """
        Prefix an id, leaving references to anything outside this SVG alone.

        Parameters
        ----------
        name : str
            The id or reference target.

        Returns
        -------
        str
            The namespaced id, or `name` unchanged if it is not defined here.
        """
        return f'{prefix}-{name}' if name in ids else name

    svg = _ID_ATTR.sub(lambda m: f'id="{rename(m.group(1))}"', svg)
    svg = _HREF.sub(lambda m: f'{m.group(1)}="#{rename(m.group(2))}"', svg)
    svg = _URL_REF.sub(lambda m: f'url(#{rename(m.group(1))})', svg)
    return svg


def _to_svg(fig: Any, figure_id: str) -> str:
    """
    Serialise a figure to inline-ready SVG markup.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to serialise. It is closed before returning.
    figure_id : str
        Slug used to namespace the figure's internal ids.

    Returns
    -------
    str
        SVG markup with the XML prologue removed and ids namespaced.
    """
    import matplotlib.pyplot as plt

    buffer = io.StringIO()
    fig.savefig(buffer, format='svg', bbox_inches='tight')
    plt.close(fig)
    svg = _XML_PROLOGUE.sub('', buffer.getvalue()).lstrip()
    # A second prologue can follow the first (declaration then doctype).
    svg = _XML_PROLOGUE.sub('', svg).lstrip()
    svg = _METADATA.sub('', svg)
    return namespace_svg_ids(svg, figure_id)


def _model_series(
    values_by_model: Dict[str, List[float]],
) -> List[tuple]:
    """
    Order per-model series and assign each a distinct colour.

    Parameters
    ----------
    values_by_model : dict
        Model name to its values.

    Returns
    -------
    list of (str, list, str)
        ``(label, values, colour)`` in plotting order.

    Notes
    -----
    Colours are taken from the report palette in its fixed order and are never
    cycled: a repeated hue in one chart would say two models are the same
    series. Past the palette's width the remaining models are folded into a
    single neutral "other models" series, which is honest about the loss in a
    way that reusing hues is not.
    """
    from tirmite.report.palette import GROUP_COLOURS, NEUTRAL_GREY

    # Busiest models first, so the fold drops the least informative ones.
    ordered = sorted(values_by_model.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    limit = len(GROUP_COLOURS)

    series = [
        (name, values, GROUP_COLOURS[i])
        for i, (name, values) in enumerate(ordered[:limit])
    ]
    if len(ordered) > limit:
        rest = [v for _, values in ordered[limit:] for v in values]
        series.append((f'other models (n={len(ordered) - limit})', rest, NEUTRAL_GREY))
    return series


def _legend_outside(ax: Any, n_entries: int) -> None:
    """
    Place a legend above the axes rather than inside them.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add the legend to.
    n_entries : int
        Number of legend entries, used to choose the column count.

    Returns
    -------
    None
        Modifies the axes in place.

    Notes
    -----
    A legend inside the axes sits on top of the data, and in a distribution
    plot it lands exactly where the tallest bars are. Anchoring it above the
    axes costs a little height and obscures nothing. ``bbox_inches='tight'``
    at save time keeps it inside the figure bounds.
    """
    ax.legend(
        loc='lower left',
        bbox_to_anchor=(0, 1.01, 1, 0.15),
        mode='expand',
        borderaxespad=0,
        ncol=min(3, max(1, n_entries)),
        handlelength=1.4,
        columnspacing=1.2,
    )


def _finish(ax: Any) -> None:
    """
    Apply the shared axis treatment.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to adjust.

    Returns
    -------
    None
        Modifies the axes in place.
    """
    ax.tick_params(length=2, width=0.5, pad=2)
    for spine in ('left', 'bottom'):
        ax.spines[spine].set_linewidth(0.5)


def _element_lengths(data: ReportData, plt: Any) -> Optional[FigureSpec]:
    """
    Plot the distribution of predicted element lengths.

    Parameters
    ----------
    data : ReportData
        The report.
    plt : module
        ``matplotlib.pyplot``.

    Returns
    -------
    FigureSpec or None
        The figure, or None if no elements were predicted.
    """
    if not data.elements:
        return None

    by_group: Dict[int, List[int]] = {}
    for element in data.elements:
        by_group.setdefault(element.group_i, []).append(element.length)

    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_IN, FIGURE_HEIGHT_IN))
    all_lengths = [e.length for e in data.elements]
    bins = min(40, max(8, len(all_lengths) // 4))

    for group_i, lengths in sorted(by_group.items()):
        group = data.groups[group_i]
        ax.hist(
            lengths,
            bins=bins,
            range=(min(all_lengths), max(all_lengths) + 1),
            histtype='step',
            linewidth=1.0,
            color=group.colour,
            label=group.label,
        )

    ax.set_xlabel('Predicted element length (bp)')
    ax.set_ylabel('Elements')
    # One series names itself in the title; two or more need a legend.
    if len(by_group) > 1:
        _legend_outside(ax, len(by_group))
    _finish(ax)

    median_length = sorted(all_lengths)[len(all_lengths) // 2]
    return FigureSpec(
        id='element-lengths',
        title='Predicted element lengths',
        caption=(
            f'{len(all_lengths):,} elements, median {median_length:,} bp, '
            f'range {min(all_lengths):,}–{max(all_lengths):,} bp.'
        ),
        svg=_to_svg(fig, 'element-lengths'),
    )


def _pairing_outcome(data: ReportData, plt: Any) -> Optional[FigureSpec]:
    """
    Plot paired against unpaired hits for each pairing group.

    Parameters
    ----------
    data : ReportData
        The report.
    plt : module
        ``matplotlib.pyplot``.

    Returns
    -------
    FigureSpec or None
        The figure, or None if there are no pairing groups.
    """
    if not data.groups:
        return None

    labels = [group.label for group in data.groups]
    paired = [group.n_pairs * 2 for group in data.groups]
    unpaired = [group.n_unpaired for group in data.groups]
    positions = range(len(labels))

    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_IN, FIGURE_HEIGHT_IN))
    ax.barh(
        list(positions),
        paired,
        height=0.55,
        color=[group.colour for group in data.groups],
        label='paired',
    )
    ax.barh(
        list(positions),
        unpaired,
        height=0.55,
        left=paired,
        color='#c9c9c3',
        label='unpaired',
    )
    ax.set_yticks(list(positions))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel('Terminus hits')
    _legend_outside(ax, 2)
    _finish(ax)

    total_paired = sum(paired)
    total = total_paired + sum(unpaired)
    rate = total_paired / total if total else 0
    return FigureSpec(
        id='pairing-outcome',
        title='Pairing outcome by group',
        caption=(
            f'{total_paired:,} of {total:,} terminus hits ({rate:.0%}) were '
            'assigned a partner.'
        ),
        svg=_to_svg(fig, 'pairing-outcome'),
    )


def _model_coverage(data: ReportData, plt: Any) -> Optional[FigureSpec]:
    """
    Plot how completely hits cover their model, per model.

    Parameters
    ----------
    data : ReportData
        The report.
    plt : module
        ``matplotlib.pyplot``.

    Returns
    -------
    FigureSpec or None
        The figure, or None if no hit carries model coordinates.
    """
    coverages: Dict[str, List[float]] = {}
    for i in range(data.hits.n):
        value = data.hits.model_cov[i]
        if value is None:
            continue
        coverages.setdefault(data.models[data.hits.model_i[i]].name, []).append(value)

    if not coverages:
        return None

    series = _model_series(coverages)
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_IN, FIGURE_HEIGHT_IN))
    for name, values, colour in series:
        ax.hist(
            values,
            bins=20,
            range=(0, 1),
            histtype='step',
            linewidth=1.0,
            color=colour,
            label=name,
        )

    mincov = data.filter_stats.get('mincov')
    if isinstance(mincov, (int, float)):
        # Drawing the filter makes the left edge of the distribution legible:
        # everything below it was removed before pairing.
        ax.axvline(float(mincov), color='#75756e', linewidth=0.75)
        ax.annotate(
            f'--mincov {mincov}',
            xy=(float(mincov), 1),
            xycoords=('data', 'axes fraction'),
            xytext=(2, -6),
            textcoords='offset points',
            fontsize=6,
            color='#75756e',
        )

    ax.set_xlabel('Fraction of model matched')
    ax.set_ylabel('Hits')
    ax.set_xlim(0, 1)
    if len(series) > 1:
        _legend_outside(ax, len(series))
    _finish(ax)

    return FigureSpec(
        id='model-coverage',
        title='Model coverage',
        caption=(
            'How much of each terminus model its hits matched, from the '
            'alignment coordinates. Not the same quantity as the span-based '
            'coverage used by --mincov, which can exceed 1.'
        ),
        svg=_to_svg(fig, 'model-coverage'),
    )


def _hits_per_contig(data: ReportData, plt: Any) -> Optional[FigureSpec]:
    """
    Plot the sequences carrying the most hits.

    Parameters
    ----------
    data : ReportData
        The report.
    plt : module
        ``matplotlib.pyplot``.

    Returns
    -------
    FigureSpec or None
        The figure, or None if fewer than two sequences carry hits.
    """
    if len(data.contigs) < 2:
        return None

    top = sorted(data.contigs, key=lambda c: -c.n_hits)[:20]
    positions = range(len(top))

    fig, ax = plt.subplots(
        figsize=(SINGLE_COLUMN_IN, max(FIGURE_HEIGHT_IN, 0.14 * len(top) + 0.5))
    )
    ax.barh(
        list(positions),
        [c.n_hits for c in top],
        height=0.6,
        color='#2a78d6',
    )
    ax.set_yticks(list(positions))
    ax.set_yticklabels([c.name for c in top])
    ax.invert_yaxis()
    ax.set_xlabel('Terminus hits')
    _finish(ax)

    shown = (
        f'top {len(top)} of {len(data.contigs):,}'
        if len(data.contigs) > len(top)
        else 'all'
    )
    return FigureSpec(
        id='hits-per-contig',
        title='Hits per sequence',
        caption=f'Sequences carrying the most terminus hits ({shown}).',
        svg=_to_svg(fig, 'hits-per-contig'),
    )


def _evalues(data: ReportData, plt: Any) -> Optional[FigureSpec]:
    """
    Plot the e-value distribution per model.

    Parameters
    ----------
    data : ReportData
        The report.
    plt : module
        ``matplotlib.pyplot``.

    Returns
    -------
    FigureSpec or None
        The figure, or None if no hit has a usable e-value.
    """
    import math

    by_model: Dict[str, List[float]] = {}
    for i in range(data.hits.n):
        value = data.hits.evalue[i]
        if value is None or value <= 0:
            continue
        by_model.setdefault(data.models[data.hits.model_i[i]].name, []).append(
            math.log10(value)
        )

    if not by_model:
        return None

    series = _model_series(by_model)
    fig, ax = plt.subplots(figsize=(SINGLE_COLUMN_IN, FIGURE_HEIGHT_IN))
    for name, values, colour in series:
        ax.hist(
            values,
            bins=24,
            histtype='step',
            linewidth=1.0,
            color=colour,
            label=name,
        )

    ax.set_xlabel('log$_{10}$(e-value)')
    ax.set_ylabel('Hits')
    if len(series) > 1:
        _legend_outside(ax, len(series))
    _finish(ax)

    return FigureSpec(
        id='evalues',
        title='Hit significance',
        caption=(
            'E-value distribution per terminus model. Hits with a '
            'non-positive or missing e-value are omitted.'
        ),
        svg=_to_svg(fig, 'evalues'),
    )


def _model_order(data: ReportData) -> Any:
    """
    Order query models by cluster, then by name, and locate the cluster blocks.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    tuple of (list, list)
        ``(models, blocks)`` where models is the axis order and blocks is a
        list of ``(label, first_index, last_index)`` for each cluster
        occupying a contiguous run.

    Notes
    -----
    A model listed in two clusters is placed under the first alphabetically.
    Any other choice is equally arbitrary; what matters is that it is stated,
    which the report does in a warning.
    """
    overlaps = data.stats.get('model_overlaps', [])
    cluster_map = data.stats.get('cluster_map', {})

    models = set()
    for entry in overlaps:
        models.add(entry['a'])
        models.add(entry['b'])
    models |= set(data.stats.get('hits_per_model_before_merge', {}))
    if not models:
        return [], []

    owner: Dict[str, str] = {}
    for cluster in sorted(cluster_map):
        for component in cluster_map[cluster]:
            owner.setdefault(component, cluster)

    # Case-insensitive, so 'hAT_cluster' does not sort after 'Mut_cluster'
    # merely because of its leading lowercase letter. Unclustered models sort
    # after every cluster, under a high sentinel key.
    ordered = sorted(
        models,
        key=lambda m: (owner.get(m, '￿').casefold(), m.casefold(), m),
    )
    return ordered, _rebuild_blocks(ordered, data)


def _model_overlap_heatmap(data: ReportData, plt: Any) -> Optional[FigureSpec]:
    """
    Plot which query models share hit loci, as a model-by-model heatmap.

    Parameters
    ----------
    data : ReportData
        The report.
    plt : module
        ``matplotlib.pyplot``.

    Returns
    -------
    FigureSpec or None
        The figure, or None if no two models share a locus.

    Notes
    -----
    A sequential single-hue ramp, because the value is a magnitude: a
    categorical or diverging scheme would imply a distinction the counts do not
    carry. The diagonal counts two hits of the *same* model at one locus, which
    is redundancy rather than confusion between models, so it is drawn but
    described separately in the caption.

    Boxes mark models that share a cluster. Overlaps inside a box are expected
    -- that is what a cluster asserts. Colour outside the boxes is the finding.
    """
    overlaps = data.stats.get('model_overlaps', [])
    if not overlaps:
        return None

    models, blocks = _model_order(data)
    if len(models) < 2:
        return None

    # Past this the cells are too small to read and the labels collide.
    limit = 40
    truncated = False
    if len(models) > limit:
        busiest: Dict[str, int] = {}
        for entry in overlaps:
            busiest[entry['a']] = busiest.get(entry['a'], 0) + entry['hits']
            busiest[entry['b']] = busiest.get(entry['b'], 0) + entry['hits']
        keep = set(sorted(busiest, key=lambda m: -busiest[m])[:limit])
        models = [m for m in models if m in keep]
        blocks = _rebuild_blocks(models, data)
        truncated = True

    index = {model: i for i, model in enumerate(models)}
    size = len(models)
    matrix = [[0 for _ in range(size)] for _ in range(size)]
    for entry in overlaps:
        i = index.get(entry['a'])
        j = index.get(entry['b'])
        if i is None or j is None:
            continue
        matrix[i][j] = entry['hits']
        matrix[j][i] = entry['hits']

    peak = max((max(row) for row in matrix), default=0)
    if not peak:
        return None

    side = min(6.5, max(SINGLE_COLUMN_IN, 0.22 * size + 1.4))
    fig, ax = plt.subplots(figsize=(side, side))

    image = ax.imshow(
        matrix,
        cmap='Blues',
        vmin=0,
        vmax=peak,
        interpolation='nearest',
    )

    ax.set_xticks(range(size))
    ax.set_yticks(range(size))
    ax.set_xticklabels(models, rotation=90, fontsize=5.5)
    ax.set_yticklabels(models, fontsize=5.5)
    ax.tick_params(length=0, pad=1)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Cell counts, while they still fit.
    if size <= 14:
        for i in range(size):
            for j in range(size):
                value = matrix[i][j]
                if not value:
                    continue
                ax.text(
                    j,
                    i,
                    str(value),
                    ha='center',
                    va='center',
                    fontsize=5.5,
                    color='white' if value > peak * 0.6 else '#1a1a19',
                )

    # A cluster's members are expected to overlap; the box says so, leaving
    # colour outside it as the thing worth looking at.
    from matplotlib.patches import Rectangle

    for _label, lo, hi in blocks:
        ax.add_patch(
            Rectangle(
                (lo - 0.5, lo - 0.5),
                hi - lo + 1,
                hi - lo + 1,
                fill=False,
                edgecolor='#eb6834',
                linewidth=1.0,
            )
        )

    bar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    bar.ax.tick_params(labelsize=6, length=2, width=0.5)
    bar.outline.set_linewidth(0.5)
    bar.set_label('Overlapping hits', fontsize=6.5)

    caption = (
        'Pairs of query models whose hits share a locus, counted before '
        'clustering. Boxes enclose models that belong to the same cluster, '
        'where overlap is expected; colour outside a box is two unrelated '
        'models claiming the same sequence. The diagonal counts overlapping '
        'hits from one model, which is redundancy rather than confusion.'
    )
    if truncated:
        caption += f' Showing the {limit} models with the most overlap.'

    return FigureSpec(
        id='model-overlaps',
        title='Shared hit loci between models',
        caption=caption,
        svg=_to_svg(fig, 'model-overlaps'),
        wide=True,
    )


def _overlap_matrix(
    models: Sequence[str], overlaps: Sequence[Dict[str, Any]]
) -> List[List[int]]:
    """
    Build a symmetric count matrix over the given models.

    Parameters
    ----------
    models : sequence of str
        Axis order.
    overlaps : sequence of dict
        Records with 'a', 'b' and 'hits'.

    Returns
    -------
    list of list of int
        Counts, symmetric, zero where two models never share a locus.
    """
    index = {model: i for i, model in enumerate(models)}
    size = len(models)
    matrix = [[0 for _ in range(size)] for _ in range(size)]
    for entry in overlaps:
        i = index.get(entry['a'])
        j = index.get(entry['b'])
        if i is None or j is None:
            continue
        matrix[i][j] = entry['hits']
        matrix[j][i] = entry['hits']
    return matrix


def _model_overlap_clustered_heatmap(
    data: ReportData, plt: Any
) -> Optional[FigureSpec]:
    """
    Plot shared hit loci with the axes ordered by hierarchical clustering.

    Parameters
    ----------
    data : ReportData
        The report.
    plt : module
        ``matplotlib.pyplot``.

    Returns
    -------
    FigureSpec or None
        The figure, or None when there is nothing to cluster.

    Notes
    -----
    The companion figure orders models by their declared cluster, which asks
    "do the hits agree with the cluster map?". This one lets the hits speak
    first: models are ordered by average-linkage clustering on how much they
    share, so groups emerge from the data whether or not a cluster map was
    supplied. Reading them together is the point -- a block here that crosses a
    box there is a cluster map at odds with its own evidence.

    Distances are Dice dissimilarities on the overlap counts, so a model that
    shares most of its few hits ranks as close as one sharing many of many.
    """
    overlaps = data.stats.get('model_overlaps', [])
    if not overlaps:
        return None

    models, _blocks = _model_order(data)
    if len(models) < 3:
        # With two models there is one join and nothing a tree can reveal.
        return None

    limit = 40
    truncated = False
    if len(models) > limit:
        busiest: Dict[str, int] = {}
        for entry in overlaps:
            busiest[entry['a']] = busiest.get(entry['a'], 0) + entry['hits']
            busiest[entry['b']] = busiest.get(entry['b'], 0) + entry['hits']
        keep = set(sorted(busiest, key=lambda m: -busiest[m])[:limit])
        models = [m for m in models if m in keep]
        truncated = True

    from tirmite.report.cluster import overlap_distances, upgma

    counts = {
        (entry['a'], entry['b']): entry['hits']
        for entry in overlaps
        if entry['a'] != entry['b']
    }
    hits_per_model = data.stats.get('hits_per_model_before_merge', {})
    distances = overlap_distances(models, counts, hits_per_model)
    tree = upgma(models, distances)

    ordered = tree.order
    matrix = _overlap_matrix(ordered, overlaps)
    peak = max((max(row) for row in matrix), default=0)
    if not peak:
        return None

    size = len(ordered)
    side = min(6.5, max(SINGLE_COLUMN_IN, 0.22 * size + 1.4))
    fig = plt.figure(figsize=(side, side * 1.2))
    grid = fig.add_gridspec(
        2, 2, height_ratios=[1, 4], width_ratios=[24, 1], hspace=0.04, wspace=0.06
    )
    tree_ax = fig.add_subplot(grid[0, 0])
    ax = fig.add_subplot(grid[1, 0], sharex=tree_ax)
    bar_ax = fig.add_subplot(grid[1, 1])

    image = ax.imshow(matrix, cmap='Blues', vmin=0, vmax=peak, interpolation='nearest')

    ax.set_xticks(range(size))
    ax.set_yticks(range(size))
    ax.set_xticklabels(ordered, rotation=90, fontsize=5.5)
    ax.set_yticklabels(ordered, fontsize=5.5)
    ax.tick_params(length=0, pad=1)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if size <= 14:
        for i in range(size):
            for j in range(size):
                value = matrix[i][j]
                if not value:
                    continue
                ax.text(
                    j,
                    i,
                    str(value),
                    ha='center',
                    va='center',
                    fontsize=5.5,
                    color='white' if value > peak * 0.6 else '#1a1a19',
                )

    # The tree, drawn as the usual brackets above its columns.
    for merge in tree.merges:
        tree_ax.plot(
            [merge.left_x, merge.left_x, merge.right_x, merge.right_x],
            [merge.left_height, merge.height, merge.height, merge.right_height],
            color='#4a4a46',
            linewidth=0.7,
            solid_joinstyle='miter',
        )
    tree_ax.set_ylim(0, max(tree.height * 1.08, 1e-6))
    tree_ax.set_xlim(-0.5, size - 0.5)
    tree_ax.set_ylabel('Distance', fontsize=6)
    tree_ax.tick_params(axis='y', labelsize=5.5, length=2, width=0.5, pad=1)
    tree_ax.tick_params(axis='x', length=0, labelbottom=False)
    for name, spine in tree_ax.spines.items():
        spine.set_visible(name == 'left')
        if name == 'left':
            spine.set_linewidth(0.5)

    bar = fig.colorbar(image, cax=bar_ax)
    bar.ax.tick_params(labelsize=6, length=2, width=0.5)
    bar.outline.set_linewidth(0.5)
    bar.set_label('Overlapping hits', fontsize=6.5)

    caption = (
        'The same counts, with models ordered by average-linkage clustering on '
        'how much they share rather than by their declared cluster. Groups here '
        'come from the hits alone, so comparing this with the boxed figure '
        'shows whether the cluster map agrees with its own evidence. Distance '
        'is a Dice dissimilarity, so a model sharing most of its few hits ranks '
        'as close as one sharing many of many.'
    )
    if truncated:
        caption += f' Showing the {limit} models with the most overlap.'

    return FigureSpec(
        id='model-overlaps-clustered',
        title='Shared hit loci, clustered',
        caption=caption,
        svg=_to_svg(fig, 'model-overlaps-clustered'),
        wide=True,
    )


def _rebuild_blocks(models: List[str], data: ReportData) -> List[Any]:
    """
    Recompute cluster blocks for a reduced model list.

    Parameters
    ----------
    models : list of str
        Models in axis order, after any truncation.
    data : ReportData
        The report, for its cluster map.

    Returns
    -------
    list of tuple
        ``(label, first_index, last_index)`` per contiguous cluster run.
    """
    cluster_map = data.stats.get('cluster_map', {})
    owner: Dict[str, str] = {}
    for cluster in sorted(cluster_map):
        for component in cluster_map[cluster]:
            owner.setdefault(component, cluster)

    blocks = []
    start = 0
    for i, model in enumerate(models):
        current = owner.get(model)
        nxt = owner.get(models[i + 1]) if i + 1 < len(models) else None
        if current != nxt:
            if current is not None:
                blocks.append((current, start, i))
            start = i + 1
    return blocks


def build_figures(data: ReportData) -> List[FigureSpec]:
    """
    Build every static figure the report can show for this run.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of FigureSpec
        Figures that had data to show, in display order. Empty when
        matplotlib is unusable.

    Notes
    -----
    A failure here costs the figures and nothing else. The rest of the report
    is already computed, and the run's real outputs were written long before,
    so an unusable plotting backend must not turn into a failed run.
    """
    try:
        import matplotlib

        # Agg before pyplot: report generation is headless, and picking an
        # interactive backend on a machine without a display raises.
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as exc:  # noqa: BLE001 - figures are optional
        logger.warning(f'matplotlib is unavailable, skipping figures: {exc}')
        return []

    # A search report has no pairs, so the pairing-outcome bar would show every
    # hit as "unpaired" -- true but vacuous. The overlap heatmap answers the
    # question a search run actually raises, and needs pre-merge data a pairing
    # run never collects.
    builders: Tuple[Callable[[ReportData, Any], Optional[FigureSpec]], ...]
    if data.kind == 'search':
        builders = (
            _model_overlap_heatmap,
            _model_overlap_clustered_heatmap,
            _model_coverage,
            _hits_per_contig,
            _evalues,
        )
    else:
        builders = (
            _element_lengths,
            _pairing_outcome,
            _model_coverage,
            _hits_per_contig,
            _evalues,
        )

    figures: List[FigureSpec] = []
    with plt.rc_context(nature_rcparams()):
        for builder in builders:
            try:
                figure = builder(data, plt)
            except Exception as exc:  # noqa: BLE001 - one bad figure is not fatal
                logger.warning(f'Could not build a report figure: {exc}')
                continue
            if figure is not None:
                figures.append(figure)

    return figures


def figure_ids(figures: Sequence[FigureSpec]) -> List[str]:
    """
    Return the ids of a figure sequence.

    Parameters
    ----------
    figures : sequence of FigureSpec
        Figures to inspect.

    Returns
    -------
    list of str
        Their ids, in order.
    """
    return [figure.id for figure in figures]
