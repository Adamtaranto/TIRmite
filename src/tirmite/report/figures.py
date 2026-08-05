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
from typing import Any, Dict, List, Optional, Sequence

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
