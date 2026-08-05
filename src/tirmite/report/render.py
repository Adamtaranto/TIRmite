"""
Render report data to a single self-contained HTML file.

The output has no external references at all -- no stylesheet links, no script
sources, no web fonts, no remote images. Reports are routinely opened from disk,
emailed and archived alongside the run that produced them, so anything fetched
over the network would be missing exactly when the report is needed. CSS and
JavaScript are therefore inlined at render time from the package's ``assets``
directory.

``jinja2`` is imported inside :func:`render_report` rather than at module scope,
so importing :mod:`tirmite.report` -- which the CLI does on every run -- does not
pay for it.
"""

from importlib import resources
import logging
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Sequence, Union

from tirmite.report.model import FigureSpec, ReportData
from tirmite.report.stats import (
    contig_table,
    filter_table,
    group_table,
    model_table,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['inline_asset', 'render_report', 'write_pair_report']

_PACKAGE = 'tirmite.report'
_ASSET_DIR = 'assets'
_TEMPLATE_DIR = 'templates'

# Loaded in order; later files may depend on earlier ones.
_SCRIPTS = ('report-core.js', 'track-browser.js', 'stats-table.js', 'msa-panel.js')
_STYLES = ('report.css',)

# Matches an end tag that would terminate the enclosing <script> or <style>
# element, in any casing and with the whitespace HTML tolerates.
_CLOSING_TAG = re.compile(r'</\s*(script|style)', re.IGNORECASE)


def inline_asset(name: str, subdir: str = _ASSET_DIR) -> str:
    """
    Read a packaged asset for inlining into the document.

    Parameters
    ----------
    name : str
        File name within the directory.
    subdir : str, optional
        Directory inside the ``tirmite.report`` package to read from.

    Returns
    -------
    str
        The file's text, with any sequence that would close the enclosing
        script or style element neutralised.

    Raises
    ------
    FileNotFoundError
        If the asset is not present, which means the wheel was built without
        it.

    Notes
    -----
    Read through :mod:`importlib.resources` rather than by path arithmetic on
    ``__file__``, so the report works from a zipped distribution as well as
    from a source checkout.
    """
    try:
        text = (
            resources.files(_PACKAGE)
            .joinpath(subdir)
            .joinpath(name)
            .read_text(encoding='utf-8')
        )
    except (FileNotFoundError, ModuleNotFoundError, OSError) as exc:
        raise FileNotFoundError(
            f'Report asset {subdir}/{name} is missing from {_PACKAGE}. The '
            'package was probably built without its templates and assets; '
            "check the 'artifacts' list in pyproject.toml."
        ) from exc
    # An asset containing the literal '</script' would end the element early.
    # No legitimate CSS or JS needs it unescaped.
    return _CLOSING_TAG.sub(r'<\\/\1', text)


def _environment() -> Any:
    """
    Build the Jinja2 environment used for report templates.

    Returns
    -------
    jinja2.Environment
        Configured with HTML autoescaping and whitespace trimming.
    """
    from jinja2 import Environment, PackageLoader, select_autoescape

    return Environment(
        loader=PackageLoader(_PACKAGE, _TEMPLATE_DIR),
        autoescape=select_autoescape(['html', 'xml', 'html.j2']),
        trim_blocks=True,
        lstrip_blocks=True,
    )


def render_report(
    data: ReportData,
    *,
    figures: Optional[Sequence[FigureSpec]] = None,
    template: str = 'pair_report.html.j2',
) -> str:
    """
    Render report data to a complete HTML document.

    Parameters
    ----------
    data : ReportData
        The report to render.
    figures : sequence of FigureSpec, optional
        Static figures to embed. Overrides any already on `data`.
    template : str, default 'pair_report.html.j2'
        Template to use.

    Returns
    -------
    str
        A complete, self-contained HTML document.
    """
    if figures is not None:
        data.figures = list(figures)

    env = _environment()
    context: Dict[str, Any] = {
        'data': data,
        'payload': data.to_json(),
        'styles': [inline_asset(name) for name in _STYLES],
        'scripts': [inline_asset(name) for name in _SCRIPTS],
        'summary': _headline_numbers(data),
        'stats_groups': group_table(data),
        'stats_models': model_table(data),
        'stats_contigs': contig_table(data),
        'stats_filters': filter_table(data.filter_stats),
    }
    # jinja2 is untyped under `no_site_packages`, so render() is Any.
    html: str = env.get_template(template).render(**context)
    return html


def _headline_numbers(data: ReportData) -> List[Dict[str, Any]]:
    """
    Build the stat tiles shown above the report.

    Parameters
    ----------
    data : ReportData
        The report.

    Returns
    -------
    list of dict
        Tiles with 'label', 'value' and 'detail' keys.

    Notes
    -----
    These are the four numbers a reader wants before scrolling anywhere: how
    many elements were predicted, how much of the hit set paired up, how many
    sequences are involved and how many pairing procedures ran.
    """
    n_elements = len(data.elements)
    n_hits = data.hits.n
    n_paired_hits = sum(1 for ix in data.hits.pair_ix if ix is not None)
    pairing_rate = (n_paired_hits / n_hits) if n_hits else 0.0

    lengths = sorted(e.length for e in data.elements)
    if lengths:
        median_length = lengths[len(lengths) // 2]
        length_detail = f'median {median_length:,} bp'
    else:
        length_detail = 'no elements predicted'

    return [
        {
            'label': 'Predicted elements',
            'value': f'{n_elements:,}',
            'detail': length_detail,
        },
        {
            'label': 'Termini paired',
            'value': f'{pairing_rate:.0%}',
            'detail': f'{n_paired_hits:,} of {n_hits:,} hits',
        },
        {
            'label': 'Sequences with hits',
            'value': f'{len(data.contigs):,}',
            'detail': f'{sum(c.n_pairs for c in data.contigs):,} carry a pair',
        },
        {
            'label': 'Pairing groups',
            'value': f'{len(data.groups):,}',
            'detail': f'{len(data.models):,} terminus model(s)',
        },
    ]


def write_pair_report(
    data: ReportData,
    outpath: Union[str, Path],
    *,
    figures: Optional[Sequence[FigureSpec]] = None,
) -> Path:
    """
    Render a pairing report and write it to disk.

    Parameters
    ----------
    data : ReportData
        The report to render.
    outpath : str or pathlib.Path
        Destination file.
    figures : sequence of FigureSpec, optional
        Static figures to embed.

    Returns
    -------
    pathlib.Path
        The path written.
    """
    path = Path(outpath)
    path.parent.mkdir(parents=True, exist_ok=True)
    html = render_report(data, figures=figures)
    path.write_text(html, encoding='utf-8')
    logger.info(f'Wrote HTML report to {path} ({len(html) / 1024:.0f} KiB)')
    return path
