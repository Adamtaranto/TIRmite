"""
HTML reporting for TIRmite.

This package turns the internal results of a TIRmite run into a single
self-contained HTML file. It is split so that the *data* side and the
*rendering* side never depend on each other:

``stats``
    Statistics shared by the plain-text summaries and the HTML report.

Modules that import third-party rendering libraries (``jinja2``,
``matplotlib``) are deliberately **not** imported here. They are pulled in
lazily inside the functions that need them so that importing ``tirmite`` --
or starting the CLI -- never pays for them.
"""

from tirmite.report.stats import (
    PairSummary,
    format_pair_summary,
    pair_summary_stats,
)

__all__ = [
    'PairSummary',
    'format_pair_summary',
    'pair_summary_stats',
]
