"""
HTML reporting for TIRmite.

This package turns the internal results of a TIRmite run into a single
self-contained HTML file. It is split so that the *data* side and the
*rendering* side never depend on each other:

``model``
    Serialisable dataclasses describing a report.
``collect``
    Turns the internals of a run into those dataclasses.
``layout``
    Stacks overlapping annotations into track rows.
``palette``
    Colours shared by every part of the report.
``stats``
    Statistics shared by the plain-text summaries and the HTML report.

Modules that import third-party rendering libraries (``jinja2``,
``matplotlib``) are deliberately **not** imported here. They are pulled in
lazily inside the functions that need them so that importing ``tirmite`` --
or starting the CLI -- never pays for them.
"""

from tirmite.report.collect import PairReportAccumulator, model_truncation
from tirmite.report.layout import assign_rows, stack_contig
from tirmite.report.model import ReportData
from tirmite.report.palette import group_colours
from tirmite.report.stats import (
    PairSummary,
    format_pair_summary,
    pair_summary_stats,
)

__all__ = [
    'PairReportAccumulator',
    'PairSummary',
    'ReportData',
    'assign_rows',
    'format_pair_summary',
    'group_colours',
    'model_truncation',
    'pair_summary_stats',
    'stack_contig',
]
