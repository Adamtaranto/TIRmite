"""Tests for tirmite.report.figures — static inline-SVG figures."""

import re

import pytest
from test_report_model import make_element, make_hit, make_index

from tirmite.core.pairing import PairingConfig
from tirmite.report.collect import PairReportAccumulator
from tirmite.report.figures import (
    build_figures,
    namespace_svg_ids,
    nature_rcparams,
)

plt = pytest.importorskip('matplotlib.pyplot')

ID_ATTR = re.compile(r'\bid="([^"]+)"')


@pytest.fixture
def report(hit_table_factory):
    """A report with two groups, several contigs and a spread of hits."""
    hits = []
    rows = []
    pairs = []
    uid = 0
    for contig, n in (('chr1', 6), ('chr2', 4)):
        for i in range(n):
            start = 1000 + i * 20_000
            end = start + 119
            hits.append(make_hit(uid, 'TIR_L', contig, start, end))
            rows.append(
                {
                    'model': 'TIR_L',
                    'target': contig,
                    'hitStart': str(start),
                    'hitEnd': str(end),
                    'hmmStart': str(1 + i),
                    'hmmEnd': '120',
                    'evalue': f'1e-{10 + i}',
                }
            )
            uid += 1
            rstart = start + 3000 + i * 400
            rend = rstart + 119
            hits.append(make_hit(uid, 'TIR_R', contig, rstart, rend, strand='-'))
            rows.append(
                {
                    'model': 'TIR_R',
                    'target': contig,
                    'hitStart': str(rstart),
                    'hitEnd': str(rend),
                    'strand': '-',
                    'hmmStart': '1',
                    'hmmEnd': '120',
                    'evalue': f'1e-{12 + i}',
                }
            )
            pairs.append({uid - 1, uid})
            uid += 1

    table = hit_table_factory(rows)
    acc = PairReportAccumulator(
        hit_table=table,
        model_lengths={'TIR_L': 120, 'TIR_R': 120},
        contig_length=lambda name: 500_000,
        filter_stats={'mincov': 0.5},
    )
    index = make_index(hits)
    by_uid = {h.idx: h for h in hits}
    elements = [
        make_element(by_uid[min(p)], by_uid[max(p)], element_id=f'E{i}')
        for i, p in enumerate(pairs)
    ]
    acc.add_group(
        left_feature='TIR_L',
        right_feature='TIR_R',
        config=PairingConfig(
            orientation='F,R', left_model='TIR_L', right_model='TIR_R'
        ),
        paired={'TIR_L': pairs},
        hit_index=index,
        elements={'TIR_L': elements},
    )
    acc.add_unpaired(index)
    return acc.finalise()


class TestNatureRcParams:
    def test_carries_the_house_style(self):
        params = nature_rcparams()
        assert params['font.size'] == 7
        assert params['axes.spines.top'] is False
        assert params['axes.spines.right'] is False
        assert params['axes.grid'] is False
        assert params['legend.frameon'] is False
        assert params['axes.linewidth'] == 0.5

    def test_lists_font_fallbacks(self):
        # CI machines rarely have Helvetica; without a fallback every figure
        # emits a font warning.
        assert 'DejaVu Sans' in nature_rcparams()['font.sans-serif']


class TestNamespaceSvgIds:
    def test_rewrites_ids_and_their_references(self):
        svg = (
            '<svg><defs><path id="pAbC" d="M0,0"/></defs>'
            '<use xlink:href="#pAbC"/><g clip-path="url(#pAbC)"/></svg>'
        )
        out = namespace_svg_ids(svg, 'fig1')
        assert 'id="fig1-pAbC"' in out
        assert 'xlink:href="#fig1-pAbC"' in out
        assert 'url(#fig1-pAbC)' in out
        assert '"#pAbC"' not in out

    def test_leaves_external_references_alone(self):
        svg = '<svg><use href="#not-mine"/></svg>'
        assert namespace_svg_ids(svg, 'fig1') == svg

    def test_no_ids_is_a_no_op(self):
        svg = '<svg><rect width="1" height="1"/></svg>'
        assert namespace_svg_ids(svg, 'fig1') == svg


class TestBuildFigures:
    def test_builds_the_expected_figures(self, report):
        figures = build_figures(report)
        ids = [figure.id for figure in figures]
        assert 'element-lengths' in ids
        assert 'pairing-outcome' in ids
        assert 'model-coverage' in ids
        assert 'evalues' in ids

    def test_each_figure_is_inline_ready_svg(self, report):
        for figure in build_figures(report):
            assert figure.svg.startswith('<svg')
            assert '<?xml' not in figure.svg
            assert '<!DOCTYPE' not in figure.svg
            assert figure.title
            assert figure.caption

    def test_figures_in_one_document_have_disjoint_ids(self, report):
        # Two matplotlib SVGs in one page silently blank each other when their
        # internal clip-path and glyph ids collide.
        figures = build_figures(report)
        assert len(figures) >= 2
        seen = set()
        for figure in figures:
            ids = set(ID_ATTR.findall(figure.svg))
            assert not (ids & seen), f'{figure.id} reuses ids from another figure'
            seen |= ids

    def test_global_rcparams_are_not_mutated(self, report):
        # TIRmite is importable as a library and must not reconfigure its
        # host's plotting.
        before = dict(plt.rcParams)
        build_figures(report)
        after = dict(plt.rcParams)
        assert before == after

    def test_mincov_is_drawn_on_the_coverage_figure(self, report):
        figure = next(f for f in build_figures(report) if f.id == 'model-coverage')
        assert '--mincov' in figure.svg or 'mincov' in figure.svg

    def test_caption_distinguishes_the_two_coverage_measures(self, report):
        figure = next(f for f in build_figures(report) if f.id == 'model-coverage')
        assert 'span-based' in figure.caption


class TestLegendPlacement:
    def test_legend_is_anchored_above_the_axes(self):
        # A legend inside the axes lands on the data -- in a distribution plot,
        # exactly where the tallest bars are.
        from tirmite.report.figures import _legend_outside

        with plt.rc_context(nature_rcparams()):
            fig, ax = plt.subplots()
            ax.plot([0, 1], [0, 1], label='series')
            _legend_outside(ax, 1)
            anchor = ax.get_legend().get_bbox_to_anchor()
            in_axes = anchor.transformed(ax.transAxes.inverted())
            plt.close(fig)

        assert in_axes.y0 >= 1.0

    def test_no_figure_places_a_legend_inside_the_axes(self, report):
        # Guards the whole set: a new figure added with a bare ax.legend()
        # would put its key back over the data.
        import tirmite.report.figures as figures_module

        placements = []
        real_legend = plt.Axes.legend

        def record(self, *args, **kwargs):
            placements.append(kwargs.get('bbox_to_anchor'))
            return real_legend(self, *args, **kwargs)

        plt.Axes.legend = record
        try:
            figures_module.build_figures(report)
        finally:
            plt.Axes.legend = real_legend

        assert placements, 'expected at least one legend'
        assert all(anchor is not None for anchor in placements)


class TestDegenerateInputs:
    def test_empty_report_builds_no_figures_and_does_not_raise(self):
        from tirmite.report.model import ReportData

        assert build_figures(ReportData()) == []

    def test_report_without_elements_skips_the_length_figure(self, hit_table_factory):
        table = hit_table_factory(
            [{'model': 'TIR', 'target': 'chr1', 'hitStart': '100', 'hitEnd': '200'}]
        )
        acc = PairReportAccumulator(
            hit_table=table,
            model_lengths={'TIR': 100},
            contig_length=lambda name: 10_000,
        )
        acc.add_unpaired(make_index([make_hit(0, 'TIR', 'chr1', 100, 200)]))
        ids = [figure.id for figure in build_figures(acc.finalise())]
        assert 'element-lengths' not in ids

    def test_unusable_backend_yields_no_figures(self, report, monkeypatch):
        import builtins

        real_import = builtins.__import__

        def fail(name, *args, **kwargs):
            if name == 'matplotlib':
                raise ImportError('no matplotlib here')
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, '__import__', fail)
        assert build_figures(report) == []


class TestRendererIntegration:
    def test_render_embeds_figures(self, report):
        from tirmite.report.render import render_report

        html = render_report(report)
        assert 'Distributions' in html
        assert 'Predicted element lengths' in html
        assert html.count('<svg') >= 2

    def test_figures_are_not_duplicated_into_the_json_payload(self, report):
        import json

        from tirmite.report.figures import build_figures as build
        from tirmite.report.render import render_report

        report.figures = build(report)
        html = render_report(report)
        payload = json.loads(
            re.search(
                r'<script id="tirmite-report-data" type="application/json">(.*?)</script>',
                html,
                re.DOTALL,
            ).group(1)
        )
        # The markup belongs in the body, not in the payload as well.
        assert payload['figures']
        assert all('svg' not in entry for entry in payload['figures'])


class TestPerModelColours:
    def _coverage_series(self, n_models, hit_table_factory):
        from tirmite.report.figures import _model_series

        return _model_series({f'M{i}': [0.5] * (n_models - i) for i in range(n_models)})

    def test_every_model_gets_a_distinct_colour(self, hit_table_factory):
        # A repeated hue in one chart would claim two models are one series.
        series = self._coverage_series(8, hit_table_factory)
        colours = [colour for _, _, colour in series]
        assert len(set(colours)) == len(colours) == 8

    def test_models_past_the_palette_fold_into_one_neutral_series(
        self, hit_table_factory
    ):
        series = self._coverage_series(12, hit_table_factory)
        assert len(series) == 9
        label, values, colour = series[-1]
        assert label == 'other models (n=4)'
        # The fold is stated in the legend rather than hidden.
        assert colour not in {c for _, _, c in series[:-1]}
        assert len(values) == sum(range(1, 5))

    def test_busiest_models_are_kept(self, hit_table_factory):
        series = self._coverage_series(12, hit_table_factory)
        assert series[0][0] == 'M0'
