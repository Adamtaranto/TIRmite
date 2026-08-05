"""Tests for the `tirmite search` HTML report."""

import json
import re

import pandas as pd
import pytest

from tirmite.cli.ensemble_search import (
    SearchFilterSummary,
    check_cross_cluster_overlaps,
)
from tirmite.report.render import render_report, write_search_report
from tirmite.report.search import SearchReportAccumulator
from tirmite.report.stats import (
    search_cluster_table,
    search_contested_table,
    search_cross_match_matrix,
    search_cross_match_table,
    search_query_table,
    search_stage_table,
)

COLUMNS = [
    'model',
    'target',
    'hitStart',
    'hitEnd',
    'strand',
    'evalue',
    'score',
    'bias',
    'hmmStart',
    'hmmEnd',
]

PAYLOAD_RE = re.compile(
    r'<script id="tirmite-report-data" type="application/json">(.*?)</script>',
    re.DOTALL,
)


def make_table(rows):
    """Build a string-valued hit table from partial row dicts."""
    defaults = {
        'model': 'TIR',
        'target': 'chr1',
        'hitStart': '100',
        'hitEnd': '200',
        'strand': '+',
        'evalue': '1e-20',
        'score': '100',
        'bias': '0.1',
        'hmmStart': '1',
        'hmmEnd': '80',
    }
    return pd.DataFrame([{**defaults, **row} for row in rows], columns=COLUMNS).astype(
        str
    )


@pytest.fixture
def summary():
    s = SearchFilterSummary()
    s.stages = [('Hits loaded', 100), ('E-value filter', 60), ('Cluster merging', 50)]
    s.cross_cluster_overlaps = [
        {
            'contig': 'chr1',
            'a_cluster': 'A',
            'a_start': 100,
            'a_end': 200,
            'b_cluster': 'B',
            'b_start': 150,
            'b_end': 260,
            'overlap': 51,
        },
        {
            'contig': 'chr2',
            'a_cluster': 'B',
            'a_start': 10,
            'a_end': 90,
            'b_cluster': 'A',
            'b_start': 50,
            'b_end': 130,
            'overlap': 41,
        },
    ]
    s.hits_per_model_before_merge = {'A_left': 30, 'A_right': 25, 'B_tir': 45}
    s.cross_model_removed = {('B', 'A'): 4}
    s.nested_removed = {'B': {'A': 2}}
    s.excluded_not_in_map = {'orphan': 9}
    return s


@pytest.fixture
def report(summary):
    table = make_table(
        [
            {'model': 'A', 'hitStart': '100', 'hitEnd': '179', 'hmmEnd': '80'},
            {'model': 'A', 'target': 'chr2', 'hitStart': '500', 'hitEnd': '579'},
            {'model': 'B', 'hitStart': '150', 'hitEnd': '229', 'hmmStart': '11'},
        ]
    )
    acc = SearchReportAccumulator(
        tirmite_version='9.9.9',
        command='tirmite search --report',
        title='Search test',
        model_lengths={'A_left': 80, 'A_right': 80, 'B_tir': 80},
        contig_length=lambda name: 10_000,
        cluster_map={'A': ['A_left', 'A_right'], 'B': ['B_tir']},
        pairing_map={'A': 'B'},
    )
    return acc.finalise(table, summary=summary)


class TestCrossClusterOverlapsReturnValue:
    def test_returns_every_overlap_it_finds(self):
        # It logs only the first ten, but the report needs them all.
        rows = []
        for i in range(15):
            base = 1000 + i * 1000
            rows.append(
                {'model': 'a1', 'hitStart': str(base), 'hitEnd': str(base + 99)}
            )
            rows.append(
                {'model': 'b1', 'hitStart': str(base + 50), 'hitEnd': str(base + 149)}
            )
        overlaps = check_cross_cluster_overlaps(
            make_table(rows), {'A': ['a1'], 'B': ['b1']}
        )
        assert len(overlaps) == 15
        assert overlaps[0]['overlap'] == 50
        assert {o['a_cluster'] for o in overlaps} == {'A'}
        assert {o['b_cluster'] for o in overlaps} == {'B'}

    def test_empty_table_returns_empty_list(self):
        empty = pd.DataFrame(columns=COLUMNS)
        assert check_cross_cluster_overlaps(empty, {'A': ['a1']}) == []

    def test_same_cluster_overlaps_are_not_reported(self):
        rows = [
            {'model': 'a1', 'hitStart': '100', 'hitEnd': '200'},
            {'model': 'a2', 'hitStart': '150', 'hitEnd': '250'},
        ]
        assert check_cross_cluster_overlaps(make_table(rows), {'A': ['a1', 'a2']}) == []


class TestSearchReportData:
    def test_shape(self, report):
        assert report.kind == 'search'
        assert report.hits.n == 3
        assert report.elements == []
        assert {g.group_id for g in report.groups} == {'A', 'B'}

    def test_groups_are_clusters_and_carry_colours(self, report):
        for group in report.groups:
            assert group.colour.startswith('#')
            assert group.colour_dark.startswith('#')
        assert len({g.colour for g in report.groups}) == 2

    def test_busiest_group_takes_the_first_palette_slot(self, report):
        assert report.groups[0].group_id == 'A'
        assert report.groups[0].n_unpaired == 2

    def test_cluster_lengths_are_resolved_from_components(self, report):
        # Hits carry cluster names after merging, but the lengths file is
        # keyed by component, so coverage would otherwise be unavailable.
        assert {m.name: m.length for m in report.models} == {'A': 80, 'B': 80}
        assert all(cov is not None for cov in report.hits.model_cov)

    def test_ambiguous_cluster_lengths_are_left_unset_with_a_warning(self):
        acc = SearchReportAccumulator(
            model_lengths={'x1': 60, 'x2': 90},
            cluster_map={'X': ['x1', 'x2']},
        )
        # Components disagree, so there is no single denominator to measure a
        # merged hit against; inventing one would misstate coverage.
        assert 'X' not in acc.model_lengths
        assert any('different lengths' in w for w in acc.warnings)

    def test_explicit_cluster_length_wins_over_components(self):
        acc = SearchReportAccumulator(
            model_lengths={'X': 100, 'x1': 60, 'x2': 90},
            cluster_map={'X': ['x1', 'x2']},
        )
        assert acc.model_lengths['X'] == 100
        assert acc.warnings == []

    def test_truncation_still_drives_jagged_edges(self, report):
        # The B hit starts at model position 11, so 10 positions are unmatched.
        index = list(report.hits.uid).index(2)
        assert report.hits.trunc_l[index] == 10

    def test_hits_are_stacked(self, report):
        assert all(row >= 0 for row in report.hits.row)

    def test_stats_payload_carries_the_search_sections(self, report):
        stats = report.stats
        assert [s['label'] for s in stats['stages']][0] == 'Hits loaded'
        assert len(stats['cross_matches']) == 2
        assert stats['cluster_map'] == {'A': ['A_left', 'A_right'], 'B': ['B_tir']}
        assert stats['hits_per_model_before_merge']['A_left'] == 30
        assert stats['cross_model_removed'] == [
            {'removed': 'B', 'kept': 'A', 'hits': 4}
        ]
        assert stats['nested_removed'] == {'B': {'A': 2}}

    def test_json_round_trips(self, report):
        payload = json.loads(report.to_json())
        assert payload['kind'] == 'search'
        assert payload['elements'] == []

    def test_no_cluster_map_groups_by_model(self):
        table = make_table([{'model': 'TIR_A'}, {'model': 'TIR_B', 'hitStart': '900'}])
        acc = SearchReportAccumulator(contig_length=lambda name: 5000)
        data = acc.finalise(table)
        assert {g.group_id for g in data.groups} == {'TIR_A', 'TIR_B'}

    def test_hit_cap_keeps_the_best_and_warns(self):
        rows = [
            {
                'model': 'TIR',
                'hitStart': str(i * 100 + 1),
                'hitEnd': str(i * 100 + 50),
                'evalue': f'1e-{i}',
            }
            for i in range(1, 11)
        ]
        acc = SearchReportAccumulator(contig_length=lambda name: 5000, max_hits=3)
        data = acc.finalise(make_table(rows))
        assert data.hits.n == 3
        assert any('omitted to keep the report' in w for w in data.warnings)

    def test_empty_hit_table_is_tolerated(self):
        acc = SearchReportAccumulator()
        data = acc.finalise(pd.DataFrame(columns=COLUMNS))
        assert data.hits.n == 0
        assert data.contigs == []


class TestSearchTables:
    def test_query_table_sums_component_counts_for_a_cluster(self, report):
        rows = {row['model']: row for row in search_query_table(report)}
        assert rows['A']['before_merge']['sort'] == 55
        assert rows['B']['before_merge']['sort'] == 45

    def test_stage_table_shows_the_change_at_each_step(self, report):
        rows = search_stage_table(report)
        assert rows[0]['change']['value'] == '–'
        assert rows[1]['change']['value'] == '-40'
        assert rows[2]['change']['value'] == '-10'

    def test_cross_match_table_is_ordered_by_overlap(self, report):
        rows = search_cross_match_table(report)
        assert [row['overlap']['sort'] for row in rows] == [51, 41]

    def test_cross_match_matrix_counts_each_pair_once(self, report):
        matrix = search_cross_match_matrix(report)
        assert matrix['labels'] == ['A', 'B']
        # Both records name the same unordered pair, in opposite orders.
        assert matrix['counts'] == {'A\tB': 2}

    def test_cluster_table_reports_what_merging_collapsed(self, report):
        rows = {row['cluster']['value']: row for row in search_cluster_table(report)}
        assert rows['A']['n_components']['sort'] == 2
        assert rows['A']['before']['sort'] == 55
        assert rows['A']['components'] == 'A_left, A_right'

    def test_contested_table_covers_all_three_filter_steps(self, report):
        steps = {row['step'] for row in search_contested_table(report)}
        assert steps == {
            'Cross-model overlap',
            'Nested in partner',
            'Not in pairing map',
        }


class TestSearchRendering:
    def test_renders_a_self_contained_document(self, report):
        html = render_report(report, template='search_report.html.j2')
        assert html.startswith('<!DOCTYPE html>')
        assert '<script src=' not in html
        assert 'rel="stylesheet"' not in html

    def test_carries_the_search_sections(self, report):
        html = render_report(report, template='search_report.html.j2')
        for heading in (
            'Query termini',
            'Cross-matches',
            'Clustering',
            'All hits',
            'Contested loci per pair of groups',
        ):
            assert heading in html

    def test_pair_only_sections_are_absent(self, report):
        html = render_report(report, template='search_report.html.j2')
        assert 'Predicted elements' not in html
        assert 'Pairing groups' not in html
        assert 'id="elements-table"' not in html

    def test_payload_matches_the_data(self, report):
        html = render_report(report, template='search_report.html.j2')
        payload = json.loads(PAYLOAD_RE.search(html).group(1))
        assert payload == json.loads(report.to_json())

    def test_empty_cross_matches_show_an_explanation(self):
        table = make_table([{'model': 'TIR'}])
        acc = SearchReportAccumulator(contig_length=lambda name: 5000)
        html = render_report(
            acc.finalise(table, summary=SearchFilterSummary()),
            template='search_report.html.j2',
        )
        assert 'No cross-matches were detected' in html
        assert 'No cluster map was used' in html

    def test_write_search_report_creates_parent_directories(self, report, tmp_path):
        target = tmp_path / 'nested' / 'search.html'
        assert write_search_report(report, target) == target
        assert target.read_text(encoding='utf-8').startswith('<!DOCTYPE html>')


class TestSearchCli:
    """End-to-end: `tirmite search --report` over precomputed nhmmer results."""

    @pytest.fixture
    def inputs(self, tmp_path):
        header = '# target accession query accession hmmfrom hmmto ...'
        lines = [header]
        for i in range(6):
            start = 1000 + i * 5000
            lines.append(
                ' '.join(
                    [
                        'chr1',
                        '-',
                        'TIR_L' if i % 2 == 0 else 'TIR_R',
                        '-',
                        '1',
                        '60',
                        str(start),
                        str(start + 59),
                        str(start),
                        str(start + 59),
                        '50000',
                        '+',
                        '1e-25',
                        '95.0',
                        '0.1',
                        '-',
                    ]
                )
            )
        hits = tmp_path / 'hits.tbl'
        hits.write_text('\n'.join(lines) + '\n')

        lengths = tmp_path / 'lengths.txt'
        lengths.write_text('TIR_L\t60\nTIR_R\t60\n')

        clusters = tmp_path / 'clusters.tsv'
        clusters.write_text('TIR_cluster\tTIR_L\tTIR_R\n')

        return {'hits': hits, 'lengths': lengths, 'clusters': clusters}

    def make_args(self, inputs, outdir, **overrides):
        import argparse

        args = argparse.Namespace(
            fasta=None,
            hmm=None,
            genome=None,
            genome_list=None,
            blast_db=None,
            blast_results=None,
            nhmmer_results=[inputs['hits']],
            lengths_file=inputs['lengths'],
            cluster_map=inputs['clusters'],
            pairing_map=None,
            max_evalue=1e-3,
            blast_max_evalue=1e-3,
            hmm_max_evalue=1e-3,
            min_identity=60,
            word_size=4,
            max_offset=None,
            orientation='F,R',
            min_score_ratio=1.5,
            merge_max_gap=1,
            outdir=outdir,
            prefix='test',
            split_paired_output=False,
            threads=1,
            blast_timeout=None,
            keep_temp=False,
            loglevel='ERROR',
            logfile=None,
            report=False,
            report_out=None,
            report_title=None,
            report_msa='off',
            report_msa_max_rows=500,
            report_max_hits=200000,
            report_max_rows=30,
        )
        for key, value in overrides.items():
            setattr(args, key, value)
        return args

    def test_report_is_written(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import main

        outdir = tmp_path / 'out'
        assert main(self.make_args(inputs, outdir, report=True)) == 0

        report = outdir / 'test_report.html'
        assert report.exists()
        payload = json.loads(
            PAYLOAD_RE.search(report.read_text(encoding='utf-8')).group(1)
        )
        assert payload['kind'] == 'search'
        assert payload['hits']['n'] == 6
        assert payload['groups'][0]['group_id'] == 'TIR_cluster'
        # The stage table records the whole pipeline, not just the end.
        assert [s['label'] for s in payload['stats']['stages']][0] == 'Hits loaded'
        assert payload['stats']['cluster_map'] == {'TIR_cluster': ['TIR_L', 'TIR_R']}

    def test_no_report_flag_writes_nothing(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import main

        outdir = tmp_path / 'out'
        assert main(self.make_args(inputs, outdir)) == 0
        assert not list(outdir.glob('*.html'))

    def test_hits_table_is_unchanged_by_the_report(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import main

        without = tmp_path / 'without'
        with_report = tmp_path / 'with'
        main(self.make_args(inputs, without))
        main(self.make_args(inputs, with_report, report=True))
        assert (with_report / 'test_hits.tab').read_bytes() == (
            without / 'test_hits.tab'
        ).read_bytes()

    def test_a_failing_report_does_not_fail_the_run(
        self, inputs, tmp_path, monkeypatch
    ):
        from tirmite.cli.ensemble_search import main
        import tirmite.report.render as render_module

        def boom(*args, **kwargs):
            raise RuntimeError('renderer exploded')

        monkeypatch.setattr(render_module, 'write_search_report', boom)

        outdir = tmp_path / 'out'
        # The hit table is already on disk; a broken visualisation must not
        # turn a completed search into a failure.
        assert main(self.make_args(inputs, outdir, report=True)) == 0
        assert (outdir / 'test_hits.tab').exists()
        assert not (outdir / 'test_report.html').exists()

    def test_report_out_implies_report(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import _validate_search_args, main

        target = tmp_path / 'elsewhere' / 'custom.html'
        args = self.make_args(inputs, tmp_path / 'out', report_out=target)
        _validate_search_args(args)
        assert args.report is True
        assert main(args) == 0
        assert target.exists()

    def test_msa_without_a_genome_is_rejected(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import (
            EnsembleSearchError,
            _validate_search_args,
        )

        args = self.make_args(
            inputs, tmp_path / 'out', report=True, report_msa='anchor'
        )
        with pytest.raises(EnsembleSearchError, match='--report-msa requires'):
            _validate_search_args(args)

    @pytest.mark.parametrize(
        'field', ['report_max_hits', 'report_max_rows', 'report_msa_max_rows']
    )
    def test_out_of_range_caps_are_rejected(self, inputs, tmp_path, field):
        from tirmite.cli.ensemble_search import (
            EnsembleSearchError,
            _validate_search_args,
        )

        args = self.make_args(inputs, tmp_path / 'out', report=True, **{field: 0})
        with pytest.raises(EnsembleSearchError, match=field.replace('_', '-')):
            _validate_search_args(args)
