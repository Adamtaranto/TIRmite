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

ID_ATTR = re.compile(r'\bid="([^"]+)"')


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
        assert 'No cluster-level cross-matches were detected' in html
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
            # One hit matches only model positions 11-60, so there is a real
            # deficit for the padding tests to act on.
            hmm_from = 11 if i == 3 else 1
            length = 60 - hmm_from + 1
            lines.append(
                ' '.join(
                    [
                        'chr1',
                        '-',
                        'TIR_L' if i % 2 == 0 else 'TIR_R',
                        '-',
                        str(hmm_from),
                        '60',
                        str(start),
                        str(start + length - 1),
                        str(start),
                        str(start + length - 1),
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

        genome = tmp_path / 'genome.fa'
        seq = ('ACGT' * 12_500)[:50_000]
        genome.write_text('>chr1\n' + seq + '\n')

        return {
            'hits': hits,
            'lengths': lengths,
            'clusters': clusters,
            'genome': genome,
        }

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
            report_pad_model=False,
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

    def test_alignment_panels_are_built_with_a_genome(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import main

        outdir = tmp_path / 'out'
        args = self.make_args(
            inputs,
            outdir,
            report=True,
            report_msa='anchor',
            genome=inputs['genome'],
        )
        assert main(args) == 0

        payload = json.loads(
            PAYLOAD_RE.search(
                (outdir / 'test_report.html').read_text(encoding='utf-8')
            ).group(1)
        )
        assert [panel['model'] for panel in payload['msa']] == ['TIR_cluster']
        panel = payload['msa'][0]
        assert panel['aligner'] == 'anchor'
        assert panel['n_rows_shown'] == 6
        assert all(row['seq'] for row in panel['rows'])

    def test_genome_gives_true_contig_lengths_even_without_panels(
        self, inputs, tmp_path
    ):
        from tirmite.cli.ensemble_search import main

        outdir = tmp_path / 'out'
        # Indexing happens for the axes too, not only for the panels.
        args = self.make_args(
            inputs, outdir, report=True, report_msa='off', genome=inputs['genome']
        )
        assert main(args) == 0
        payload = json.loads(
            PAYLOAD_RE.search(
                (outdir / 'test_report.html').read_text(encoding='utf-8')
            ).group(1)
        )
        assert payload['contigs'][0]['length_source'] == 'source'
        assert payload['contigs'][0]['length'] == 50_000

    def test_no_temporary_files_are_left_in_the_output_directory(
        self, inputs, tmp_path
    ):
        from tirmite.cli.ensemble_search import main

        outdir = tmp_path / 'out'
        args = self.make_args(
            inputs,
            outdir,
            report=True,
            report_msa='anchor',
            genome=inputs['genome'],
        )
        assert main(args) == 0
        # MAFFT's scratch space belongs in a temporary directory, not beside
        # the user's results.
        assert sorted(p.name for p in outdir.iterdir()) == [
            'test_hits.tab',
            'test_report.html',
        ]

    def test_report_pad_model_shows_sequence_in_the_panels(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import main

        outdir = tmp_path / 'out'
        args = self.make_args(
            inputs,
            outdir,
            report=True,
            report_msa='anchor',
            genome=inputs['genome'],
            report_pad_model=True,
        )
        assert main(args) == 0
        payload = json.loads(
            PAYLOAD_RE.search(
                (outdir / 'test_report.html').read_text(encoding='utf-8')
            ).group(1)
        )
        # Grey runs only exist when the neighbouring bases were fetched.
        assert any(
            run[2] == 'm'
            for panel in payload['msa']
            for row in panel['rows']
            for run in row['pad']
        )

    def test_panels_use_gaps_without_the_flag(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import main

        outdir = tmp_path / 'out'
        args = self.make_args(
            inputs,
            outdir,
            report=True,
            report_msa='anchor',
            genome=inputs['genome'],
        )
        assert main(args) == 0
        payload = json.loads(
            PAYLOAD_RE.search(
                (outdir / 'test_report.html').read_text(encoding='utf-8')
            ).group(1)
        )
        assert not any(
            run[2] == 'm'
            for panel in payload['msa']
            for row in panel['rows']
            for run in row['pad']
        )

    def test_pad_model_without_panels_is_rejected(self, inputs, tmp_path):
        from tirmite.cli.ensemble_search import (
            EnsembleSearchError,
            _validate_search_args,
        )

        # The panels are off by default here, so the flag alone does nothing.
        args = self.make_args(
            inputs,
            tmp_path / 'out',
            report=True,
            report_pad_model=True,
            report_msa='off',
        )
        with pytest.raises(EnsembleSearchError, match='--report-pad-model'):
            _validate_search_args(args)

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


class TestMultiFastaSource:
    """A run over several genomes must be able to read every hit."""

    def make_genome(self, tmp_path, name, contigs):
        from tirmite.utils.utils import indexGenome

        path = tmp_path / name
        path.write_text(''.join(f'>{sid}\n{seq}\n' for sid, seq in contigs.items()))
        genome, _ = indexGenome(path)
        return genome

    def test_dispatches_by_sequence_name(self, tmp_path):
        from tirmite.utils.extract import MultiFastaSource

        a = self.make_genome(tmp_path, 'a.fa', {'chrA': 'ACGT' * 10})
        b = self.make_genome(tmp_path, 'b.fa', {'chrB': 'TTTT' * 10})
        source = MultiFastaSource([a, b])

        assert source.contig_length('chrA') == 40
        assert source.contig_length('chrB') == 40
        assert source.fetch_raw('chrA', 1, 4) == 'ACGT'
        assert source.fetch_raw('chrB', 1, 4) == 'TTTT'

    def test_unknown_sequence_returns_none(self, tmp_path):
        from tirmite.utils.extract import MultiFastaSource

        a = self.make_genome(tmp_path, 'a.fa', {'chrA': 'ACGT'})
        source = MultiFastaSource([a])
        assert source.contig_length('nope') is None
        assert source.fetch_raw('nope', 1, 2) is None

    def test_duplicate_name_takes_the_first_genome_and_warns(self, tmp_path, caplog):
        from tirmite.utils.extract import MultiFastaSource

        a = self.make_genome(tmp_path, 'a.fa', {'shared': 'AAAA'})
        b = self.make_genome(tmp_path, 'b.fa', {'shared': 'CCCC'})
        source = MultiFastaSource([a, b])
        assert source.fetch_raw('shared', 1, 4) == 'AAAA'
        assert 'appears in more than one genome' in caplog.text

    def test_works_through_fetch_region_padded(self, tmp_path):
        from tirmite.utils.extract import MultiFastaSource, fetch_region_padded

        a = self.make_genome(tmp_path, 'a.fa', {'chrA': 'ACGTACGTAC'})
        source = MultiFastaSource([a])
        # Off the start, so the region is padded rather than short: this is the
        # path the alignment panels take.
        region = fetch_region_padded(source, 'chrA', -1, 4, pad_char='-')
        assert region.seq == '--ACGT'
        assert region.left_pad == 2


class TestReportGenomeResolution:
    def test_single_genome(self, tmp_path):
        import argparse

        from tirmite.cli.ensemble_search import _report_genome_paths

        genome = tmp_path / 'g.fa'
        genome.write_text('>chr1\nACGT\n')
        args = argparse.Namespace(genome=genome, genome_list=None)
        assert _report_genome_paths(args) == [genome]

    def test_genome_list_wins_and_is_expanded(self, tmp_path):
        import argparse

        from tirmite.cli.ensemble_search import _report_genome_paths

        first = tmp_path / 'a.fa'
        second = tmp_path / 'b.fa'
        for path in (first, second):
            path.write_text('>chr\nACGT\n')
        listing = tmp_path / 'genomes.txt'
        listing.write_text(f'{first}\n{second}\n')

        args = argparse.Namespace(genome=tmp_path / 'ignored.fa', genome_list=listing)
        assert _report_genome_paths(args) == [first, second]

    def test_missing_files_are_dropped(self, tmp_path):
        import argparse

        from tirmite.cli.ensemble_search import _report_genome_paths

        args = argparse.Namespace(genome=tmp_path / 'absent.fa', genome_list=None)
        assert _report_genome_paths(args) == []

    def test_no_genome_gives_no_source(self):
        from tirmite.cli.ensemble_search import _index_report_source

        assert _index_report_source([]) is None

    def test_several_genomes_give_a_multi_source(self, tmp_path):
        from tirmite.cli.ensemble_search import _index_report_source
        from tirmite.utils.extract import MultiFastaSource

        paths = []
        for name, sid in (('a.fa', 'chrA'), ('b.fa', 'chrB')):
            path = tmp_path / name
            path.write_text(f'>{sid}\nACGTACGT\n')
            paths.append(path)

        source = _index_report_source(paths)
        assert isinstance(source, MultiFastaSource)
        assert source.contig_length('chrA') == 8
        assert source.contig_length('chrB') == 8

    def test_unreadable_genome_is_skipped_not_fatal(self, tmp_path):
        from tirmite.cli.ensemble_search import _index_report_source

        good = tmp_path / 'good.fa'
        good.write_text('>chrA\nACGT\n')
        broken = tmp_path / 'broken.fa'
        broken.write_text('not a fasta at all')

        source = _index_report_source([broken, good])
        assert source is not None
        assert source.contig_length('chrA') == 4


class TestCountModelOverlaps:
    def test_counts_pairs_of_models_sharing_a_locus(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        rows = [
            {'model': 'A', 'hitStart': '100', 'hitEnd': '200'},
            {'model': 'B', 'hitStart': '150', 'hitEnd': '250'},
            {'model': 'C', 'hitStart': '900', 'hitEnd': '950'},
        ]
        assert count_model_overlaps(make_table(rows)) == [
            {'a': 'A', 'b': 'B', 'hits': 1}
        ]

    def test_pairs_are_unordered_and_counted_once(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        rows = [
            {'model': 'B', 'hitStart': '100', 'hitEnd': '200'},
            {'model': 'A', 'hitStart': '150', 'hitEnd': '250'},
            {'model': 'A', 'hitStart': '1000', 'hitEnd': '1100'},
            {'model': 'B', 'hitStart': '1050', 'hitEnd': '1150'},
        ]
        # Both loci name the same unordered pair, so one row of two.
        assert count_model_overlaps(make_table(rows)) == [
            {'a': 'A', 'b': 'B', 'hits': 2}
        ]

    def test_same_model_overlaps_are_reported_as_redundancy(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        rows = [
            {'model': 'A', 'hitStart': '100', 'hitEnd': '200'},
            {'model': 'A', 'hitStart': '150', 'hitEnd': '250'},
        ]
        assert count_model_overlaps(make_table(rows)) == [
            {'a': 'A', 'b': 'A', 'hits': 1}
        ]

    def test_different_contigs_never_overlap(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        rows = [
            {'model': 'A', 'target': 'chr1', 'hitStart': '100', 'hitEnd': '200'},
            {'model': 'B', 'target': 'chr2', 'hitStart': '100', 'hitEnd': '200'},
        ]
        assert count_model_overlaps(make_table(rows)) == []

    def test_abutting_hits_do_not_overlap(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        rows = [
            {'model': 'A', 'hitStart': '100', 'hitEnd': '200'},
            {'model': 'B', 'hitStart': '201', 'hitEnd': '300'},
        ]
        assert count_model_overlaps(make_table(rows)) == []

    def test_single_shared_base_counts(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        # Coordinates are inclusive, so sharing one base is an overlap of 1.
        rows = [
            {'model': 'A', 'hitStart': '100', 'hitEnd': '200'},
            {'model': 'B', 'hitStart': '200', 'hitEnd': '300'},
        ]
        assert count_model_overlaps(make_table(rows))[0]['hits'] == 1

    def test_busiest_pair_comes_first(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        rows = [
            {'model': 'X', 'hitStart': '10', 'hitEnd': '99'},
            {'model': 'Y', 'hitStart': '50', 'hitEnd': '150'},
        ]
        for i in range(3):
            base = 1000 + i * 500
            rows.append({'model': 'A', 'hitStart': str(base), 'hitEnd': str(base + 99)})
            rows.append(
                {'model': 'B', 'hitStart': str(base + 50), 'hitEnd': str(base + 149)}
            )
        result = count_model_overlaps(make_table(rows))
        assert result[0] == {'a': 'A', 'b': 'B', 'hits': 3}

    def test_empty_table(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        assert count_model_overlaps(pd.DataFrame(columns=COLUMNS)) == []


class TestMultiClusterWarning:
    def test_model_in_two_clusters_is_warned_about(self):
        acc = SearchReportAccumulator(
            cluster_map={'A': ['shared', 'a1'], 'B': ['shared', 'b1']},
        )
        assert any('more than one cluster' in w for w in acc.warnings)
        assert any('shared in A, B' in w for w in acc.warnings)

    def test_clean_cluster_map_warns_about_nothing(self):
        acc = SearchReportAccumulator(
            cluster_map={'A': ['a1', 'a2'], 'B': ['b1']},
        )
        assert acc.warnings == []

    def test_the_warning_reaches_the_report(self):
        acc = SearchReportAccumulator(
            cluster_map={'A': ['shared'], 'B': ['shared']},
            contig_length=lambda name: 5000,
        )
        data = acc.finalise(make_table([{'model': 'shared'}]))
        assert any('more than one cluster' in w for w in data.warnings)


class TestOverlapHeatmap:
    def build_report(self, overlaps, cluster_map, models=None):
        summary = SearchFilterSummary()
        summary.model_overlaps = overlaps
        summary.hits_per_model_before_merge = dict.fromkeys(models or [], 1)
        acc = SearchReportAccumulator(
            cluster_map=cluster_map, contig_length=lambda name: 5000
        )
        rows = [
            {'model': m, 'hitStart': str(100 + i * 500), 'hitEnd': str(180 + i * 500)}
            for i, m in enumerate(models or ['x'])
        ]
        return acc.finalise(make_table(rows), summary=summary)

    def test_models_sort_by_cluster_then_name(self):
        from tirmite.report.figures import _model_order

        data = self.build_report(
            [{'a': 'b2', 'b': 'a1', 'hits': 3}],
            {'zebra': ['a1', 'a2'], 'alpha': ['b1', 'b2']},
            models=['a1', 'a2', 'b1', 'b2', 'loner'],
        )
        order, blocks = _model_order(data)
        # Cluster 'alpha' precedes 'zebra'; unclustered models come last.
        assert order == ['b1', 'b2', 'a1', 'a2', 'loner']
        assert blocks == [('alpha', 0, 1), ('zebra', 2, 3)]

    def test_cluster_order_is_case_insensitive(self):
        from tirmite.report.figures import _model_order

        data = self.build_report(
            [{'a': 'h1', 'b': 'M1', 'hits': 1}],
            {'Mut_cluster': ['M1'], 'hAT_cluster': ['h1']},
            models=['h1', 'M1'],
        )
        order, _ = _model_order(data)
        # 'hAT_cluster' must not sort after 'Mut_cluster' on case alone.
        assert order == ['h1', 'M1']

    def test_a_model_in_two_clusters_lands_in_the_first(self):
        from tirmite.report.figures import _model_order

        data = self.build_report(
            [{'a': 'shared', 'b': 'a1', 'hits': 1}],
            {'B_cluster': ['shared', 'b1'], 'A_cluster': ['a1']},
            models=['a1', 'b1', 'shared'],
        )
        order, blocks = _model_order(data)
        assert order.index('shared') > order.index('a1')
        assert ('A_cluster', 0, 0) in blocks

    def test_heatmap_is_built_for_a_search_report(self):
        from tirmite.report.figures import build_figures

        data = self.build_report(
            [{'a': 'a1', 'b': 'a2', 'hits': 4}, {'a': 'a1', 'b': 'b1', 'hits': 1}],
            {'A': ['a1', 'a2'], 'B': ['b1']},
            models=['a1', 'a2', 'b1'],
        )
        figure = next(
            (f for f in build_figures(data) if f.id == 'model-overlaps'), None
        )
        assert figure is not None
        assert figure.svg.startswith('<svg')
        assert figure.wide is True
        assert 'Boxes enclose models' in figure.caption

    def test_no_overlaps_means_no_heatmap(self):
        from tirmite.report.figures import build_figures

        data = self.build_report([], {'A': ['a1']}, models=['a1'])
        assert not [f for f in build_figures(data) if f.id == 'model-overlaps']

    def test_pairing_outcome_is_not_built_for_a_search_report(self):
        from tirmite.report.figures import build_figures

        data = self.build_report(
            [{'a': 'a1', 'b': 'a2', 'hits': 2}],
            {'A': ['a1', 'a2']},
            models=['a1', 'a2'],
        )
        # Every search hit is "unpaired", so that chart would say nothing.
        assert not [f for f in build_figures(data) if f.id == 'pairing-outcome']

    def test_table_twin_labels_the_relationship(self):
        from tirmite.report.stats import search_model_overlap_table

        data = self.build_report(
            [
                {'a': 'a1', 'b': 'a2', 'hits': 5},
                {'a': 'a1', 'b': 'b1', 'hits': 2},
                {'a': 'a1', 'b': 'a1', 'hits': 1},
            ],
            {'A': ['a1', 'a2'], 'B': ['b1']},
            models=['a1', 'a2', 'b1'],
        )
        rows = {
            (r['a'], r['b']): r['relation'] for r in search_model_overlap_table(data)
        }
        assert rows[('a1', 'a2')] == 'same cluster (A)'
        assert rows[('a1', 'b1')] == 'different clusters'
        assert rows[('a1', 'a1')] == 'same model (redundant hits)'


class TestWithinPairOverlapsAreExcluded:
    """The two termini of one element are not a cross-match."""

    def rows(self):
        return [
            {'model': 'TIR_L', 'hitStart': '100', 'hitEnd': '200', 'score': '150'},
            {'model': 'TIR_R', 'hitStart': '150', 'hitEnd': '250', 'score': '50'},
        ]

    def test_partners_are_not_counted(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        # Where a left and right model share a short stretch of identity they
        # hit each other's loci; filter_best_model_per_locus resolves that,
        # and counting it here would light up exactly the pairs the heatmap
        # is meant to leave alone.
        assert (
            count_model_overlaps(
                make_table(self.rows()), pairing_map={'TIR_L': 'TIR_R'}
            )
            == []
        )

    def test_the_same_overlap_counts_without_a_pairing_map(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        assert count_model_overlaps(make_table(self.rows())) == [
            {'a': 'TIR_L', 'b': 'TIR_R', 'hits': 1}
        ]

    def test_partnership_is_symmetric(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        # The map names L->R; an R/L encounter must be excluded too.
        assert (
            count_model_overlaps(
                make_table(self.rows()), pairing_map={'TIR_R': 'TIR_L'}
            )
            == []
        )

    def test_a_sequence_of_pairs_is_accepted(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        # expand_pairing_map_to_components returns pairs, not a dict.
        assert (
            count_model_overlaps(
                make_table(self.rows()), pairing_map=[('TIR_L', 'TIR_R')]
            )
            == []
        )

    def test_unrelated_models_still_count(self):
        from tirmite.cli.ensemble_search import count_model_overlaps

        rows = self.rows() + [
            {'model': 'Other', 'hitStart': '120', 'hitEnd': '220'},
        ]
        result = count_model_overlaps(make_table(rows), pairing_map={'TIR_L': 'TIR_R'})
        pairs = {(entry['a'], entry['b']) for entry in result}
        assert ('Other', 'TIR_L') in pairs
        assert ('Other', 'TIR_R') in pairs
        assert ('TIR_L', 'TIR_R') not in pairs

    def test_within_pair_hits_are_still_filtered_by_the_pipeline(self):
        from tirmite.cli.ensemble_search import filter_best_model_per_locus

        # Excluding them from the heatmap does not mean leaving them in the
        # hit table: the stronger model still takes the locus.
        kept = filter_best_model_per_locus(make_table(self.rows()), {'TIR_L': 'TIR_R'})
        assert kept['model'].tolist() == ['TIR_L']


class TestClusteredHeatmap:
    def build_report(self, overlaps, cluster_map, models):
        summary = SearchFilterSummary()
        summary.model_overlaps = overlaps
        summary.hits_per_model_before_merge = dict.fromkeys(models, 20)
        acc = SearchReportAccumulator(
            cluster_map=cluster_map, contig_length=lambda name: 50_000
        )
        rows = [
            {'model': m, 'hitStart': str(100 + i * 500), 'hitEnd': str(180 + i * 500)}
            for i, m in enumerate(models)
        ]
        return acc.finalise(make_table(rows), summary=summary)

    def sample(self):
        models = ['a1', 'a2', 'a3', 'b1', 'b2']
        overlaps = [
            {'a': 'a1', 'b': 'a2', 'hits': 18},
            {'a': 'a1', 'b': 'a3', 'hits': 16},
            {'a': 'a2', 'b': 'a3', 'hits': 17},
            {'a': 'b1', 'b': 'b2', 'hits': 19},
            {'a': 'a1', 'b': 'b1', 'hits': 1},
        ]
        return self.build_report(
            overlaps, {'A': ['a1', 'a2', 'a3'], 'B': ['b1', 'b2']}, models
        )

    def test_both_heatmaps_are_built(self):
        from tirmite.report.figures import build_figures

        ids = [f.id for f in build_figures(self.sample())]
        assert 'model-overlaps' in ids
        assert 'model-overlaps-clustered' in ids
        # The clustered view sits directly after the declared-cluster one, so
        # the two can be read against each other.
        assert ids.index('model-overlaps-clustered') == ids.index('model-overlaps') + 1

    def test_the_clustered_figure_is_wide_inline_svg(self):
        from tirmite.report.figures import build_figures

        figure = next(
            f
            for f in build_figures(self.sample())
            if f.id == 'model-overlaps-clustered'
        )
        assert figure.svg.startswith('<svg')
        assert figure.wide is True
        assert 'average-linkage' in figure.caption

    def test_the_two_heatmaps_have_disjoint_svg_ids(self):
        from tirmite.report.figures import build_figures

        figures = {f.id: f for f in build_figures(self.sample())}
        a = set(ID_ATTR.findall(figures['model-overlaps'].svg))
        b = set(ID_ATTR.findall(figures['model-overlaps-clustered'].svg))
        assert not (a & b)

    def test_two_models_are_too_few_to_cluster(self):
        from tirmite.report.figures import build_figures

        data = self.build_report(
            [{'a': 'a1', 'b': 'b1', 'hits': 4}],
            {'A': ['a1'], 'B': ['b1']},
            ['a1', 'b1'],
        )
        # One join reveals nothing a tree could show.
        ids = [f.id for f in build_figures(data)]
        assert 'model-overlaps' in ids
        assert 'model-overlaps-clustered' not in ids

    def test_no_overlaps_means_no_clustered_figure(self):
        from tirmite.report.figures import build_figures

        data = self.build_report([], {'A': ['a1']}, ['a1', 'a2', 'a3'])
        assert not [
            f for f in build_figures(data) if f.id == 'model-overlaps-clustered'
        ]


TEXT_NODE = re.compile(r'<text[^>]*>(.*?)</text>', re.DOTALL)


def svg_text(svg):
    """Return the readable text nodes of an inline SVG."""
    return [re.sub(r'<[^>]+>', '', node).strip() for node in TEXT_NODE.findall(svg)]


class TestClusterSymbols:
    def test_each_cluster_gets_a_distinct_symbol(self):
        from tirmite.report.figures import _cluster_symbols

        symbols = _cluster_symbols(['a', 'b', 'c'])
        assert len(set(symbols.values())) == 3

    def test_colour_advances_before_shape(self):
        from tirmite.report.figures import _cluster_symbols
        from tirmite.report.palette import GROUP_COLOURS

        clusters = [f'c{i}' for i in range(len(GROUP_COLOURS))]
        symbols = _cluster_symbols(clusters)
        # Within one palette pass the shape is constant and the hue varies.
        assert len({colour for colour, _ in symbols.values()}) == len(GROUP_COLOURS)
        assert len({marker for _, marker in symbols.values()}) == 1

    def test_shape_distinguishes_clusters_past_the_palette(self):
        from tirmite.report.figures import _cluster_symbols
        from tirmite.report.palette import GROUP_COLOURS

        clusters = [f'c{i}' for i in range(len(GROUP_COLOURS) + 2)]
        symbols = _cluster_symbols(clusters)
        # The ninth cluster reuses a hue, so its shape has to differ.
        assert len(set(symbols.values())) == len(clusters)
        first, ninth = symbols['c0'], symbols[f'c{len(GROUP_COLOURS)}']
        assert first[0] == ninth[0]
        assert first[1] != ninth[1]


class TestClusteredHeatmapAnnotation:
    def build(self, cluster_map, models, overlaps, hits=None):
        summary = SearchFilterSummary()
        summary.model_overlaps = overlaps
        summary.hits_per_model_before_merge = hits or dict.fromkeys(models, 20)
        acc = SearchReportAccumulator(
            cluster_map=cluster_map, contig_length=lambda name: 50_000
        )
        rows = [
            {'model': m, 'hitStart': str(100 + i * 500), 'hitEnd': str(180 + i * 500)}
            for i, m in enumerate(models)
        ]
        return acc.finalise(make_table(rows), summary=summary)

    def figure(self, data):
        from tirmite.report.figures import build_figures

        return next(
            f for f in build_figures(data) if f.id == 'model-overlaps-clustered'
        )

    def sample(self):
        models = ['a1', 'a2', 'b1', 'b2']
        overlaps = [
            {'a': 'a1', 'b': 'a2', 'hits': 18},
            {'a': 'b1', 'b': 'b2', 'hits': 19},
            {'a': 'a1', 'b': 'b1', 'hits': 1},
        ]
        return self.build(
            {'Alpha': ['a1', 'a2'], 'beta': ['b1', 'b2']},
            models,
            overlaps,
            hits={'a1': 20, 'a2': 21, 'b1': 22, 'b2': 23},
        )

    def test_the_key_names_every_cluster_present(self):
        text = svg_text(self.figure(self.sample()).svg)
        assert 'Cluster' in text
        assert 'Alpha' in text
        assert 'beta' in text

    def test_row_labels_carry_total_hit_counts(self):
        text = svg_text(self.figure(self.sample()).svg)
        # Totals make an off-diagonal count interpretable; they are not on the
        # diagonal because they would swamp the colour scale.
        assert any('a1' in t and '(20)' in t for t in text)
        assert any('b2' in t and '(23)' in t for t in text)

    def test_column_labels_stay_bare(self):
        text = svg_text(self.figure(self.sample()).svg)
        # Repeating the totals on both axes would be noise.
        assert 'a1' in text

    def test_the_key_is_ordered_case_insensitively(self):
        from tirmite.report.figures import _cluster_symbols

        data = self.sample()
        cluster_map = data.stats['cluster_map']
        present = sorted(cluster_map, key=str.casefold)
        # 'Alpha' before 'beta' despite the leading capital.
        assert present == ['Alpha', 'beta']
        assert list(_cluster_symbols(present)) == present

    def test_unclustered_models_are_labelled_in_the_key(self):
        data = self.build(
            {'Alpha': ['a1', 'a2']},
            ['a1', 'a2', 'loner'],
            [
                {'a': 'a1', 'b': 'a2', 'hits': 18},
                {'a': 'a1', 'b': 'loner', 'hits': 2},
            ],
        )
        assert 'no cluster' in svg_text(self.figure(data).svg)

    def test_no_cluster_map_still_renders(self):
        data = self.build(
            {},
            ['m1', 'm2', 'm3'],
            [{'a': 'm1', 'b': 'm2', 'hits': 5}, {'a': 'm2', 'b': 'm3', 'hits': 3}],
        )
        figure = self.figure(data)
        assert figure.svg.startswith('<svg')
        # Nothing to key, so no cluster legend is drawn.
        assert 'Cluster' not in svg_text(figure.svg)

    def test_the_diagonal_is_self_overlap_not_the_total(self):
        # Two hits of one model at one locus is redundancy; the model's total
        # hit count is a different quantity and stays off this colour scale.
        data = self.build(
            {'Alpha': ['a1', 'a2']},
            ['a1', 'a2', 'a3'],
            [
                {'a': 'a1', 'b': 'a2', 'hits': 4},
                {'a': 'a2', 'b': 'a3', 'hits': 3},
                {'a': 'a1', 'b': 'a1', 'hits': 2},
            ],
            hits={'a1': 500, 'a2': 500, 'a3': 500},
        )
        text = svg_text(self.figure(data).svg)
        assert '2' in text
        # The totals appear only as labels, never as a cell value.
        assert '500' not in text
        assert any('(500)' in t for t in text)
