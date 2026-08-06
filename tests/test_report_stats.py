"""
Tests for tirmite.report.stats.

The golden-string tests here exist to guard a refactor: the counting and the
formatting of the pair summary were split out of
``hmm_pair._write_pair_summary`` into this module, and the ``*_summary.txt``
file must remain byte-identical.
"""

from tirmite.cli.hmm_pair import _write_pair_summary
from tirmite.report.stats import (
    PairSummary,
    format_pair_summary,
    pair_summary_stats,
)


def make_hit_index(counts):
    """Build a minimal hitIndex with `counts[model]` entries per model."""
    index = {}
    uid = 0
    for model, n in counts.items():
        index[model] = {}
        for _ in range(n):
            index[model][uid] = {'rec': None, 'partner': None, 'candidates': []}
            uid += 1
    return index


FULL_FILTER_STATS = {
    'initial_hits': 100,
    'mincov': 0.5,
    'coverage_excluded': 7,
    'maxeval': 0.001,
    'evalue_excluded': 3,
    'max_offset': 10,
    'anchor_excluded': 2,
    'pairing_map_models_ignored': ['ZZZ', 'AAA'],
    'pairing_map_hits_ignored': 5,
    'after_filtering': 83,
}


class TestPairSummaryStats:
    def test_counts_hits_per_model(self):
        summary = pair_summary_stats(
            left_feature='TIR_L',
            right_feature='TIR_R',
            pair_paired={},
            pair_hitIndex=make_hit_index({'TIR_L': 4, 'TIR_R': 6}),
            total_pairs=0,
            total_elements=0,
        )
        assert summary.hits_per_model == {'TIR_L': 4, 'TIR_R': 6}
        assert summary.total_hits == 10

    def test_model_absent_from_index_counts_zero(self):
        # A zero row must be present rather than silently missing, so the
        # report shows that one side of the pairing found nothing.
        summary = pair_summary_stats(
            left_feature='TIR_L',
            right_feature='TIR_R',
            pair_paired={},
            pair_hitIndex=make_hit_index({'TIR_L': 3}),
            total_pairs=0,
            total_elements=0,
        )
        assert summary.hits_per_model == {'TIR_L': 3, 'TIR_R': 0}
        assert summary.total_hits == 3

    def test_unpaired_is_total_minus_distinct_paired_ids(self):
        summary = pair_summary_stats(
            left_feature='TIR_L',
            right_feature='TIR_R',
            pair_paired={'TIR_L': [{0, 4}, {1, 5}]},
            pair_hitIndex=make_hit_index({'TIR_L': 4, 'TIR_R': 4}),
            total_pairs=2,
            total_elements=2,
        )
        assert summary.total_hits == 8
        assert summary.total_unpaired == 4

    def test_symmetric_pairing_labels(self):
        summary = pair_summary_stats(
            left_feature='TIR',
            right_feature='TIR',
            pair_paired={},
            pair_hitIndex=make_hit_index({'TIR': 2}),
            total_pairs=0,
            total_elements=0,
        )
        # The human-facing label collapses, but the group id must not: a run
        # can mix a symmetric row and an asymmetric row that share a model.
        assert summary.pair_label == 'TIR'
        assert summary.group_id == 'TIR__TIR'
        assert summary.hits_per_model == {'TIR': 2}

    def test_asymmetric_pairing_labels(self):
        summary = pair_summary_stats(
            left_feature='TIR_L',
            right_feature='TIR_R',
            pair_paired={},
            pair_hitIndex={},
            total_pairs=0,
            total_elements=0,
        )
        assert summary.pair_label == 'TIR_L_TIR_R'
        assert summary.group_id == 'TIR_L__TIR_R'

    def test_filter_stats_are_copied_not_aliased(self):
        stats = dict(FULL_FILTER_STATS)
        summary = pair_summary_stats(
            left_feature='A',
            right_feature='B',
            pair_paired={},
            pair_hitIndex={},
            total_pairs=0,
            total_elements=0,
            filter_stats=stats,
        )
        stats['initial_hits'] = 999
        assert summary.filter_stats['initial_hits'] == 100

    def test_no_filter_stats_gives_empty_dict(self):
        summary = pair_summary_stats(
            left_feature='A',
            right_feature='B',
            pair_paired={},
            pair_hitIndex={},
            total_pairs=0,
            total_elements=0,
        )
        assert summary.filter_stats == {}


class TestFormatPairSummary:
    def test_golden_full(self):
        summary = pair_summary_stats(
            left_feature='TIR_L',
            right_feature='TIR_R',
            pair_paired={'TIR_L': [{0, 4}]},
            pair_hitIndex=make_hit_index({'TIR_L': 4, 'TIR_R': 4}),
            total_pairs=1,
            total_elements=1,
            filter_stats=FULL_FILTER_STATS,
        )
        expected = (
            'TIRmite Pair Summary Report\n'
            '===========================\n'
            '\n'
            'Model pair: TIR_L <-> TIR_R\n'
            '\n'
            'Filtering criteria applied\n'
            '--------------------------\n'
            '  Initial hits imported: 100\n'
            '  Pairing-map model filter: 2 model(s) ignored (AAA, ZZZ), '
            '5 hits excluded\n'
            '  Coverage filter (min coverage >= 0.5): 7 hit(s) excluded\n'
            '  E-value filter (max e-value <= 0.001): 3 hit(s) excluded\n'
            '  Anchor offset filter (max offset <= 10): 2 hit(s) excluded\n'
            '  Hits remaining after all filters: 83\n'
            '\n'
            'Hits per model (after all filters)\n'
            '----------------------------------\n'
            '  TIR_L: 4\n'
            '  TIR_R: 4\n'
            '\n'
            'Total hits for this pair: 8\n'
            'Total pairs found: 1\n'
            'Total paired elements extracted: 1\n'
            'Total unpaired hits: 6\n'
        )
        assert format_pair_summary(summary) == expected

    def test_golden_no_filter_stats(self):
        summary = pair_summary_stats(
            left_feature='TIR',
            right_feature='TIR',
            pair_paired={},
            pair_hitIndex=make_hit_index({'TIR': 2}),
            total_pairs=0,
            total_elements=0,
        )
        expected = (
            'TIRmite Pair Summary Report\n'
            '===========================\n'
            '\n'
            'Model pair: TIR <-> TIR\n'
            '\n'
            'Hits per model (after all filters)\n'
            '----------------------------------\n'
            '  TIR: 2\n'
            '\n'
            # Symmetric pairing counts the single model once, not twice.
            'Total hits for this pair: 2\n'
            'Total pairs found: 0\n'
            'Total paired elements extracted: 0\n'
            'Total unpaired hits: 2\n'
        )
        assert format_pair_summary(summary) == expected

    def test_anchor_excluded_unknown(self):
        # max_offset was applied but the exclusion count was never recorded.
        summary = PairSummary(
            left_feature='A',
            right_feature='B',
            hits_per_model={'A': 1, 'B': 1},
            total_hits=2,
            total_pairs=0,
            total_elements=0,
            total_unpaired=2,
            filter_stats={'max_offset': 5, 'anchor_excluded': None},
        )
        text = format_pair_summary(summary)
        assert 'Anchor offset filter (max offset <= 5): unknown hit(s) excluded' in text

    def test_optional_filter_lines_are_omitted_when_absent(self):
        summary = PairSummary(
            left_feature='A',
            right_feature='B',
            hits_per_model={'A': 1, 'B': 1},
            total_hits=2,
            total_pairs=0,
            total_elements=0,
            total_unpaired=2,
            filter_stats={'initial_hits': 2},
        )
        text = format_pair_summary(summary)
        assert 'Initial hits imported: 2' in text
        assert 'Coverage filter' not in text
        assert 'E-value filter' not in text
        assert 'Anchor offset filter' not in text
        assert 'Pairing-map model filter' not in text
        assert 'Hits remaining after all filters' not in text


class TestWritePairSummaryWrapper:
    def test_writes_formatted_summary_and_returns_it(self, tmp_path):
        summary = _write_pair_summary(
            outdir=str(tmp_path),
            prefix='run1',
            left_feature='TIR_L',
            right_feature='TIR_R',
            hitTable=None,
            pair_paired={'TIR_L': [{0, 4}]},
            pair_hitIndex=make_hit_index({'TIR_L': 4, 'TIR_R': 4}),
            total_pairs=1,
            total_elements=1,
            filter_stats=FULL_FILTER_STATS,
        )
        written = (tmp_path / 'run1_TIR_L_TIR_R_summary.txt').read_text()
        assert written == format_pair_summary(summary)
        assert isinstance(summary, PairSummary)

    def test_filename_without_prefix(self, tmp_path):
        _write_pair_summary(
            outdir=str(tmp_path),
            prefix=None,
            left_feature='TIR',
            right_feature='TIR',
            hitTable=None,
            pair_paired={},
            pair_hitIndex=make_hit_index({'TIR': 1}),
            total_pairs=0,
            total_elements=0,
        )
        assert (tmp_path / 'TIR_summary.txt').exists()
