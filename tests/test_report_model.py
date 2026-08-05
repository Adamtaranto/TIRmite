"""Tests for tirmite.report.model and tirmite.report.collect."""

from collections import namedtuple
import json

import pytest

from tirmite.core.pairing import PairingConfig
from tirmite.report.collect import PairReportAccumulator, model_truncation
from tirmite.report.model import SCHEMA_VERSION, to_json
from tirmite.report.stats import pair_summary_stats

hitTup = namedtuple(
    'hitTup',
    ['model', 'target', 'hitStart', 'hitEnd', 'strand', 'idx', 'evalue'],
)

gffTup = namedtuple(
    'gffTup',
    [
        'model',
        'chromosome',
        'start',
        'end',
        'strand',
        'type',
        'id',
        'leftHit',
        'rightHit',
        'seq',
        'evalue',
    ],
)


def make_hit(idx, model, target, start, end, strand='+', evalue='1e-20'):
    return hitTup(model, target, start, end, strand, idx, evalue)


def make_index(hits):
    """Build a hitIndex from hitTups, keyed by model then row index."""
    index = {}
    for hit in hits:
        index.setdefault(hit.model, {})[hit.idx] = {
            'rec': hit,
            'partner': None,
            'candidates': [],
        }
    return index


def make_element(left, right, element_id='Element_1', seq='ACGTACGT'):
    return gffTup(
        left.model,
        left.target,
        left.hitStart,
        right.hitEnd,
        left.strand,
        'Element',
        element_id,
        left,
        right,
        seq,
        left.evalue,
    )


class TestModelTruncation:
    def test_full_length_hit_has_no_deficit(self):
        assert model_truncation(1, 100, 100, '+') == (0, 0)

    def test_plus_strand_start_deficit_is_at_the_lower_end(self):
        # 10 model positions missing before the alignment starts.
        assert model_truncation(11, 100, 100, '+') == (10, 0)

    def test_plus_strand_end_deficit_is_at_the_higher_end(self):
        assert model_truncation(1, 90, 100, '+') == (0, 10)

    def test_minus_strand_mirrors_the_mapping(self):
        # On the minus strand hmmStart aligns to the higher genomic coordinate.
        assert model_truncation(11, 100, 100, '-') == (0, 10)
        assert model_truncation(1, 90, 100, '-') == (10, 0)

    def test_missing_coordinates_give_none(self):
        assert model_truncation(None, 100, 100, '+') == (None, None)
        assert model_truncation(1, None, 100, '+') == (None, None)
        assert model_truncation(1, 100, None, '+') == (None, None)

    def test_negative_deficit_is_clamped(self):
        # An alignment past the declared model end means a wrong model length,
        # not a negative truncation.
        assert model_truncation(1, 120, 100, '+') == (0, 0)


class TestAccumulatorHits:
    def _accumulate(self, hit_table_factory, **kwargs):
        hits = [
            make_hit(0, 'TIR_L', 'chr1', 100, 200),
            make_hit(1, 'TIR_R', 'chr1', 5000, 5100, strand='-'),
        ]
        table = hit_table_factory(
            [
                {
                    'model': 'TIR_L',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'hmmStart': '11',
                    'hmmEnd': '100',
                    'score': '88.5',
                },
                {
                    'model': 'TIR_R',
                    'target': 'chr1',
                    'hitStart': '5000',
                    'hitEnd': '5100',
                    'strand': '-',
                    'hmmStart': '1',
                    'hmmEnd': '100',
                    'score': '90.0',
                },
            ]
        )
        acc = PairReportAccumulator(
            hit_table=table,
            model_lengths={'TIR_L': 100, 'TIR_R': 100},
            contig_length=lambda name: 100_000,
            **kwargs,
        )
        index = make_index(hits)
        config = PairingConfig(
            orientation='F,R', left_model='TIR_L', right_model='TIR_R'
        )
        summary = pair_summary_stats(
            left_feature='TIR_L',
            right_feature='TIR_R',
            pair_paired={'TIR_L': [{0, 1}]},
            pair_hitIndex=index,
            total_pairs=1,
            total_elements=1,
        )
        acc.add_group(
            left_feature='TIR_L',
            right_feature='TIR_R',
            config=config,
            paired={'TIR_L': [{0, 1}]},
            hit_index=index,
            elements={'TIR_L': [make_element(hits[0], hits[1])]},
            summary=summary,
        )
        return acc

    def test_basic_shape(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        assert data.kind == 'pair'
        assert data.schema_version == SCHEMA_VERSION
        assert data.hits.n == 2
        assert len(data.groups) == 1
        assert len(data.contigs) == 1
        assert len(data.elements) == 1

    def test_model_coverage_and_truncation(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        # Hits are sorted by (contig, start), so TIR_L is index 0.
        assert data.hits.hmm_start[0] == 11
        assert data.hits.model_cov[0] == pytest.approx(0.9)
        # Plus strand, 10 model positions missing before the alignment.
        assert data.hits.trunc_l[0] == 10
        assert data.hits.trunc_r[0] == 0
        # Minus strand hit fully covers its model.
        assert data.hits.trunc_l[1] == 0
        assert data.hits.trunc_r[1] == 0

    def test_score_is_read_from_the_hit_table(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        assert data.hits.score[0] == pytest.approx(88.5)

    def test_strand_is_packed_as_one_char_per_hit(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        assert data.hits.strand == '+-'

    def test_pair_members_get_roles_and_element(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        assert data.hits.role == [1, 2]  # left, right
        assert data.hits.pair_ix == [0, 0]
        element = data.elements[0]
        assert element.start == 100
        assert element.end == 5100
        assert element.length == 5001
        assert element.left_uid == 0
        assert element.right_uid == 1
        assert element.inner_distance == 4799

    def test_pair_members_share_a_track_row(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        assert data.hits.row[0] == data.hits.row[1]

    def test_contig_slice_covers_its_hits(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        contig = data.contigs[0]
        assert contig.name == 'chr1'
        assert (contig.hit_lo, contig.hit_hi) == (0, 2)
        assert contig.length == 100_000
        assert contig.length_source == 'source'
        assert contig.n_pairs == 1

    def test_group_colour_is_assigned(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        assert data.groups[0].colour.startswith('#')
        assert data.groups[0].group_id == 'TIR_L__TIR_R'

    def test_sequences_are_embedded(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        assert data.sequences.embedded is True
        assert data.sequences.seq == {'TIR_L__TIR_R:0': 'ACGTACGT'}

    def test_sequences_can_be_disabled(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise(embed_sequences=False)
        assert data.sequences.embedded is False
        assert data.sequences.seq == {}

    def test_sequence_budget_omits_and_warns(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise(max_seq_bytes=1)
        assert data.sequences.seq == {}
        assert data.sequences.omitted == 1
        assert any('too large to embed' in w for w in data.warnings)

    def test_json_round_trips(self, hit_table_factory):
        data = self._accumulate(hit_table_factory).finalise()
        payload = json.loads(data.to_json())
        assert payload['schema_version'] == SCHEMA_VERSION
        assert payload['hits']['n'] == 2
        assert payload['contigs'][0]['name'] == 'chr1'


class TestContigLength:
    def test_inferred_when_no_source(self, hit_table_factory):
        table = hit_table_factory(
            [{'model': 'TIR', 'target': 'chr1', 'hitStart': '100', 'hitEnd': '200'}]
        )
        acc = PairReportAccumulator(hit_table=table, model_lengths={'TIR': 100})
        acc.add_unpaired(make_index([make_hit(0, 'TIR', 'chr1', 100, 200)]))
        data = acc.finalise()
        contig = data.contigs[0]
        assert contig.length_source == 'inferred'
        assert contig.length > 200
        assert any('estimated from the hits' in w for w in data.warnings)

    def test_source_failure_falls_back_to_inferred(self, hit_table_factory):
        def boom(_name):
            raise OSError('index missing')

        table = hit_table_factory(
            [{'model': 'TIR', 'target': 'chr1', 'hitStart': '100', 'hitEnd': '200'}]
        )
        acc = PairReportAccumulator(
            hit_table=table, model_lengths={'TIR': 100}, contig_length=boom
        )
        acc.add_unpaired(make_index([make_hit(0, 'TIR', 'chr1', 100, 200)]))
        assert acc.finalise().contigs[0].length_source == 'inferred'


class TestClipping:
    def test_truncation_at_contig_start_is_flagged_as_clipped(self, hit_table_factory):
        # A hit 5 bp from the contig start, missing 10 model positions on that
        # side: the model cannot be covered because the sequence ran out.
        table = hit_table_factory(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '5',
                    'hitEnd': '95',
                    'strand': '+',
                    'hmmStart': '11',
                    'hmmEnd': '100',
                }
            ]
        )
        acc = PairReportAccumulator(
            hit_table=table,
            model_lengths={'TIR': 100},
            contig_length=lambda name: 1000,
        )
        acc.add_unpaired(make_index([make_hit(0, 'TIR', 'chr1', 5, 95)]))
        data = acc.finalise()
        assert data.hits.trunc_l == [10]
        assert data.hits.clip_l == [1]
        assert data.hits.clip_r == [0]

    def test_truncation_away_from_contig_end_is_not_clipped(self, hit_table_factory):
        table = hit_table_factory(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '500',
                    'hitEnd': '590',
                    'strand': '+',
                    'hmmStart': '11',
                    'hmmEnd': '100',
                }
            ]
        )
        acc = PairReportAccumulator(
            hit_table=table,
            model_lengths={'TIR': 100},
            contig_length=lambda name: 1000,
        )
        acc.add_unpaired(make_index([make_hit(0, 'TIR', 'chr1', 500, 590)]))
        data = acc.finalise()
        assert data.hits.trunc_l == [10]
        assert data.hits.clip_l == [0]


class TestMissingModelCoordinates:
    def test_bed_input_gives_none_coverage_and_no_jag(self, hit_table_factory):
        # import_BED writes the literal string 'NA' for the model coordinates.
        table = hit_table_factory(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'hmmStart': 'NA',
                    'hmmEnd': 'NA',
                }
            ]
        )
        acc = PairReportAccumulator(
            hit_table=table,
            model_lengths={'TIR': 100},
            contig_length=lambda name: 1000,
        )
        acc.add_unpaired(make_index([make_hit(0, 'TIR', 'chr1', 100, 200)]))
        data = acc.finalise()
        assert data.hits.hmm_start == [None]
        assert data.hits.model_cov == [None]
        assert data.hits.trunc_l == [None]
        assert data.hits.clip_l == [0]

    def test_non_numeric_evalue_becomes_none(self, hit_table_factory):
        table = hit_table_factory(
            [{'model': 'TIR', 'target': 'chr1', 'evalue': 'NA', 'score': 'NA'}]
        )
        acc = PairReportAccumulator(hit_table=table)
        acc.add_unpaired(
            make_index([make_hit(0, 'TIR', 'chr1', 100, 200, evalue='NA')])
        )
        data = acc.finalise()
        assert data.hits.evalue == [None]
        assert data.hits.score == [None]


class TestMultipleGroups:
    def test_hit_in_two_pairing_rows_lists_both_groups(self, hit_table_factory):
        # TIR_M appears on the right of one pairing-map row and the left of
        # another, which load_pairing_map explicitly warns about but permits.
        hits = [
            make_hit(0, 'TIR_L', 'chr1', 100, 200),
            make_hit(1, 'TIR_M', 'chr1', 5000, 5100),
            make_hit(2, 'TIR_R', 'chr1', 9000, 9100),
        ]
        table = hit_table_factory(
            [
                {'model': 'TIR_L', 'hitStart': '100', 'hitEnd': '200'},
                {'model': 'TIR_M', 'hitStart': '5000', 'hitEnd': '5100'},
                {'model': 'TIR_R', 'hitStart': '9000', 'hitEnd': '9100'},
            ]
        )
        acc = PairReportAccumulator(hit_table=table, contig_length=lambda name: 100_000)
        index = make_index(hits)
        for left, right, pair in (
            ('TIR_L', 'TIR_M', {0, 1}),
            ('TIR_M', 'TIR_R', {1, 2}),
        ):
            acc.add_group(
                left_feature=left,
                right_feature=right,
                config=PairingConfig(
                    orientation='F,R', left_model=left, right_model=right
                ),
                paired={left: [pair]},
                hit_index=index,
            )
        data = acc.finalise()
        assert [g.group_id for g in data.groups] == ['TIR_L__TIR_M', 'TIR_M__TIR_R']
        # TIR_M is at index 1 after sorting by start.
        assert data.hits.group_ix[1] == [0, 1]
        assert data.hits.group_ix[0] == [0]
        assert data.groups[0].colour != data.groups[1].colour

    def test_symmetric_group_id_does_not_collapse(self, hit_table_factory):
        table = hit_table_factory([{'model': 'TIR'}])
        acc = PairReportAccumulator(hit_table=table)
        acc.add_group(
            left_feature='TIR',
            right_feature='TIR',
            config=PairingConfig(orientation='F,R', single_model='TIR'),
            paired={},
            hit_index=make_index([make_hit(0, 'TIR', 'chr1', 100, 200)]),
        )
        group = acc.finalise().groups[0]
        assert group.group_id == 'TIR__TIR'
        assert group.label == 'TIR'


class TestHitCap:
    def test_cap_keeps_paired_hits_and_warns(self, hit_table_factory):
        hits = [
            make_hit(i, 'TIR', 'chr1', i * 100 + 1, i * 100 + 50) for i in range(10)
        ]
        table = hit_table_factory(
            [
                {
                    'model': 'TIR',
                    'hitStart': str(i * 100 + 1),
                    'hitEnd': str(i * 100 + 50),
                }
                for i in range(10)
            ]
        )
        acc = PairReportAccumulator(
            hit_table=table, contig_length=lambda name: 100_000, max_hits=4
        )
        index = make_index(hits)
        acc.add_group(
            left_feature='TIR',
            right_feature='TIR',
            config=PairingConfig(orientation='F,R', single_model='TIR'),
            paired={'TIR': [{8, 9}]},
            hit_index=index,
        )
        data = acc.finalise()
        assert data.hits.n == 4
        assert {8, 9} <= set(data.hits.uid)
        assert any('omitted to keep the report responsive' in w for w in data.warnings)


class TestJsonEscaping:
    def test_script_close_in_a_model_name_is_escaped(self):
        payload = to_json({'model': 'a</script>b'})
        assert '</script>' not in payload
        assert '\\u003c' in payload
        assert json.loads(payload) == {'model': 'a</script>b'}

    def test_line_separators_are_escaped(self):
        # U+2028/U+2029 are legal in JSON strings but terminate a line in
        # older JavaScript parsers, which would break the embedded payload.
        raw = 'a\u2028b\u2029c'
        payload = to_json({'name': raw})
        assert '\u2028' not in payload
        assert '\u2029' not in payload
        assert json.loads(payload) == {'name': raw}

    def test_report_with_hostile_model_name_stays_parseable(self, hit_table_factory):
        table = hit_table_factory([{'model': 'a</script>b', 'target': 'chr1'}])
        acc = PairReportAccumulator(hit_table=table)
        acc.add_unpaired(make_index([make_hit(0, 'a</script>b', 'chr1', 100, 200)]))
        payload = acc.finalise().to_json()
        assert '</script>' not in payload
        assert json.loads(payload)['models'][0]['name'] == 'a</script>b'
