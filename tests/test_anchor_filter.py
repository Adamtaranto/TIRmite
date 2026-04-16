"""Tests for the --max-offset anchor filter in tirmite pair."""

import pandas as pd

from tirmite.cli.hmm_pair import compute_outer_edge_offset, filter_hits_by_anchor

# ---------------------------------------------------------------------------
# compute_outer_edge_offset tests
# ---------------------------------------------------------------------------


class TestComputeOuterEdgeOffset:
    """Tests for compute_outer_edge_offset."""

    def test_left_plus_at_edge(self):
        """Left terminus, + strand: hmmStart=1 → offset 0."""
        assert compute_outer_edge_offset(1, 100, 100, '+', 'left') == 0

    def test_left_plus_offset_5(self):
        """Left terminus, + strand: hmmStart=6 → offset 5."""
        assert compute_outer_edge_offset(6, 100, 100, '+', 'left') == 5

    def test_left_minus_at_edge(self):
        """Left terminus, - strand: hmmEnd=model_len → offset 0."""
        assert compute_outer_edge_offset(1, 100, 100, '-', 'left') == 0

    def test_left_minus_offset_5(self):
        """Left terminus, - strand: hmmEnd=95 with model_len=100 → offset 5."""
        assert compute_outer_edge_offset(1, 95, 100, '-', 'left') == 5

    def test_right_plus_at_edge(self):
        """Right terminus, + strand: hmmEnd=model_len → offset 0."""
        assert compute_outer_edge_offset(1, 100, 100, '+', 'right') == 0

    def test_right_plus_offset_5(self):
        """Right terminus, + strand: hmmEnd=95, model_len=100 → offset 5."""
        assert compute_outer_edge_offset(1, 95, 100, '+', 'right') == 5

    def test_right_minus_at_edge(self):
        """Right terminus, - strand: hmmStart=1 → offset 0."""
        assert compute_outer_edge_offset(1, 100, 100, '-', 'right') == 0

    def test_right_minus_offset_5(self):
        """Right terminus, - strand: hmmStart=6 → offset 5."""
        assert compute_outer_edge_offset(6, 100, 100, '-', 'right') == 5


# ---------------------------------------------------------------------------
# filter_hits_by_anchor tests
# ---------------------------------------------------------------------------


def _make_hit_table(rows):
    """Helper to build a hit DataFrame."""
    return pd.DataFrame(rows)


class TestFilterHitsByAnchorFR:
    """Tests for F,R orientation (canonical TIR): + strand = left, - strand = right."""

    def test_left_plus_passes_within_offset(self):
        """F,R: left(+) hit with hmmStart=6 (offset=5) passes max_offset=10."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '6',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 1

    def test_left_plus_removed_exceeds_offset(self):
        """F,R: left(+) hit with hmmStart=16 (offset=15) is removed at max_offset=10."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '16',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 0

    def test_right_minus_passes_within_offset(self):
        """F,R: right(-) hit with hmmStart=1 (offset=0) passes max_offset=10."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '-',
                    'evalue': '1e-10',
                    'hmmStart': '1',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 1

    def test_right_minus_removed_exceeds_offset(self):
        """F,R: right(-) hit with hmmStart=20 (offset=19) removed at max_offset=10."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '-',
                    'evalue': '1e-10',
                    'hmmStart': '20',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 0


class TestFilterHitsByAnchorAsymmetric:
    """Tests for asymmetric pairing (model name determines terminus type)."""

    def test_asymmetric_left_model_filtered(self):
        """Asymmetric: left model with large offset removed."""
        df = _make_hit_table(
            [
                {
                    'model': 'LeftModel',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '20',  # offset = 19 > max_offset=10
                    'hmmEnd': '100',
                },
                {
                    'model': 'RightModel',
                    'target': 'chr1',
                    'hitStart': '300',
                    'hitEnd': '400',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '1',  # offset = 0 ≤ max_offset=10
                    'hmmEnd': '100',
                },
            ]
        )
        result = filter_hits_by_anchor(
            df,
            {'LeftModel': 100, 'RightModel': 100},
            max_offset=10,
            orientation='F,F',
            pairing_map=[('LeftModel', 'RightModel')],
        )
        assert len(result) == 1
        assert result.iloc[0]['model'] == 'RightModel'


class TestFilterHitsByAnchorSymmetricSameStrand:
    """Tests for symmetric same-strand (F,F or R,R): must cover both ends."""

    def test_ff_both_ends_within_offset_passes(self):
        """F,F symmetric: hit covering near both ends passes."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '3',  # offset_start = 2
                    'hmmEnd': '98',  # offset_end = 2
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation='F,F'
        )
        assert len(result) == 1

    def test_ff_start_exceeds_offset_removed(self):
        """F,F symmetric: hmmStart too far from position 1 → removed."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '20',  # offset_start = 19 > 5
                    'hmmEnd': '100',  # offset_end = 0
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation='F,F'
        )
        assert len(result) == 0

    def test_ff_end_exceeds_offset_removed(self):
        """F,F symmetric: hmmEnd too far from model_len → removed."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '1',  # offset_start = 0
                    'hmmEnd': '80',  # offset_end = 20 > 5
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation='F,F'
        )
        assert len(result) == 0

    def test_rr_both_ends_within_offset_passes(self):
        """R,R symmetric: hit covering near both ends passes."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '-',
                    'evalue': '1e-10',
                    'hmmStart': '2',  # offset_start = 1
                    'hmmEnd': '99',  # offset_end = 1
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation='R,R'
        )
        assert len(result) == 1

    def test_rr_both_ends_exceed_offset_removed(self):
        """R,R symmetric: neither end within offset → removed."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '-',
                    'evalue': '1e-10',
                    'hmmStart': '20',  # offset_start = 19 > 5
                    'hmmEnd': '80',  # offset_end = 20 > 5
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation='R,R'
        )
        assert len(result) == 0


class TestFilterHitsByAnchorDefault:
    """Tests for default behavior (no --max-offset)."""

    def test_empty_table_returns_empty(self):
        """Empty DataFrame should pass through."""
        df = pd.DataFrame(
            columns=[
                'model',
                'target',
                'hitStart',
                'hitEnd',
                'strand',
                'evalue',
                'hmmStart',
                'hmmEnd',
            ]
        )
        result = filter_hits_by_anchor(df, {'TIR': 100}, max_offset=10)
        assert len(result) == 0

    def test_missing_model_length_keeps_hit(self):
        """Hits for models without known length are kept."""
        df = _make_hit_table(
            [
                {
                    'model': 'Unknown',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '50',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(df, {}, max_offset=5, orientation='F,R')
        assert len(result) == 1


class TestAnchorFilterCLI:
    """Tests for the --max-offset CLI argument."""

    def test_parser_has_max_offset_argument(self):
        """Verify --max-offset is accepted by pair parser."""
        from tirmite.cli.hmm_pair import create_pair_parser

        parser = create_pair_parser()
        args = parser.parse_args(['--max-offset', '10', '--nhmmer-file', 'test.tbl'])
        assert args.max_offset == 10

    def test_parser_max_offset_default_none(self):
        """Verify --max-offset defaults to None."""
        from tirmite.cli.hmm_pair import create_pair_parser

        parser = create_pair_parser()
        args = parser.parse_args(['--nhmmer-file', 'test.tbl'])
        assert args.max_offset is None
