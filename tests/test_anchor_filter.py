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


# ---------------------------------------------------------------------------
# Reverse-inserted asymmetric elements, and symmetric pairing-map rows
# ---------------------------------------------------------------------------


def _asym_hit(model, strand, hmm_start, hmm_end):
    """One asymmetric hit with the given model-space alignment."""
    return {
        'model': model,
        'target': 'chr1',
        'hitStart': '100',
        'hitEnd': '200',
        'strand': strand,
        'evalue': '1e-10',
        'hmmStart': str(hmm_start),
        'hmmEnd': str(hmm_end),
    }


LENGTHS = {'L': 100, 'R': 100}
PMAP = [('L', 'R')]


class TestAnchorFilterReverseInsertion:
    """
    --max-offset must measure the same model edge in both insertion directions.

    The pairing map gives a terminus ROLE, but the offset is computed from the
    genomic side the outer edge faces. Those diverge for a reverse insertion,
    where the left model's hit is on the opposite strand. Measuring the wrong
    edge silently discarded valid reverse-oriented hits before pairing.
    """

    def test_left_model_at_outer_edge_kept_on_both_strands(self):
        """hmmStart=1 reaches the left model's outer edge either way."""
        for strand in ('+', '-'):
            df = _make_hit_table([_asym_hit('L', strand, 1, 90)])
            result = filter_hits_by_anchor(
                df, LENGTHS, max_offset=5, orientation='F,R', pairing_map=PMAP
            )
            assert len(result) == 1, f'left model on {strand} strand was dropped'

    def test_left_model_short_of_outer_edge_dropped_on_both_strands(self):
        """hmmStart=20 misses it by 19 either way."""
        for strand in ('+', '-'):
            df = _make_hit_table([_asym_hit('L', strand, 20, 100)])
            result = filter_hits_by_anchor(
                df, LENGTHS, max_offset=5, orientation='F,R', pairing_map=PMAP
            )
            assert len(result) == 0, f'left model on {strand} strand was kept'

    def test_right_model_at_outer_edge_kept_on_both_strands(self):
        for strand in ('+', '-'):
            df = _make_hit_table([_asym_hit('R', strand, 1, 90)])
            result = filter_hits_by_anchor(
                df, LENGTHS, max_offset=5, orientation='F,R', pairing_map=PMAP
            )
            assert len(result) == 1, f'right model on {strand} strand was dropped'

    def test_right_model_short_of_outer_edge_dropped_on_both_strands(self):
        for strand in ('+', '-'):
            df = _make_hit_table([_asym_hit('R', strand, 11, 100)])
            result = filter_hits_by_anchor(
                df, LENGTHS, max_offset=5, orientation='F,R', pairing_map=PMAP
            )
            assert len(result) == 0, f'right model on {strand} strand was kept'

    def test_forward_and_reverse_agree_in_every_orientation(self):
        """Insertion direction must never change the verdict."""
        for orientation in ('F,R', 'R,F', 'F,F', 'R,R'):
            for model in ('L', 'R'):
                verdicts = set()
                for strand in ('+', '-'):
                    df = _make_hit_table([_asym_hit(model, strand, 1, 90)])
                    result = filter_hits_by_anchor(
                        df,
                        LENGTHS,
                        max_offset=5,
                        orientation=orientation,
                        pairing_map=PMAP,
                    )
                    verdicts.add(len(result))
                assert len(verdicts) == 1, (
                    f'{orientation} {model}: verdict depends on strand ({verdicts})'
                )


class TestAnchorFilterSymmetricMapRow:
    """
    A pairing-map row naming the same feature twice describes a symmetric
    element, so it has no fixed terminus role and must get the both-ends test.
    Previously the second assignment overwrote the first and such models were
    always treated as right termini.
    """

    def test_symmetric_row_uses_both_ends_test(self):
        """
        Under a same-strand orientation both model ends must be reached.

        A model listed symmetrically has no fixed terminus role, so it must
        fall through to the both-ends rule. Assigning it a role instead lets a
        hit through on the strength of one edge alone.
        """
        keep = _make_hit_table([_asym_hit('LTR', '+', 1, 100)])
        assert (
            len(
                filter_hits_by_anchor(
                    keep,
                    {'LTR': 100},
                    max_offset=5,
                    orientation='F,F',
                    pairing_map=[('LTR', 'LTR')],
                )
            )
            == 1
        )

        # Reaches the model END but misses the START by 19. A model forced to
        # the 'right' role would measure only model_len - hmmEnd = 0 and keep
        # this hit.
        drop_start = _make_hit_table([_asym_hit('LTR', '+', 20, 100)])
        assert (
            len(
                filter_hits_by_anchor(
                    drop_start,
                    {'LTR': 100},
                    max_offset=5,
                    orientation='F,F',
                    pairing_map=[('LTR', 'LTR')],
                )
            )
            == 0
        )

        # Mirror: reaches the START but misses the END by 20.
        drop_end = _make_hit_table([_asym_hit('LTR', '+', 1, 80)])
        assert (
            len(
                filter_hits_by_anchor(
                    drop_end,
                    {'LTR': 100},
                    max_offset=5,
                    orientation='F,F',
                    pairing_map=[('LTR', 'LTR')],
                )
            )
            == 0
        )

    def test_symmetric_row_alongside_asymmetric_rows(self):
        """A symmetric row does not disturb asymmetric rows in the same map."""
        df = _make_hit_table(
            [
                _asym_hit('LTR', '+', 20, 100),  # symmetric: fails both-ends
                _asym_hit('L', '+', 1, 90),  # asymmetric: reaches outer edge
            ]
        )
        result = filter_hits_by_anchor(
            df,
            {'LTR': 100, 'L': 100, 'R': 100},
            max_offset=5,
            orientation='F,F',
            pairing_map=[('LTR', 'LTR'), ('L', 'R')],
        )
        assert list(result['model']) == ['L']
