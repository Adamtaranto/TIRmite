"""Tests for the --max-offset anchor filter.

Covers both `tirmite pair` and `tirmite search`, which share a single
implementation in tirmite.core.hit_filters.
"""

import pandas as pd
import pytest

from tirmite.cli import ensemble_search as _search_mod
from tirmite.cli import hmm_pair as _pair_mod
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


# ---------------------------------------------------------------------------
# Regression tests for the unification of the two anchor-filter copies
#
# compute_outer_edge_offset and filter_hits_by_anchor were implemented twice,
# independently, in hmm_pair.py and ensemble_search.py. The hmm_pair copy
# received two corrections the ensemble_search copy never did. These tests run
# against BOTH import paths, so the two can never diverge again.
# ---------------------------------------------------------------------------

# Both subcommands must resolve to one implementation.
_ANCHOR_IMPLEMENTATIONS = [
    pytest.param(_pair_mod.filter_hits_by_anchor, id='hmm_pair'),
    pytest.param(_search_mod.filter_hits_by_anchor, id='ensemble_search'),
]


class TestAnchorFilterIsSharedBetweenSubcommands:
    """tirmite pair and tirmite search must use the same implementation."""

    def test_filter_is_the_same_object(self):
        """Both modules re-export the single core implementation."""
        from tirmite.core.hit_filters import filter_hits_by_anchor as canonical

        assert _pair_mod.filter_hits_by_anchor is canonical
        assert _search_mod.filter_hits_by_anchor is canonical

    def test_offset_helper_is_the_same_object(self):
        """The offset helper is shared too."""
        from tirmite.core.hit_filters import compute_outer_edge_offset as canonical

        assert _pair_mod.compute_outer_edge_offset is canonical
        assert _search_mod.compute_outer_edge_offset is canonical


class TestReverseInsertedElements:
    """Reverse-inserted asymmetric elements must survive the anchor filter.

    Under orientation F,R a forward-inserted element puts its LeftA hit on '+'
    and its RightA hit on '-'. When the whole element is inserted in reverse,
    both strands flip: LeftA appears on '-' and RightA on '+'. The model name
    still gives the pairing ROLE, but the genomic SIDE its outer edge faces has
    swapped, so the offset must be measured against the other model end.

    The ensemble_search copy omitted that swap and measured against the hit's
    INNER edge, silently discarding valid reverse-inserted hits. Fixed for
    tirmite pair in 1.5.0; tirmite search never received the fix until the two
    were unified.
    """

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    def test_reverse_inserted_right_model_on_plus_is_kept(self, filter_fn):
        """RightA on '+' is reverse-inserted; its outer edge is at hmmStart."""
        df = _make_hit_table(
            [
                {
                    'model': 'RightA',
                    'target': 'chr1',
                    'hitStart': '1000',
                    'hitEnd': '1100',
                    'strand': '+',
                    'evalue': '1e-10',
                    # Anchored at the model start, which for a reverse-inserted
                    # right terminus IS the outer edge.
                    'hmmStart': '1',
                    'hmmEnd': '60',
                }
            ]
        )
        result = filter_fn(
            df,
            {'LeftA': 100, 'RightA': 100},
            max_offset=5,
            orientation='F,R',
            pairing_map=[('LeftA', 'RightA')],
        )
        assert len(result) == 1, (
            'Reverse-inserted RightA hit anchored at its outer edge was discarded'
        )

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    def test_reverse_inserted_right_model_off_edge_is_removed(self, filter_fn):
        """The swap must not make the filter permissive: interior hits still go."""
        df = _make_hit_table(
            [
                {
                    'model': 'RightA',
                    'target': 'chr1',
                    'hitStart': '1000',
                    'hitEnd': '1100',
                    'strand': '+',
                    'evalue': '1e-10',
                    # 29 positions short of the model start.
                    'hmmStart': '30',
                    'hmmEnd': '60',
                }
            ]
        )
        result = filter_fn(
            df,
            {'LeftA': 100, 'RightA': 100},
            max_offset=5,
            orientation='F,R',
            pairing_map=[('LeftA', 'RightA')],
        )
        assert len(result) == 0

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    def test_reverse_inserted_left_model_on_minus_is_kept(self, filter_fn):
        """LeftA on '-' is reverse-inserted and still anchors at model start.

        The point of the strand swap is that a left-role model's outer edge is
        model position 1 regardless of insertion direction. Without the swap
        this hit would be measured against ``model_len - hmmEnd`` instead.
        """
        df = _make_hit_table(
            [
                {
                    'model': 'LeftA',
                    'target': 'chr1',
                    'hitStart': '1000',
                    'hitEnd': '1100',
                    'strand': '-',
                    'evalue': '1e-10',
                    'hmmStart': '1',
                    'hmmEnd': '60',
                }
            ]
        )
        result = filter_fn(
            df,
            {'LeftA': 100, 'RightA': 100},
            max_offset=5,
            orientation='F,R',
            pairing_map=[('LeftA', 'RightA')],
        )
        assert len(result) == 1

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    def test_left_model_anchor_edge_is_direction_independent(self, filter_fn):
        """A left-role hit off the model start is removed on either strand.

        Pinning both directions together is what stops the swap being
        reintroduced in only one of them.
        """
        rows = [
            {
                'model': 'LeftA',
                'target': 'chr1',
                'hitStart': '1000',
                'hitEnd': '1100',
                'strand': strand,
                'evalue': '1e-10',
                'hmmStart': '40',
                'hmmEnd': '100',
            }
            for strand in ('+', '-')
        ]
        result = filter_fn(
            _make_hit_table(rows),
            {'LeftA': 100, 'RightA': 100},
            max_offset=5,
            orientation='F,R',
            pairing_map=[('LeftA', 'RightA')],
        )
        assert len(result) == 0

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    def test_forward_inserted_still_filtered_normally(self, filter_fn):
        """A forward insertion is unaffected by the reverse-insertion branch."""
        df = _make_hit_table(
            [
                # LeftA on '+' anchored at model start: kept.
                {
                    'model': 'LeftA',
                    'target': 'chr1',
                    'hitStart': '1000',
                    'hitEnd': '1100',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '1',
                    'hmmEnd': '60',
                },
                # LeftA on '+' starting deep inside the model: removed.
                {
                    'model': 'LeftA',
                    'target': 'chr1',
                    'hitStart': '2000',
                    'hitEnd': '2100',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '40',
                    'hmmEnd': '100',
                },
            ]
        )
        result = filter_fn(
            df,
            {'LeftA': 100, 'RightA': 100},
            max_offset=5,
            orientation='F,R',
            pairing_map=[('LeftA', 'RightA')],
        )
        assert len(result) == 1
        assert result.iloc[0]['hitStart'] == '1000'


class TestSymmetricModelInPairingMap:
    """A self-paired model has no fixed terminus role.

    The ensemble_search copy set model_terminus[left]='left' and then
    immediately overwrote it with model_terminus[right]='right' for a row like
    (SymTIR, SymTIR), so every SymTIR hit was treated as a right terminus.
    """

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    def test_self_paired_model_uses_both_ends_under_same_strand(self, filter_fn):
        """Under F,F a self-paired model must satisfy BOTH ends."""
        df = _make_hit_table(
            [
                # Reaches both model ends: kept.
                {
                    'model': 'SymTIR',
                    'target': 'chr1',
                    'hitStart': '1000',
                    'hitEnd': '1100',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '1',
                    'hmmEnd': '100',
                },
                # Reaches the end but not the start: removed. Labelling this
                # model 'right' would have wrongly kept it.
                {
                    'model': 'SymTIR',
                    'target': 'chr1',
                    'hitStart': '2000',
                    'hitEnd': '2100',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '40',
                    'hmmEnd': '100',
                },
            ]
        )
        result = filter_fn(
            df,
            {'SymTIR': 100},
            max_offset=5,
            orientation='F,F',
            pairing_map=[('SymTIR', 'SymTIR')],
        )
        assert len(result) == 1
        assert result.iloc[0]['hitStart'] == '1000'

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    def test_symmetric_model_is_symmetric_everywhere(self, filter_fn):
        """A model listed symmetrically anywhere is symmetric in every row."""
        df = _make_hit_table(
            [
                {
                    'model': 'SymTIR',
                    'target': 'chr1',
                    'hitStart': '2000',
                    'hitEnd': '2100',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '40',
                    'hmmEnd': '100',
                },
            ]
        )
        # SymTIR appears both self-paired and as the left of another pair.
        result = filter_fn(
            df,
            {'SymTIR': 100, 'OtherR': 100},
            max_offset=5,
            orientation='F,F',
            pairing_map=[('SymTIR', 'SymTIR'), ('SymTIR', 'OtherR')],
        )
        assert len(result) == 0


class TestPairingMapShapes:
    """Both pairing-map shapes must behave identically.

    tirmite pair supplies a list of tuples; tirmite search supplies a dict.
    """

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    @pytest.mark.parametrize(
        'pairing_map',
        [
            pytest.param([('LeftA', 'RightA')], id='list-of-tuples'),
            pytest.param({'LeftA': 'RightA'}, id='dict'),
        ],
    )
    def test_both_shapes_give_the_same_result(self, filter_fn, pairing_map):
        """A dict and a tuple list describing the same pairing agree."""
        df = _make_hit_table(
            [
                {
                    'model': 'RightA',
                    'target': 'chr1',
                    'hitStart': '1000',
                    'hitEnd': '1100',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '1',
                    'hmmEnd': '60',
                }
            ]
        )
        result = filter_fn(
            df,
            {'LeftA': 100, 'RightA': 100},
            max_offset=5,
            orientation='F,R',
            pairing_map=pairing_map,
        )
        assert len(result) == 1

    @pytest.mark.parametrize('filter_fn', _ANCHOR_IMPLEMENTATIONS)
    def test_none_pairing_map_falls_back_to_strand(self, filter_fn):
        """Without a pairing map, F,R assigns terminus by strand."""
        df = _make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '1000',
                    'hitEnd': '1100',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '1',
                    'hmmEnd': '60',
                }
            ]
        )
        result = filter_fn(df, {'TIR': 100}, max_offset=5, orientation='F,R')
        assert len(result) == 1
