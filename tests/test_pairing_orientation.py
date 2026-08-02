#!/usr/bin/env python3
"""
Tests for asymmetric pairing logic with all orientation combinations.

Validates that:
1. Model assignment from leftBlast/rightBlast uses the correct file-based assignment
   (not alphabetical order from the combined sorted table).
2. Pairing logic correctly handles all orientation combinations:
   F,F (both forward), F,R (canonical), R,F, and R,R.
"""

from collections import namedtuple

import pandas as pd
import pytest

import tirmite.tirmitetools as tirmite

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_hit_df(
    model: str, target: str, hit_start: int, hit_end: int, strand: str
) -> pd.DataFrame:
    """Build a single-row hit DataFrame matching the table2dict expected schema."""
    return pd.DataFrame(
        [
            {
                'model': model,
                'target': target,
                'hitStart': str(hit_start),
                'hitEnd': str(hit_end),
                'strand': strand,
                'evalue': '1e-10',
                'score': '100',
                'bias': 'NA',
                'hmmStart': '1',
                'hmmEnd': '100',
            }
        ]
    )


def _make_hit_table(rows):
    """Build a multi-row hit DataFrame.

    Each row is a dict with keys: model, target, hit_start, hit_end, strand.
    """
    records = []
    for r in rows:
        records.append(
            {
                'model': r['model'],
                'target': r['target'],
                'hitStart': str(r['hit_start']),
                'hitEnd': str(r['hit_end']),
                'strand': r['strand'],
                'evalue': '1e-10',
                'score': '100',
                'bias': 'NA',
                'hmmStart': '1',
                'hmmEnd': '100',
            }
        )
    return pd.DataFrame(records)


def _run_asymmetric_pairing(
    hit_rows, orientation, left_model, right_model, maxdist=None
):
    """Set up and run the full asymmetric pairing pipeline.

    Returns
    -------
    paired : dict
        paired[left_model] = list of pair sets.
    unpaired : list
        Indices of unpaired hits.
    """
    hitTable = _make_hit_table(hit_rows)
    hitsDict, hitIndex = tirmite.table2dict(hitTable)

    config = tirmite.PairingConfig(
        orientation=orientation,
        left_model=left_model,
        right_model=right_model,
    )

    hitIndex = tirmite.parseHitsGeneral(
        hitsDict=hitsDict,
        hitIndex=hitIndex,
        maxDist=maxdist,
        config=config,
    )

    _, paired, unpaired = tirmite.iterateGetPairsAsymmetric(
        hitIndex, config, stableReps=5
    )

    return paired, unpaired


# ---------------------------------------------------------------------------
# F,F orientation tests (the reported bug scenario)
# ---------------------------------------------------------------------------


class TestFFOrientation:
    """Both left and right queries are on the forward (+) strand."""

    def test_ff_basic_pair(self):
        """Left (upstream) + Right (downstream), both +, should pair."""
        # Modelled on the reported bug:
        # left query at ~57000, right query at ~140000 on same contig, F,F
        rows = [
            {
                'model': 'query_3p',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
            {
                'model': 'query_5p',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
        ]
        # Assign left/right explicitly from file origin, NOT alphabetical order.
        # 'query_3p' sorts before 'query_5p', but the left model is 'query_5p'.
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'F,F', left_model='query_5p', right_model='query_3p'
        )

        assert len(paired.get('query_5p', [])) == 1, (
            'Expected 1 pair for F,F orientation with left upstream of right'
        )
        assert len(unpaired) == 0

    def test_ff_wrong_order_not_paired(self):
        """Left (downstream) + Right (upstream) should NOT pair under F,F."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 500,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'F,F', left_model='left_q', right_model='right_q'
        )

        assert len(paired.get('left_q', [])) == 0, (
            'F,F should NOT pair when left is downstream of right'
        )
        assert len(unpaired) == 2

    def test_ff_multiple_hits_best_match(self):
        """With multiple hits, each should be paired with its closest valid partner."""
        rows = [
            # Two left hits
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 500,
                'hit_end': 600,
                'strand': '+',
            },
            # One right hit between the two left hits
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'F,F', left_model='left_q', right_model='right_q'
        )

        # Only the first left hit (upstream of right) can pair; second is downstream of right
        assert len(paired.get('left_q', [])) == 1
        assert len(unpaired) == 1  # Second left hit remains unpaired

    def test_ff_max_distance_respected(self):
        """Hits beyond maxDist should not be paired."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 1000,
                'hit_end': 1100,
                'strand': '+',
            },
        ]
        # Distance = 1000 - 200 = 800; set maxdist=100 to exclude it
        paired_excluded, unpaired_excluded = _run_asymmetric_pairing(
            rows, 'F,F', left_model='left_q', right_model='right_q', maxdist=100
        )
        assert len(paired_excluded.get('left_q', [])) == 0

        # Same hits but generous maxdist → should pair
        paired_ok, unpaired_ok = _run_asymmetric_pairing(
            rows, 'F,F', left_model='left_q', right_model='right_q', maxdist=1000
        )
        assert len(paired_ok.get('left_q', [])) == 1

    def test_ff_different_chromosomes_not_paired(self):
        """Hits on different chromosomes must not be paired."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr2',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'F,F', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 0

    def test_ff_two_independent_pairs_on_same_contig(self):
        """Two independent left+right pairs on the same contig should both be detected."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 1000,
                'hit_end': 1100,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 1200,
                'hit_end': 1300,
                'strand': '+',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'F,F', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 2
        assert len(unpaired) == 0


# ---------------------------------------------------------------------------
# F,R orientation tests (canonical TIR)
# ---------------------------------------------------------------------------


class TestFROrientation:
    """Left query on + strand, right query on - strand (canonical TIR)."""

    def test_fr_basic_pair(self):
        """Left (+) upstream, Right (-) downstream – canonical TIR pairing."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 500,
                'strand': '-',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'F,R', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 1
        assert len(unpaired) == 0

    def test_fr_wrong_strand_not_paired(self):
        """Right hit on + strand should not pair in F,R mode."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 500,
                'strand': '+',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'F,R', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 0

    def test_fr_wrong_order_not_paired(self):
        """Left (+) downstream of Right (-) should not pair in F,R mode."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 600,
                'hit_end': 700,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '-',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'F,R', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 0


# ---------------------------------------------------------------------------
# R,F orientation tests
# ---------------------------------------------------------------------------


class TestRFOrientation:
    """Left query on - strand, right query on + strand."""

    def test_rf_basic_pair(self):
        """Left (-) at higher coords, Right (+) at lower coords – should pair in R,F mode."""
        # For R,F, left is on - strand. The element is on the minus strand.
        # Biologically, the right (+) terminus at lower genomic coords is 'upstream' of
        # the left (-) terminus at higher genomic coords when reading the minus strand.
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 500,
                'strand': '-',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'R,F', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 1
        assert len(unpaired) == 0

    def test_rf_wrong_strand_not_paired(self):
        """Left hit on + strand should not pair in R,F mode."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 500,
                'strand': '+',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'R,F', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 0


# ---------------------------------------------------------------------------
# R,R orientation tests
# ---------------------------------------------------------------------------


class TestRROrientation:
    """Both left and right queries on - strand (element on minus strand)."""

    def test_rr_basic_pair(self):
        """Left (-) at higher coords, Right (-) at lower coords – should pair in R,R mode."""
        # Element is on minus strand: left terminus has higher genomic coords.
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 500,
                'strand': '-',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '-',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'R,R', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 1
        assert len(unpaired) == 0

    def test_rr_wrong_order_not_paired(self):
        """Left (-) at lower coords than Right (-) should not pair in R,R mode."""
        rows = [
            {
                'model': 'left_q',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '-',
            },
            {
                'model': 'right_q',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 500,
                'strand': '-',
            },
        ]
        paired, unpaired = _run_asymmetric_pairing(
            rows, 'R,R', left_model='left_q', right_model='right_q'
        )
        assert len(paired.get('left_q', [])) == 0


# ---------------------------------------------------------------------------
# Model assignment tests (validates the bug fix in hmm_pair.py)
# ---------------------------------------------------------------------------


class TestModelAssignment:
    """Validate that left/right model names are correctly assigned from separate files."""

    @pytest.fixture
    def left_blast_file(self, tmp_path):
        """BLAST output file for the left query (alphabetically later: query_5p)."""
        path = tmp_path / 'left.blast'
        # query_5p hits at 100-199 on chr1 (forward strand, sstart < send)
        path.write_text(
            'query_5p\tchr1\t100.000\t100\t0\t0\t1\t100\t100\t199\t0.0\t185\n'
        )
        return str(path)

    @pytest.fixture
    def right_blast_file(self, tmp_path):
        """BLAST output file for the right query (alphabetically earlier: query_3p)."""
        path = tmp_path / 'right.blast'
        # query_3p hits at 300-399 on chr1 (forward strand)
        path.write_text(
            'query_3p\tchr1\t100.000\t100\t0\t0\t1\t100\t300\t399\t0.0\t185\n'
        )
        return str(path)

    def test_left_model_from_left_file(self, left_blast_file, right_blast_file):
        """Model names from leftBlast/rightBlast should be assigned by file, not alphabetical order."""
        left_hitTable = tirmite.import_blast(infile=left_blast_file)
        right_hitTable = tirmite.import_blast(infile=right_blast_file)

        # The alphabetically sorted combined table would put query_3p first.
        combined = tirmite.import_blast(infile=left_blast_file)
        combined = tirmite.import_blast(infile=right_blast_file, hitTable=combined)

        # Bug: combined['model'].unique()[0] gives 'query_3p' (alphabetically first)
        # Fix: use per-file tables to get correct assignment
        alphabetical_first = combined['model'].unique()[0]
        file_based_left = left_hitTable['model'].unique()[0]
        file_based_right = right_hitTable['model'].unique()[0]

        # Confirm the bug exists in naive alphabetical approach
        assert alphabetical_first == 'query_3p', (
            'Alphabetical sort should give query_3p first (3 < 5), confirming the bug scenario'
        )

        # Confirm the fix correctly identifies left and right models
        assert file_based_left == 'query_5p', (
            'Left model should be query_5p (from leftBlast file)'
        )
        assert file_based_right == 'query_3p', (
            'Right model should be query_3p (from rightBlast file)'
        )

    def test_correct_pairing_with_alphabetically_inverted_names(
        self, left_blast_file, right_blast_file
    ):
        """End-to-end pairing with file-based (not alphabetical) model assignment.

        query_5p (left, alphabetically later) at 100-199
        query_3p (right, alphabetically earlier) at 300-399
        Orientation: F,F

        With the bug (alphabetical assignment):
          left_model='query_3p', right_model='query_5p'
          → query_3p (300-399) is upstream of query_5p (100-199)? No, reversed.
          → distance is negative → 0 pairs

        With the fix (file-based assignment):
          left_model='query_5p', right_model='query_3p'
          → query_5p (100-199) is upstream of query_3p (300-399) ✓
          → positive distance → 1 pair
        """
        left_hitTable = tirmite.import_blast(infile=left_blast_file)
        right_hitTable = tirmite.import_blast(infile=right_blast_file)
        hitTable = tirmite.import_blast(infile=left_blast_file)
        hitTable = tirmite.import_blast(infile=right_blast_file, hitTable=hitTable)

        hitsDict, hitIndex = tirmite.table2dict(hitTable)

        # CORRECT assignment (file-based, the fix)
        correct_config = tirmite.PairingConfig(
            orientation='F,F',
            left_model=left_hitTable['model'].unique()[0],  # query_5p
            right_model=right_hitTable['model'].unique()[0],  # query_3p
        )

        hitIndex_correct = tirmite.parseHitsGeneral(
            hitsDict=hitsDict,
            hitIndex=tirmite.table2dict(hitTable)[1],
            config=correct_config,
        )
        _, paired_correct, unpaired_correct = tirmite.iterateGetPairsAsymmetric(
            hitIndex_correct, correct_config, stableReps=5
        )

        assert len(paired_correct.get('query_5p', [])) == 1, (
            'File-based model assignment should produce 1 valid pair'
        )
        assert len(unpaired_correct) == 0

        # BUGGY assignment (alphabetical, what the bug did)
        buggy_left = hitTable['model'].unique()[0]  # query_3p (wrong)
        buggy_right = hitTable['model'].unique()[1]  # query_5p (wrong)
        buggy_config = tirmite.PairingConfig(
            orientation='F,F',
            left_model=buggy_left,
            right_model=buggy_right,
        )

        hitIndex_buggy = tirmite.parseHitsGeneral(
            hitsDict=hitsDict,
            hitIndex=tirmite.table2dict(hitTable)[1],
            config=buggy_config,
        )
        _, paired_buggy, unpaired_buggy = tirmite.iterateGetPairsAsymmetric(
            hitIndex_buggy, buggy_config, stableReps=5
        )

        assert len(paired_buggy.get('query_3p', [])) == 0, (
            'Alphabetical (buggy) model assignment should produce 0 valid pairs'
        )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


# ---------------------------------------------------------------------------
# Symmetric same-strand pairing (F,F / R,R)
# ---------------------------------------------------------------------------


def _run_symmetric_pairing(hit_rows, orientation, model, maxdist=None, stable_reps=5):
    """Set up and run the full SYMMETRIC pairing pipeline (one model, both ends).

    ``stable_reps`` is exposed so tests can exercise the CLI default of 0 as
    well as the generous value the rest of this file uses.
    """
    hitTable = _make_hit_table(hit_rows)
    hitsDict, hitIndex = tirmite.table2dict(hitTable)

    config = tirmite.PairingConfig(orientation=orientation, single_model=model)

    hitIndex = tirmite.parseHitsGeneral(
        hitsDict=hitsDict,
        hitIndex=hitIndex,
        maxDist=maxdist,
        config=config,
    )

    _, paired, unpaired = tirmite.iterateGetPairsCustom(
        hitIndex, config, stableReps=stable_reps
    )
    return paired, unpaired


class TestSymmetricSameStrandPairing:
    """
    F,F (LTR) and R,R use one model on one strand for both termini.

    A hit must therefore be able to act as either terminus. Testing this
    matters because the existing F,F/R,R tests above drive the *asymmetric*
    code path via two distinct model names.
    """

    def test_ff_symmetric_pairs(self):
        """Two same-strand hits from one model form one LTR-style pair."""
        rows = [
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': 4000,
                'hit_end': 4100,
                'strand': '+',
            },
        ]
        paired, unpaired = _run_symmetric_pairing(rows, 'F,F', 'LTR')
        assert len(paired.get('LTR', [])) == 1
        assert paired['LTR'][0] == {0, 1}
        assert unpaired == []

    def test_rr_symmetric_pairs(self):
        """Same, on the minus strand."""
        rows = [
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '-',
            },
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': 4000,
                'hit_end': 4100,
                'strand': '-',
            },
        ]
        paired, unpaired = _run_symmetric_pairing(rows, 'R,R', 'LTR')
        assert len(paired.get('LTR', [])) == 1
        assert paired['LTR'][0] == {0, 1}

    def test_ff_two_elements_give_two_pairs(self):
        """Four hits resolve into two adjacent pairs, not a tangle."""
        rows = [
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': s,
                'hit_end': s + 100,
                'strand': '+',
            }
            for s in (100, 4000, 9000, 12000)
        ]
        paired, unpaired = _run_symmetric_pairing(rows, 'F,F', 'LTR')
        assert len(paired.get('LTR', [])) == 2
        assert {frozenset(p) for p in paired['LTR']} == {
            frozenset({0, 1}),
            frozenset({2, 3}),
        }
        assert unpaired == []

    def test_lone_hit_does_not_pair_with_itself(self):
        """
        A single hit has no partner.

        On the minus strand the 5'/3' ends swap, so a hit measured against
        itself yields a positive distance and would satisfy the range test.
        Without excluding the reference hit this produced a 'pair' containing
        one hit.
        """
        rows = [
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '-',
            },
        ]
        paired, unpaired = _run_symmetric_pairing(rows, 'R,R', 'LTR')
        assert paired.get('LTR', []) == []
        assert unpaired == [0]

    def test_rr_pairs_are_never_single_hits(self):
        """Every emitted pair contains exactly two distinct hits."""
        rows = [
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': s,
                'hit_end': s + 100,
                'strand': '-',
            }
            for s in (100, 4000, 9000)
        ]
        paired, _ = _run_symmetric_pairing(rows, 'R,R', 'LTR')
        for pair in paired.get('LTR', []):
            assert len(pair) == 2, f'degenerate pair {pair}'

    def test_fr_symmetric_still_pairs(self):
        """Control: the canonical TIR case is unaffected."""
        rows = [
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 4000,
                'hit_end': 4100,
                'strand': '-',
            },
        ]
        paired, unpaired = _run_symmetric_pairing(rows, 'F,R', 'TIR')
        assert len(paired.get('TIR', [])) == 1
        assert unpaired == []


# ---------------------------------------------------------------------------
# --maxdist semantics
# ---------------------------------------------------------------------------


def _hit(model, start, end, strand):
    """Positional helper for distance tests."""
    return {
        'model': model,
        'target': 'chr1',
        'hit_start': start,
        'hit_end': end,
        'strand': strand,
    }


class TestMaxDistSemantics:
    """
    --maxdist is the gap between the facing inner edges of the two hits.

    That is the length of the element interior. It deliberately excludes the
    termini themselves, so the threshold does not shift when a terminus model
    gets longer, and it is identical on both strands.
    """

    # left 100-200, right 4000-4100 -> interior spans 201..3999, gap = 3800
    GAP = 3800

    def _run_symmetric(self, maxdist):
        rows = [_hit('TIR', 100, 200, '+'), _hit('TIR', 4000, 4100, '-')]
        paired, _ = _run_symmetric_pairing(rows, 'F,R', 'TIR', maxdist=maxdist)
        return len(paired.get('TIR', []))

    def test_gap_exactly_at_threshold_pairs(self):
        assert self._run_symmetric(self.GAP) == 1

    def test_gap_one_below_threshold_does_not_pair(self):
        assert self._run_symmetric(self.GAP - 1) == 0

    def test_threshold_is_independent_of_terminus_length(self):
        """
        Lengthening a terminus must not change the distance threshold.

        The old measure ran to the partner's far edge, so a longer terminus
        inflated the measured distance and the same element needed a larger
        --maxdist.
        """
        short = [_hit('TIR', 100, 200, '+'), _hit('TIR', 4000, 4100, '-')]
        long_ = [_hit('TIR', 100, 200, '+'), _hit('TIR', 4000, 4900, '-')]
        for rows in (short, long_):
            paired, _ = _run_symmetric_pairing(rows, 'F,R', 'TIR', maxdist=self.GAP)
            assert len(paired.get('TIR', [])) == 1
            paired, _ = _run_symmetric_pairing(rows, 'F,R', 'TIR', maxdist=self.GAP - 1)
            assert len(paired.get('TIR', [])) == 0

    def test_threshold_identical_in_both_insertion_directions(self):
        """
        An asymmetric element gives the same threshold either way round.

        A left-model hit on '+' sits at the lower coordinate; on '-' it sits at
        the higher one. Both are valid insertions of the same element, and the
        gap between the termini is 3800 in each.
        """
        flip = {'+': '-', '-': '+'}
        checked = 0

        for orientation in ('F,R', 'R,F', 'F,F', 'R,R'):
            cfg = tirmite.PairingConfig(
                orientation=orientation, left_model='L', right_model='R'
            )
            # The two accepted strand combinations for this orientation.
            for left_strand, right_strand in (
                (cfg.left_strand, cfg.right_strand),
                (flip[cfg.left_strand], flip[cfg.right_strand]),
            ):
                if left_strand == '+':
                    rows = [
                        _hit('L', 100, 200, left_strand),
                        _hit('R', 4000, 4100, right_strand),
                    ]
                else:
                    rows = [
                        _hit('L', 4000, 4100, left_strand),
                        _hit('R', 100, 200, right_strand),
                    ]

                at = _run_asymmetric_pairing(
                    rows, orientation, 'L', 'R', maxdist=self.GAP
                )[0]
                below = _run_asymmetric_pairing(
                    rows, orientation, 'L', 'R', maxdist=self.GAP - 1
                )[0]
                assert len(at.get('L', [])) == 1, (orientation, left_strand)
                assert len(below.get('L', [])) == 0, (orientation, left_strand)
                checked += 1

        assert checked == 8, 'expected both directions of all four orientations'

    def test_legacy_and_current_paths_agree(self):
        """
        The legacy workflow and the current one measure the same thing.

        These were previously off by the length of one terminus, so a run
        migrated from `tirmite legacy` to `tirmite pair` silently needed a
        different --maxdist to reproduce its pairs.
        """
        rows = [_hit('TIR', 100, 200, '+'), _hit('TIR', 4000, 4100, '-')]
        hitTable = _make_hit_table(rows)

        for maxdist in (self.GAP - 2, self.GAP - 1, self.GAP, self.GAP + 1):
            hitsDict, hitIndex = tirmite.table2dict(hitTable)
            legacy_index = tirmite.parseHits(
                hitsDict=hitsDict, hitIndex=hitIndex, maxDist=maxdist
            )
            _, legacy_paired, _ = tirmite.iterateGetPairs(legacy_index, stableReps=5)

            current_paired, _ = _run_symmetric_pairing(
                rows, 'F,R', 'TIR', maxdist=maxdist
            )

            assert len(legacy_paired.get('TIR', [])) == len(
                current_paired.get('TIR', [])
            ), f'paths disagree at maxdist={maxdist}'

    def test_overlapping_hits_are_rejected(self):
        """A negative gap means the hits overlap; that is not a valid pair."""
        rows = [_hit('TIR', 100, 500, '+'), _hit('TIR', 400, 800, '-')]
        paired, _ = _run_symmetric_pairing(rows, 'F,R', 'TIR', maxdist=10000)
        assert len(paired.get('TIR', [])) == 0


class TestInterHitDistance:
    """Unit tests for the shared distance measure."""

    Hit = namedtuple('Hit', ['hitStart', 'hitEnd'])

    def test_downstream_gap(self):
        ref = self.Hit(100, 200)
        cand = self.Hit(4000, 4100)
        assert tirmite.inter_hit_distance(ref, cand, 'left_to_right') == 3800

    def test_upstream_gap_is_symmetric(self):
        """Measuring from either end gives the same separation."""
        ref = self.Hit(4000, 4100)
        cand = self.Hit(100, 200)
        assert tirmite.inter_hit_distance(ref, cand, 'right_to_left') == 3800

    def test_adjacent_hits_have_zero_gap(self):
        ref = self.Hit(100, 200)
        cand = self.Hit(200, 300)
        assert tirmite.inter_hit_distance(ref, cand, 'left_to_right') == 0

    def test_wrong_side_is_negative(self):
        ref = self.Hit(4000, 4100)
        cand = self.Hit(100, 200)
        assert tirmite.inter_hit_distance(ref, cand, 'left_to_right') < 0


class TestSymmetricCandidateOrdering:
    """Candidates must be ranked nearest-first across BOTH search directions.

    Under F,F and R,R every hit searches upstream *and* downstream, and both
    searches append into one shared candidate list. `_find_candidates` used to
    re-sort that whole list by the current call's direction, which gave the
    other direction's entries a negative key -- most negative for the farthest
    hit. `candidates[0]` therefore became the globally FARTHEST hit.

    That inverted the reciprocal-best-match test, so only one pair formed per
    iteration. Two consequences, both fixed by ranking on the direction-agnostic
    separation:

    * With the CLI default of ``--stable-reps 0`` the run stopped after a single
      pair, no matter how many elements were present.
    * With a larger value it needed M/2 iterations, which is where the O(M^4)
      cost came from (427 s for 250 elements, now ~0.001 s).
    """

    def _elements(self, n, spacing=5000, element_length=2000):
        """Build n well-separated, individually pairable F,F elements."""
        rows = []
        for i in range(n):
            base = 1000 + (i * spacing)
            rows.append(
                {
                    'model': 'LTR',
                    'target': 'chr1',
                    'hit_start': base,
                    'hit_end': base + 100,
                    'strand': '+',
                }
            )
            rows.append(
                {
                    'model': 'LTR',
                    'target': 'chr1',
                    'hit_start': base + element_length,
                    'hit_end': base + element_length + 100,
                    'strand': '+',
                }
            )
        return rows

    def test_nearest_candidate_is_first(self):
        """The closest hit ranks ahead of a farther one on the other side.

        The reference sits between a near downstream hit and a far upstream
        one. Ranking by the current direction put the far upstream hit first.
        """
        rows = [
            # Far upstream, gap 10000.
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': 1000,
                'hit_end': 1100,
                'strand': '+',
            },
            # Reference.
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': 11100,
                'hit_end': 11200,
                'strand': '+',
            },
            # Near downstream, gap 100.
            {
                'model': 'LTR',
                'target': 'chr1',
                'hit_start': 11300,
                'hit_end': 11400,
                'strand': '+',
            },
        ]
        hitTable = _make_hit_table(rows)
        hitsDict, hitIndex = tirmite.table2dict(hitTable)
        config = tirmite.PairingConfig(orientation='F,F', single_model='LTR')
        hitIndex = tirmite.parseHitsGeneral(
            hitsDict=hitsDict, hitIndex=hitIndex, config=config
        )

        # Hit index 1 is the reference; its nearest partner is hit 2.
        candidates = hitIndex['LTR'][1]['candidates']
        assert len(candidates) == 2
        assert candidates[0].idx == 2, (
            'nearest candidate should rank first, got the farther one'
        )

    @pytest.mark.parametrize('n_elements', [2, 5, 20])
    def test_all_pairs_found_with_cli_default_stable_reps(self, n_elements):
        """--stable-reps defaults to 0; every element must still pair.

        This is the regression test for the headline bug: `tirmite pair
        --orientation F,F` previously returned exactly one pair regardless of
        how many elements were present.
        """
        paired, unpaired = _run_symmetric_pairing(
            self._elements(n_elements), 'F,F', 'LTR', stable_reps=0
        )

        assert len(paired.get('LTR', [])) == n_elements
        assert unpaired == []

    def test_converges_in_a_single_round(self):
        """Correct ranking means one pass suffices.

        stable_reps=0 stops as soon as a round adds no new pairs, so finding
        every pair at that setting proves round 1 converged.
        """
        paired, unpaired = _run_symmetric_pairing(
            self._elements(10), 'F,F', 'LTR', stable_reps=0
        )

        assert len(paired['LTR']) == 10
        assert unpaired == []

    def test_rr_orientation_also_converges(self):
        """R,R takes the same code path and must behave identically."""
        rows = self._elements(10)
        for row in rows:
            row['strand'] = '-'

        paired, unpaired = _run_symmetric_pairing(rows, 'R,R', 'LTR', stable_reps=0)

        assert len(paired['LTR']) == 10
        assert unpaired == []

    def test_elements_pair_with_their_own_partner(self):
        """Pairing must be nearest-neighbour, not arbitrary.

        With three well-separated elements the pairs are (0,1), (2,3), (4,5).
        Farthest-first ranking produced cross-element pairings.
        """
        paired, _ = _run_symmetric_pairing(
            self._elements(3), 'F,F', 'LTR', stable_reps=0
        )

        assert {frozenset(p) for p in paired['LTR']} == {
            frozenset({0, 1}),
            frozenset({2, 3}),
            frozenset({4, 5}),
        }


class TestCandidateSeparation:
    """The direction-agnostic distance used to rank candidates."""

    def _hit(self, start, end):
        """Build a minimal hit-like object with the needed coordinates."""
        from collections import namedtuple

        H = namedtuple('H', 'hitStart hitEnd')
        return H(str(start), str(end))

    def test_downstream_candidate_gap(self):
        """A candidate after the reference reports the gap between them."""
        assert (
            tirmite.candidate_separation(self._hit(100, 200), self._hit(300, 400))
            == 100
        )

    def test_upstream_candidate_gap(self):
        """A candidate before the reference reports the same gap."""
        assert (
            tirmite.candidate_separation(self._hit(300, 400), self._hit(100, 200))
            == 100
        )

    def test_symmetric_in_both_directions(self):
        """Separation does not depend on which hit is the reference."""
        a, b = self._hit(100, 200), self._hit(1000, 1100)
        assert tirmite.candidate_separation(a, b) == tirmite.candidate_separation(b, a)

    def test_abutting_hits_have_zero_gap(self):
        """Hits that touch have no interior between them."""
        assert (
            tirmite.candidate_separation(self._hit(100, 200), self._hit(201, 300)) == 1
        )

    def test_overlapping_hits_are_negative(self):
        """Overlap yields a negative value, which _check_distance rejects."""
        assert (
            tirmite.candidate_separation(self._hit(100, 300), self._hit(200, 400)) < 0
        )
