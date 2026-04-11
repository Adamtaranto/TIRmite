#!/usr/bin/env python3
"""
Tests for asymmetric pairing logic with all orientation combinations.

Validates that:
1. Model assignment from leftBlast/rightBlast uses the correct file-based assignment
   (not alphabetical order from the combined sorted table).
2. Pairing logic correctly handles all orientation combinations:
   F,F (both forward), F,R (canonical), R,F, and R,R.
"""

import os
import tempfile

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
    def left_blast_file(self):
        """BLAST output file for the left query (alphabetically later: query_5p)."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.blast', delete=False) as f:
            # query_5p hits at 100-200 on chr1 (forward strand, sstart < send)
            f.write('query_5p\tchr1\t100.000\t100\t0\t0\t1\t100\t100\t199\t0.0\t185\n')
            fname = f.name
        yield fname
        os.unlink(fname)

    @pytest.fixture
    def right_blast_file(self):
        """BLAST output file for the right query (alphabetically earlier: query_3p)."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.blast', delete=False) as f:
            # query_3p hits at 300-399 on chr1 (forward strand)
            f.write('query_3p\tchr1\t100.000\t100\t0\t0\t1\t100\t300\t399\t0.0\t185\n')
            fname = f.name
        yield fname
        os.unlink(fname)

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
