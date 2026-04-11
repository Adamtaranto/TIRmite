#!/usr/bin/env python3
"""
Tests for --pairing-map integration with tirmite pair.

Covers:
- Multi-model BLAST files: independent pairing per map entry
- Missing model warnings and skip
- Models not in map are ignored (not paired)
- --orientation respected for all map entries
- Warning when multiple models in file and no pairing map
- Pairing map file loading with expected formats
"""

import os
import tempfile

import pandas as pd
import pytest

from tirmite.cli.hmm_pair import load_pairing_map
import tirmite.tirmitetools as tirmite

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_hit_table(rows):
    """Build a hit DataFrame from a list of dicts with pairing-relevant fields."""
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


def _run_pairing_with_map(hit_rows, pairing_list, orientation='F,F'):
    """
    Run the full pairing-map loop (mirrors hmm_pair.py logic).

    Parameters
    ----------
    hit_rows : list of dict
    pairing_list : list of (left_model, right_model) tuples
    orientation : str

    Returns
    -------
    all_paired : dict  {model: [set of hit indices]}
    all_paired_hits : set  hit indices that were paired
    unpaired_hits : list  hit indices that remained unpaired
    warned_skipped : list of str  feature names that were warned as missing
    """
    hitTable = _make_hit_table(hit_rows)
    hitsDict, hitIndex = tirmite.table2dict(hitTable)

    all_paired = {}
    all_paired_hits = set()
    original_hitIndex = hitIndex
    warned_skipped = []

    for left_feature, right_feature in pairing_list:
        if left_feature not in hitsDict:
            warned_skipped.append(left_feature)
            continue
        if right_feature not in hitsDict:
            warned_skipped.append(right_feature)
            continue

        if left_feature == right_feature:
            pair_config = tirmite.PairingConfig(
                orientation=orientation, single_model=left_feature
            )
        else:
            pair_config = tirmite.PairingConfig(
                orientation=orientation,
                left_model=left_feature,
                right_model=right_feature,
            )

        pair_hitIndex = tirmite.parseHitsGeneral(
            hitsDict=hitsDict,
            hitIndex=tirmite.table2dict(hitTable)[1],  # fresh index each time
            config=pair_config,
        )

        if pair_config.is_asymmetric:
            _, pair_paired, _ = tirmite.iterateGetPairsAsymmetric(
                pair_hitIndex, pair_config, stableReps=5
            )
        else:
            _, pair_paired, _ = tirmite.iterateGetPairsCustom(
                pair_hitIndex, pair_config, stableReps=5
            )

        for model, pairs in pair_paired.items():
            if model not in all_paired:
                all_paired[model] = []
            all_paired[model].extend(pairs)
            for pair_set in pairs:
                all_paired_hits.update(pair_set)

    # Collect unpaired hits
    unpaired_hits = []
    for model in hitsDict:
        if model in original_hitIndex:
            for hit_id in original_hitIndex[model]:
                if hit_id not in all_paired_hits:
                    unpaired_hits.append(hit_id)

    return all_paired, all_paired_hits, unpaired_hits, warned_skipped


# ---------------------------------------------------------------------------
# Pairing map: multi-model independent pairing
# ---------------------------------------------------------------------------


class TestPairingMapMultiModel:
    """Multi-model blast files with pairing map → independent pairing per entry."""

    def test_two_independent_pairs_from_map(self):
        """Two pairs in map → each produces 1 pair on its own contig."""
        rows = [
            # Pair A: left_A upstream of right_A on chr1
            {
                'model': 'left_A',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_A',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
            # Pair B: left_B upstream of right_B on chr2
            {
                'model': 'left_B',
                'target': 'chr2',
                'hit_start': 500,
                'hit_end': 600,
                'strand': '+',
            },
            {
                'model': 'right_B',
                'target': 'chr2',
                'hit_start': 700,
                'hit_end': 800,
                'strand': '+',
            },
        ]
        pairing_list = [('left_A', 'right_A'), ('left_B', 'right_B')]

        all_paired, all_paired_hits, unpaired, warned = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        assert len(all_paired.get('left_A', [])) == 1, (
            'Expected 1 pair for left_A → right_A'
        )
        assert len(all_paired.get('left_B', [])) == 1, (
            'Expected 1 pair for left_B → right_B'
        )
        assert len(unpaired) == 0, 'All hits should be paired'
        assert warned == []

    def test_pairing_map_multiple_hits_per_model(self):
        """Two left_A hits + two right_A hits → 2 pairs from map."""
        rows = [
            {
                'model': 'left_A',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_A',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
            {
                'model': 'left_A',
                'target': 'chr1',
                'hit_start': 1000,
                'hit_end': 1100,
                'strand': '+',
            },
            {
                'model': 'right_A',
                'target': 'chr1',
                'hit_start': 1200,
                'hit_end': 1300,
                'strand': '+',
            },
        ]
        pairing_list = [('left_A', 'right_A')]

        all_paired, all_paired_hits, unpaired, warned = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        assert len(all_paired.get('left_A', [])) == 2, 'Expected 2 pairs'
        assert len(unpaired) == 0

    def test_pairing_map_fr_orientation(self):
        """Pairing map respects F,R orientation (left on +, right on -)."""
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
                'strand': '-',
            },
        ]
        pairing_list = [('left_q', 'right_q')]

        all_paired, _, unpaired, warned = _run_pairing_with_map(
            rows, pairing_list, orientation='F,R'
        )

        assert len(all_paired.get('left_q', [])) == 1
        assert len(unpaired) == 0

    def test_pairing_map_rr_orientation(self):
        """Pairing map respects R,R orientation (both on - strand)."""
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
        pairing_list = [('left_q', 'right_q')]

        all_paired, _, unpaired, warned = _run_pairing_with_map(
            rows, pairing_list, orientation='R,R'
        )

        assert len(all_paired.get('left_q', [])) == 1
        assert len(unpaired) == 0


# ---------------------------------------------------------------------------
# Pairing map: skipping missing models
# ---------------------------------------------------------------------------


class TestPairingMapMissingModels:
    """Map entries with models absent from hits → warn and skip that entry."""

    def test_missing_left_model_is_skipped(self):
        """Pairing entry with missing left model is skipped; others proceed."""
        rows = [
            # Only Pair A is present
            {
                'model': 'left_A',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_A',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
        ]
        pairing_list = [
            ('left_B', 'right_B'),  # Both missing
            ('left_A', 'right_A'),  # Present
        ]

        all_paired, _, unpaired, warned = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        # left_B is not in hitsDict → should be in warned and skipped
        assert 'left_B' in warned, 'left_B should be warned as missing'
        # Pair A should still be found
        assert len(all_paired.get('left_A', [])) == 1, 'Pair A should still be paired'

    def test_missing_right_model_is_skipped(self):
        """Pairing entry with missing right model is skipped."""
        rows = [
            {
                'model': 'left_A',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_A',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
        ]
        pairing_list = [
            ('left_A', 'right_missing'),  # right_missing absent
            ('left_A', 'right_A'),
        ]

        all_paired, _, _, warned = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        assert 'right_missing' in warned
        assert len(all_paired.get('left_A', [])) == 1

    def test_all_map_models_missing_produces_no_pairs(self):
        """If every entry in map references missing models, zero pairs produced."""
        rows = [
            {
                'model': 'model_X',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
        ]
        pairing_list = [('left_A', 'right_A'), ('left_B', 'right_B')]

        all_paired, _, unpaired, warned = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        total_pairs = sum(len(p) for p in all_paired.values())
        assert total_pairs == 0, 'No pairs when all map models are absent'
        # Each entry warns for the FIRST missing model it encounters (left before right)
        assert len(warned) == len(
            pairing_list
        )  # One warning per entry (left model missing)


# ---------------------------------------------------------------------------
# Pairing map: models not in map are ignored
# ---------------------------------------------------------------------------


class TestPairingMapIgnoresUnlistedModels:
    """Only models listed in the map are paired; extra models are not."""

    def test_extra_model_not_paired(self):
        """Model 'extra' not in map → its hits remain unpaired."""
        rows = [
            {
                'model': 'left_A',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_A',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
            # 'extra' model not mentioned in map
            {
                'model': 'extra',
                'target': 'chr1',
                'hit_start': 500,
                'hit_end': 600,
                'strand': '+',
            },
        ]
        pairing_list = [('left_A', 'right_A')]

        all_paired, all_paired_hits, unpaired, _ = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        assert len(all_paired.get('left_A', [])) == 1
        # The extra model hit should not appear in any pair
        for pair_set in all_paired.get('left_A', []):
            for _hit_idx in pair_set:
                # Get that hit's model from the original rows (by index)
                pass
        # At minimum: total paired hits == 2 (left_A + right_A), extra is unpaired
        assert len(all_paired_hits) == 2

    def test_map_defines_only_subset_of_models(self):
        """With 4 models and map for 1 pair, only that pair's hits are paired."""
        rows = [
            {
                'model': 'left_A',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_A',
                'target': 'chr1',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
            {
                'model': 'left_B',
                'target': 'chr2',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'right_B',
                'target': 'chr2',
                'hit_start': 300,
                'hit_end': 400,
                'strand': '+',
            },
        ]
        # Only pair A is in the map
        pairing_list = [('left_A', 'right_A')]

        all_paired, all_paired_hits, unpaired, _ = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        assert len(all_paired.get('left_A', [])) == 1
        assert len(all_paired_hits) == 2  # left_A + right_A only
        assert len(unpaired) == 2  # left_B and right_B are unpaired


# ---------------------------------------------------------------------------
# Pairing map file format tests (for the pair command context)
# ---------------------------------------------------------------------------


class TestPairingMapFileFormat:
    """Tests for pairing map file loading specific to hmm_pair context."""

    def test_blast_pairing_map_two_column_tsv(self):
        """Two-column TSV pairing map with left/right blast query names loads correctly."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('# Pairing map for blast queries\n')
            f.write('hphi_5p_head_1000\thphi_3p_tail_1000\n')
            f.write('ltr_left_query\tltr_right_query\n')
            fname = f.name
        try:
            pairings = load_pairing_map(fname)
            assert len(pairings) == 2
            assert pairings[0] == ('hphi_5p_head_1000', 'hphi_3p_tail_1000')
            assert pairings[1] == ('ltr_left_query', 'ltr_right_query')
        finally:
            os.unlink(fname)

    def test_pairing_map_comments_and_blank_lines_skipped(self):
        """Comment lines and blank lines are skipped during loading."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('# Header comment\n')
            f.write('\n')
            f.write('left1\tright1\n')
            f.write('\n')
            f.write('# Another comment\n')
            f.write('left2\tright2\n')
            fname = f.name
        try:
            pairings = load_pairing_map(fname)
            assert len(pairings) == 2
            assert pairings[0] == ('left1', 'right1')
            assert pairings[1] == ('left2', 'right2')
        finally:
            os.unlink(fname)


# ---------------------------------------------------------------------------
# Warning when multiple models without pairing map (end-to-end logic check)
# ---------------------------------------------------------------------------


class TestMultipleModelsWarning:
    """Multiple models in hit files without pairing map → warning issued."""

    def test_multiple_models_in_left_file_triggers_warning(self, caplog):
        """Simulates the warning path: check_multiple_models detects >1 model."""
        from tirmite.cli.hmm_pair import check_multiple_models

        left_table = _make_hit_table(
            [
                {
                    'model': 'query_5p',
                    'target': 'chr1',
                    'hit_start': 100,
                    'hit_end': 200,
                    'strand': '+',
                },
                {
                    'model': 'extra_left',
                    'target': 'chr1',
                    'hit_start': 400,
                    'hit_end': 500,
                    'strand': '+',
                },
            ]
        )

        models = check_multiple_models(left_table)
        assert len(models) == 2
        assert 'query_5p' in models
        assert 'extra_left' in models

    def test_single_model_per_file_no_warning_needed(self, caplog):
        """Single model per file → check_multiple_models returns 1 model."""
        from tirmite.cli.hmm_pair import check_multiple_models

        left_table = _make_hit_table(
            [
                {
                    'model': 'query_5p',
                    'target': 'chr1',
                    'hit_start': 100,
                    'hit_end': 200,
                    'strand': '+',
                },
                {
                    'model': 'query_5p',
                    'target': 'chr2',
                    'hit_start': 100,
                    'hit_end': 200,
                    'strand': '+',
                },
            ]
        )

        models = check_multiple_models(left_table)
        assert len(models) == 1
        assert models[0] == 'query_5p'


# ---------------------------------------------------------------------------
# Integration: realistic multi-model scenario (simulates issue report)
# ---------------------------------------------------------------------------


class TestMultiModelPairingIntegration:
    """
    Simulates the reported bug scenario with two BLAST files each containing
    two query names.  The pairing map routes each left/right pair correctly.
    """

    def test_two_queries_two_contigs_with_map(self):
        """
        Left file: hphi_5p + ltr_5p
        Right file: hphi_3p + ltr_3p
        Pairing map: hphi_5p→hphi_3p, ltr_5p→ltr_3p
        F,F orientation: left upstream of right on same contig.
        """
        rows = [
            # hphi pair on JANCMO
            {
                'model': 'hphi_5p',
                'target': 'JANCMO01',
                'hit_start': 57230,
                'hit_end': 58229,
                'strand': '+',
            },
            {
                'model': 'hphi_3p',
                'target': 'JANCMO01',
                'hit_start': 140875,
                'hit_end': 141874,
                'strand': '+',
            },
            # ltr pair on RHLL
            {
                'model': 'ltr_5p',
                'target': 'RHLL01',
                'hit_start': 98017,
                'hit_end': 99016,
                'strand': '+',
            },
            {
                'model': 'ltr_3p',
                'target': 'RHLL01',
                'hit_start': 174581,
                'hit_end': 175580,
                'strand': '+',
            },
        ]
        pairing_list = [('hphi_5p', 'hphi_3p'), ('ltr_5p', 'ltr_3p')]

        all_paired, all_paired_hits, unpaired, warned = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        assert len(all_paired.get('hphi_5p', [])) == 1, (
            'hphi pair should produce 1 pair'
        )
        assert len(all_paired.get('ltr_5p', [])) == 1, 'ltr pair should produce 1 pair'
        assert len(unpaired) == 0, 'All hits should be paired'
        assert warned == []

    def test_alphabetical_name_inversion_with_map(self):
        """
        When the right-query name sorts before the left-query name,
        a pairing map ensures correct left/right assignment.

        query_3p (alphabetically first, right)
        query_5p (alphabetically later, left)
        """
        rows = [
            {
                'model': 'query_5p',
                'target': 'chr1',
                'hit_start': 100,
                'hit_end': 200,
                'strand': '+',
            },
            {
                'model': 'query_3p',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 500,
                'strand': '+',
            },
        ]
        pairing_list = [('query_5p', 'query_3p')]  # explicit left=5p, right=3p

        all_paired, _, unpaired, _ = _run_pairing_with_map(
            rows, pairing_list, orientation='F,F'
        )

        assert len(all_paired.get('query_5p', [])) == 1, (
            'Pairing map should correctly assign query_5p as left and query_3p as right'
        )
        assert len(unpaired) == 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
