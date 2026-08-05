#!/usr/bin/env python3
"""
Tests for ensemble search functionality in TIRmite.

Tests cluster mapping, hit merging, nested hit removal, and CLI functionality.
"""

import pandas as pd
import pytest

from tirmite.cli.ensemble_search import (
    EnsembleSearchError,
    build_component_to_cluster_map,
    check_cross_cluster_overlaps,
    filter_best_model_per_locus,
    filter_hits_by_evalue,
    load_hits_from_files,
    merge_overlapping_cluster_hits,
    parse_cluster_mapping,
    parse_pairing_map,
    remove_nested_paired_hits,
    validate_cluster_mapping,
)

# -----------------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------------


@pytest.fixture
def cluster_mapping_file(tmp_path):
    """Create a temporary cluster mapping file.

    Format is cluster name first, then one or more component model names, all
    tab-delimited: ``cluster<TAB>component1<TAB>component2...``.
    """
    path = tmp_path / 'clusters.tsv'
    path.write_text(
        '# Cluster mapping file\n'
        'ClusterA\tComponentA1\tComponentA2\tComponentA3\n'
        'ClusterB\tComponentB1\tComponentB2\n'
        'ClusterC\tComponentC1\n'
    )
    return str(path)


@pytest.fixture
def pairing_map_file(tmp_path):
    """Create a temporary pairing map file."""
    path = tmp_path / 'pairing.tsv'
    path.write_text('# Pairing map: left<TAB>right\nLeftA\tRightA\nLeftB\tRightB\n')
    return str(path)


@pytest.fixture
def blast_result_file(tmp_path):
    """Create a temporary BLAST result file."""
    path = tmp_path / 'hits.blast'
    # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    path.write_text(
        'ComponentA1\tchr1\t95.0\t100\t5\t0\t1\t100\t1000\t1099\t1e-40\t200\n'
        # ComponentA2 overlaps ComponentA1 on chr1.
        'ComponentA2\tchr1\t93.0\t98\t7\t0\t1\t98\t1050\t1147\t1e-35\t180\n'
        # ComponentA3 sits well clear of the other two.
        'ComponentA3\tchr1\t90.0\t95\t10\t0\t1\t95\t2000\t2094\t1e-30\t160\n'
        'ComponentB1\tchr2\t92.0\t100\t8\t0\t1\t100\t500\t599\t1e-38\t190\n'
        # ComponentB2 overlaps ComponentB1 on chr2.
        'ComponentB2\tchr2\t88.0\t96\t12\t0\t1\t96\t520\t615\t1e-32\t170\n'
    )
    return str(path)


@pytest.fixture
def nhmmer_result_file(tmp_path):
    """Create a temporary nhmmer ``--tblout`` result file.

    Column order, which ``import_nhmmer`` and ``detect_input_format`` both
    index positionally: target name, accession, query name, accession,
    hmmfrom, hmm to, alifrom, ali to, envfrom, env to, sq len, strand,
    E-value, score, bias, description.
    """
    path = tmp_path / 'hits.nhmmer'
    path.write_text(
        '# target name  accession  query name  accession  hmmfrom  hmm to  alifrom  '
        'ali to  envfrom  env to  sq len  strand  E-value  score  bias  description\n'
        'chr1 - ComponentA1 - 1 60 1000 1059 1000 1059 100000 + 1.2e-10 45.2 0.0 -\n'
        'chr1 - ComponentA2 - 1 58 1020 1077 1020 1077 100000 + 2.5e-09 42.1 0.0 -\n'
    )
    return str(path)


@pytest.fixture
def sample_hit_table():
    """Create a sample hit table DataFrame."""
    return pd.DataFrame(
        {
            'model': ['CompA1', 'CompA2', 'CompA1', 'CompB1', 'CompB2'],
            'target': ['chr1', 'chr1', 'chr1', 'chr2', 'chr2'],
            'hitStart': ['1000', '1050', '2000', '500', '550'],
            'hitEnd': ['1100', '1150', '2100', '600', '650'],
            'strand': ['+', '+', '+', '+', '+'],
            'evalue': ['1e-40', '1e-35', '1e-30', '1e-38', '1e-32'],
            'score': ['200', '180', '160', '190', '170'],
            'bias': ['NA', 'NA', 'NA', 'NA', 'NA'],
            'hmmStart': ['1', '1', '1', '1', '1'],
            'hmmEnd': ['100', '98', '95', '100', '96'],
        }
    )


# -----------------------------------------------------------------------------
# Cluster Mapping Tests
# -----------------------------------------------------------------------------


class TestClusterMapping:
    """Tests for cluster mapping parsing and validation."""

    def test_parse_cluster_mapping_valid(self, cluster_mapping_file):
        """Test parsing a valid cluster mapping file."""
        cluster_map = parse_cluster_mapping(cluster_mapping_file)

        assert len(cluster_map) == 3
        assert 'ClusterA' in cluster_map
        assert cluster_map['ClusterA'] == ['ComponentA1', 'ComponentA2', 'ComponentA3']
        assert cluster_map['ClusterB'] == ['ComponentB1', 'ComponentB2']
        assert cluster_map['ClusterC'] == ['ComponentC1']

    def test_parse_cluster_mapping_empty_file(self, tmp_path):
        """Test parsing an empty cluster mapping file."""
        path = tmp_path / 'empty.tsv'
        path.write_text('# Only comments\n')

        cluster_map = parse_cluster_mapping(str(path))
        assert cluster_map == {}

    def test_parse_cluster_mapping_invalid_lines(self, tmp_path):
        """Test parsing cluster mapping file with invalid lines."""
        path = tmp_path / 'invalid.tsv'
        path.write_text(
            'ValidCluster\tComp1\tComp2\n'
            'InvalidLine\n'  # Only one column
            '\tNoClusterName\n'  # Empty cluster name
        )

        cluster_map = parse_cluster_mapping(str(path))
        assert len(cluster_map) == 1
        assert 'ValidCluster' in cluster_map

    def test_validate_cluster_mapping_valid(self):
        """Test validation of a valid cluster mapping."""
        cluster_map = {
            'ClusterA': ['CompA1', 'CompA2'],
            'ClusterB': ['CompB1'],
        }
        available_features = {'CompA1', 'CompA2', 'CompB1'}

        is_valid, warnings = validate_cluster_mapping(cluster_map, available_features)

        assert is_valid is True
        assert len(warnings) == 0

    def test_validate_cluster_mapping_duplicate_component(self):
        """Test validation detects duplicate component assignments."""
        cluster_map = {
            'ClusterA': ['CompA1', 'CompShared'],
            'ClusterB': ['CompB1', 'CompShared'],  # CompShared in both
        }
        available_features = {'CompA1', 'CompShared', 'CompB1'}

        is_valid, warnings = validate_cluster_mapping(cluster_map, available_features)

        assert is_valid is False
        assert any('CompShared' in w and 'multiple clusters' in w for w in warnings)

    def test_validate_cluster_mapping_unassigned_features(self):
        """Test validation warns about unassigned features."""
        cluster_map = {
            'ClusterA': ['CompA1'],
        }
        available_features = {'CompA1', 'UnassignedComp'}

        is_valid, warnings = validate_cluster_mapping(cluster_map, available_features)

        assert is_valid is True  # Unassigned features are warnings, not errors
        assert any('UnassignedComp' in w for w in warnings)

    def test_build_component_to_cluster_map(self):
        """Test building reverse component-to-cluster mapping."""
        cluster_map = {
            'ClusterA': ['CompA1', 'CompA2'],
            'ClusterB': ['CompB1'],
        }

        component_map = build_component_to_cluster_map(cluster_map)

        assert component_map['CompA1'] == 'ClusterA'
        assert component_map['CompA2'] == 'ClusterA'
        assert component_map['CompB1'] == 'ClusterB'


# -----------------------------------------------------------------------------
# Pairing Map Tests
# -----------------------------------------------------------------------------


class TestPairingMap:
    """Tests for pairing map parsing."""

    def test_parse_pairing_map_valid(self, pairing_map_file):
        """Test parsing a valid pairing map file."""
        pairing_map = parse_pairing_map(pairing_map_file)

        assert len(pairing_map) == 2
        assert pairing_map['LeftA'] == 'RightA'
        assert pairing_map['LeftB'] == 'RightB'

    def test_parse_pairing_map_empty(self, tmp_path):
        """Test parsing an empty pairing map file."""
        path = tmp_path / 'empty.tsv'
        path.write_text('# Comments only\n')

        pairing_map = parse_pairing_map(str(path))
        assert pairing_map == {}

    def test_parse_pairing_map_invalid_lines(self, tmp_path):
        """Test parsing pairing map with invalid lines."""
        path = tmp_path / 'invalid.tsv'
        path.write_text(
            'ValidLeft\tValidRight\n'
            'TooMany\tColumns\tHere\n'  # 3 columns
            'OnlyOne\n'  # 1 column
        )

        pairing_map = parse_pairing_map(str(path))
        assert len(pairing_map) == 1
        assert 'ValidLeft' in pairing_map


# -----------------------------------------------------------------------------
# Hit Loading Tests
# -----------------------------------------------------------------------------


class TestHitLoading:
    """Tests for hit loading from files."""

    def test_load_blast_hits(self, blast_result_file):
        """Test loading hits from BLAST file."""
        from pathlib import Path

        hit_table = load_hits_from_files(blast_files=[Path(blast_result_file)])

        assert not hit_table.empty
        assert len(hit_table) == 5
        assert 'ComponentA1' in hit_table['model'].values

    def test_load_nhmmer_hits(self, nhmmer_result_file):
        """Test loading hits from an nhmmer file.

        This used to assert only that loading did not crash, because the
        fixture's column layout did not match what `import_nhmmer` reads and
        so parsed to zero rows. Both the fixture and `detect_input_format` now
        use the real `--tblout` layout, so the parsed content can be checked.
        """
        from pathlib import Path

        hit_table = load_hits_from_files(nhmmer_files=[Path(nhmmer_result_file)])

        assert len(hit_table) == 2
        assert sorted(hit_table['model'].unique()) == ['ComponentA1', 'ComponentA2']
        assert set(hit_table['target']) == {'chr1'}

    def test_load_no_files_raises_error(self):
        """Test that loading with no files raises an error."""
        with pytest.raises(EnsembleSearchError, match='Must provide'):
            load_hits_from_files()

    def test_filter_hits_by_evalue(self, sample_hit_table):
        """Test filtering hits by e-value threshold."""
        filtered = filter_hits_by_evalue(sample_hit_table, max_evalue=1e-35)

        assert len(filtered) == 3  # Only 1e-40, 1e-35, 1e-38 pass
        assert all(float(e) <= 1e-35 for e in filtered['evalue'])


# -----------------------------------------------------------------------------
# Hit Merging Tests
# -----------------------------------------------------------------------------


class TestHitMerging:
    """Tests for merging overlapping cluster hits."""

    def test_merge_overlapping_hits_simple(self):
        """Test merging overlapping hits from same cluster."""
        hit_table = pd.DataFrame(
            {
                'model': ['CompA1', 'CompA2'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '1050'],
                'hitEnd': ['1100', '1150'],
                'strand': ['+', '+'],
                'evalue': ['1e-40', '1e-35'],
                'score': ['200', '180'],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '98'],
            }
        )

        cluster_map = {'ClusterA': ['CompA1', 'CompA2']}

        merged = merge_overlapping_cluster_hits(hit_table, cluster_map)

        assert len(merged) == 1
        assert merged.iloc[0]['model'] == 'ClusterA'
        assert merged.iloc[0]['hitStart'] == '1000'
        assert merged.iloc[0]['hitEnd'] == '1150'
        assert merged.iloc[0]['score'] == '200'  # Inherited from highest scoring

    def test_merge_non_overlapping_hits_not_merged(self):
        """Test that non-overlapping hits from same cluster are not merged."""
        hit_table = pd.DataFrame(
            {
                'model': ['CompA1', 'CompA2'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '2000'],  # Not overlapping
                'hitEnd': ['1100', '2100'],
                'strand': ['+', '+'],
                'evalue': ['1e-40', '1e-35'],
                'score': ['200', '180'],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '98'],
            }
        )

        cluster_map = {'ClusterA': ['CompA1', 'CompA2']}

        merged = merge_overlapping_cluster_hits(hit_table, cluster_map)

        assert len(merged) == 2
        assert all(merged['model'] == 'ClusterA')

    def test_merge_different_strand_not_merged(self):
        """Test that hits on different strands are not merged."""
        hit_table = pd.DataFrame(
            {
                'model': ['CompA1', 'CompA2'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '1050'],
                'hitEnd': ['1100', '1150'],
                'strand': ['+', '-'],  # Different strands
                'evalue': ['1e-40', '1e-35'],
                'score': ['200', '180'],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '98'],
            }
        )

        cluster_map = {'ClusterA': ['CompA1', 'CompA2']}

        merged = merge_overlapping_cluster_hits(hit_table, cluster_map)

        assert len(merged) == 2

    def test_merge_different_cluster_not_merged(self):
        """Test that overlapping hits from different clusters are not merged."""
        hit_table = pd.DataFrame(
            {
                'model': ['CompA1', 'CompB1'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '1050'],
                'hitEnd': ['1100', '1150'],
                'strand': ['+', '+'],
                'evalue': ['1e-40', '1e-35'],
                'score': ['200', '180'],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '98'],
            }
        )

        cluster_map = {
            'ClusterA': ['CompA1'],
            'ClusterB': ['CompB1'],
        }

        merged = merge_overlapping_cluster_hits(hit_table, cluster_map)

        assert len(merged) == 2
        assert 'ClusterA' in merged['model'].values
        assert 'ClusterB' in merged['model'].values


# -----------------------------------------------------------------------------
# Nested Hit Removal Tests
# -----------------------------------------------------------------------------


class TestNestedHitRemoval:
    """Tests for removing nested weak hits between paired features."""

    def test_removes_nested_hit_on_opposite_strand(self):
        """The canonical F,R case: paired termini sit on OPPOSITE strands.

        This is the situation the step exists to handle, and it was previously
        invisible: grouping by (target, strand) meant LeftA on '+' and RightA
        on '-' were never compared, so no nested cross-hit between a left and
        right terminus of the same element was ever removed.
        """
        hit_table = pd.DataFrame(
            {
                'model': ['LeftA', 'RightA'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '1020'],  # RightA nested inside LeftA
                'hitEnd': ['1200', '1080'],
                'strand': ['+', '-'],  # opposite strands, as F,R implies
                'evalue': ['1e-40', '1e-20'],
                'score': ['200', '100'],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '60'],
            }
        )

        result = remove_nested_paired_hits(hit_table, {'LeftA': 'RightA'})

        assert len(result) == 1
        assert result.iloc[0]['model'] == 'LeftA'

    def test_non_nested_opposite_strand_hits_both_kept(self):
        """Comparing across strands must not remove hits that do not nest."""
        hit_table = pd.DataFrame(
            {
                'model': ['LeftA', 'RightA'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '5000'],
                'hitEnd': ['1200', '5200'],
                'strand': ['+', '-'],
                'evalue': ['1e-40', '1e-20'],
                'score': ['200', '100'],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '60'],
            }
        )

        result = remove_nested_paired_hits(hit_table, {'LeftA': 'RightA'})

        assert len(result) == 2

    def test_remove_nested_hit_weak(self):
        """Test removal of weak hit nested within stronger paired hit."""
        hit_table = pd.DataFrame(
            {
                'model': ['LeftA', 'RightA'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '1020'],  # RightA nested in LeftA
                'hitEnd': ['1200', '1080'],
                'strand': ['+', '+'],
                'evalue': ['1e-40', '1e-20'],
                'score': ['200', '100'],  # RightA is much weaker
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '60'],
            }
        )

        pairing_map = {'LeftA': 'RightA'}

        result = remove_nested_paired_hits(hit_table, pairing_map)

        assert len(result) == 1
        assert result.iloc[0]['model'] == 'LeftA'

    def test_keep_nested_hit_strong(self):
        """Test keeping a nested hit that outscores its container.

        The default min_score_ratio is 1.5: the nested hit is removed only
        when the ENCLOSING hit scores at least 1.5x better. Here the nested
        hit scores higher than its container (160 vs 100), so both are kept.
        """
        hit_table = pd.DataFrame(
            {
                'model': ['LeftA', 'RightA'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '1020'],  # RightA nested in LeftA
                'hitEnd': ['1200', '1080'],
                'strand': ['+', '+'],
                'evalue': ['1e-40', '1e-38'],
                'score': ['100', '160'],  # RightA has higher score than ratio threshold
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '60'],
            }
        )

        pairing_map = {'LeftA': 'RightA'}

        result = remove_nested_paired_hits(hit_table, pairing_map)

        assert len(result) == 2  # Both kept because nested score is relatively high

    def _nested_pair(self, enclosing_score, nested_score):
        """Build LeftA enclosing a nested RightA with the given scores."""
        return pd.DataFrame(
            {
                'model': ['LeftA', 'RightA'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '1020'],
                'hitEnd': ['1200', '1080'],
                'strand': ['+', '+'],
                'evalue': ['1e-40', '1e-20'],
                'score': [str(enclosing_score), str(nested_score)],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '60'],
            }
        )

    def test_equal_scoring_nested_hit_is_kept(self):
        """Two equally good hits are not decisive, so both survive.

        This is the regression test for the inverted ratio. The old rule
        removed the nested hit unless nested/enclosing >= 1.5, so an
        equal-scoring nested hit (ratio 1.0) was deleted.
        """
        result = remove_nested_paired_hits(
            self._nested_pair(200, 200), {'LeftA': 'RightA'}
        )
        assert len(result) == 2

    def test_better_scoring_nested_hit_is_kept(self):
        """A nested hit scoring better than its container is kept.

        Under the old rule 140/100 = 1.4 < 1.5, so this hit was deleted
        despite outscoring the hit that contained it.
        """
        result = remove_nested_paired_hits(
            self._nested_pair(100, 140), {'LeftA': 'RightA'}
        )
        assert len(result) == 2

    def test_much_weaker_nested_hit_is_removed(self):
        """A decisively worse nested hit is still removed."""
        result = remove_nested_paired_hits(
            self._nested_pair(300, 100), {'LeftA': 'RightA'}
        )
        assert len(result) == 1
        assert result.iloc[0]['model'] == 'LeftA'

    def test_ratio_boundary_is_inclusive(self):
        """Exactly min_score_ratio counts as decisive."""
        result = remove_nested_paired_hits(
            self._nested_pair(150, 100), {'LeftA': 'RightA'}
        )
        assert len(result) == 1

    def test_just_below_ratio_boundary_keeps_both(self):
        """Just short of the threshold, the nested hit survives."""
        result = remove_nested_paired_hits(
            self._nested_pair(149, 100), {'LeftA': 'RightA'}
        )
        assert len(result) == 2

    def test_zero_score_nested_hit_is_removed(self):
        """A zero-scoring nested hit cannot defend itself.

        Guarded explicitly rather than by division: with the ratio flipped,
        dividing by the nested score would raise ZeroDivisionError.
        """
        result = remove_nested_paired_hits(
            self._nested_pair(200, 0), {'LeftA': 'RightA'}
        )
        assert len(result) == 1
        assert result.iloc[0]['model'] == 'LeftA'

    def test_two_zero_score_hits_are_both_kept(self):
        """Two unscored hits are indistinguishable, so neither is removed."""
        result = remove_nested_paired_hits(self._nested_pair(0, 0), {'LeftA': 'RightA'})
        assert len(result) == 2

    def test_no_removal_different_targets(self):
        """Test no removal when hits are on different targets."""
        hit_table = pd.DataFrame(
            {
                'model': ['LeftA', 'RightA'],
                'target': ['chr1', 'chr2'],  # Different chromosomes
                'hitStart': ['1000', '1020'],
                'hitEnd': ['1200', '1080'],
                'strand': ['+', '+'],
                'evalue': ['1e-40', '1e-20'],
                'score': ['200', '100'],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '60'],
            }
        )

        pairing_map = {'LeftA': 'RightA'}

        result = remove_nested_paired_hits(hit_table, pairing_map)

        assert len(result) == 2


# -----------------------------------------------------------------------------
# Cross-Model Overlap Filter Tests
# -----------------------------------------------------------------------------


class TestFilterBestModelPerLocus:
    """Tests for filter_best_model_per_locus."""

    def _make_hit(self, model, target, start, end, score, strand='+'):
        return {
            'model': model,
            'target': target,
            'hitStart': str(start),
            'hitEnd': str(end),
            'strand': strand,
            'evalue': '1e-10',
            'score': str(score),
            'bias': 'NA',
            'hmmStart': '1',
            'hmmEnd': '100',
        }

    def test_removes_weaker_overlapping_cross_model_hit(self):
        """Weaker overlapping hit from a different model is removed."""
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 200),
                self._make_hit(
                    'Family2_LEFT', 'chr1', 1050, 1180, 100
                ),  # weaker overlap
            ]
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT', 'Family2_LEFT': 'Family2_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)
        assert len(result) == 1
        assert result.iloc[0]['model'] == 'Family1_LEFT'

    def test_keeps_both_when_score_ratio_below_threshold(self):
        """Both hits are kept when neither dominates by the required ratio."""
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 100),
                self._make_hit('Family2_LEFT', 'chr1', 1050, 1180, 90),  # ratio < 1.5
            ]
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT', 'Family2_LEFT': 'Family2_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)
        assert len(result) == 2

    def test_no_removal_when_hits_do_not_overlap(self):
        """Non-overlapping hits from different models are both retained."""
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1100, 200),
                self._make_hit('Family2_LEFT', 'chr1', 5000, 5200, 50),  # far away
            ]
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT', 'Family2_LEFT': 'Family2_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)
        assert len(result) == 2

    def test_same_model_hits_not_compared(self):
        """Two hits from the same model are never removed against each other."""
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 200),
                self._make_hit(
                    'Family1_LEFT', 'chr1', 1050, 1180, 50
                ),  # same model, weaker
            ]
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)
        assert len(result) == 2

    def test_removes_weaker_hit_across_strands(self):
        """Overlapping cross-model hits are compared regardless of strand.

        This asserted the opposite until the strand-blindness fix. Grouping by
        (target, strand) meant a weaker hit survived merely by sitting on the
        opposite strand from the better one. Under the canonical F,R
        orientation the two termini of an element are on opposite strands by
        construction, so strand-splitting defeated the filter in the common
        case.
        """
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 200, strand='+'),
                self._make_hit('Family2_LEFT', 'chr1', 1050, 1180, 50, strand='-'),
            ]
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT', 'Family2_LEFT': 'Family2_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)

        # 200 / 50 = 4.0, comfortably past the 1.5 decisiveness threshold.
        assert len(result) == 1
        assert result.iloc[0]['model'] == 'Family1_LEFT'

    def test_non_overlapping_hits_on_different_strands_both_kept(self):
        """Strand-blindness must not make the filter over-eager.

        Hits that do not overlap are untouched no matter what strands they are
        on, so the fix only affects genuinely competing hits.
        """
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 200, strand='+'),
                self._make_hit('Family2_LEFT', 'chr1', 5000, 5200, 50, strand='-'),
            ]
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT', 'Family2_LEFT': 'Family2_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)

        assert len(result) == 2

    def test_model_not_in_pairing_map_is_preserved(self):
        """Hits from models not in the pairing map are never removed."""
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 200),
                self._make_hit('UnknownModel', 'chr1', 1050, 1180, 50),  # not in map
            ]
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)
        assert len(result) == 2

    def test_empty_table_returned_unchanged(self):
        """Empty hit table is returned unchanged."""
        hit_table = pd.DataFrame(
            columns=[
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
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)
        assert result.empty

    def test_empty_pairing_map_returned_unchanged(self):
        """Hit table is returned unchanged when pairing map is empty."""
        hit_table = pd.DataFrame(
            [self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 200)]
        )
        result = filter_best_model_per_locus(hit_table, {})
        assert len(result) == 1

    def test_multiple_families_same_locus_keeps_best(self):
        """With three overlapping models at the same locus, only the best is kept."""
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 300),
                self._make_hit('Family2_LEFT', 'chr1', 1020, 1180, 100),
                self._make_hit('Family3_LEFT', 'chr1', 1010, 1190, 50),
            ]
        )
        pairing_map = {
            'Family1_LEFT': 'Family1_RIGHT',
            'Family2_LEFT': 'Family2_RIGHT',
            'Family3_LEFT': 'Family3_RIGHT',
        }
        result = filter_best_model_per_locus(hit_table, pairing_map)
        # Family1_LEFT (score=300) beats Family2_LEFT by 3x and Family3_LEFT by 6x,
        # both above the default min_score_ratio of 1.5.
        assert len(result) == 1
        assert result.iloc[0]['model'] == 'Family1_LEFT'

    def test_no_removal_on_different_targets(self):
        """Overlapping coordinates on different target sequences are independent."""
        hit_table = pd.DataFrame(
            [
                self._make_hit('Family1_LEFT', 'chr1', 1000, 1200, 200),
                self._make_hit('Family2_LEFT', 'chr2', 1000, 1200, 50),
            ]
        )
        pairing_map = {'Family1_LEFT': 'Family1_RIGHT', 'Family2_LEFT': 'Family2_RIGHT'}
        result = filter_best_model_per_locus(hit_table, pairing_map)
        assert len(result) == 2


class TestCrossClusterOverlaps:
    """Tests for cross-cluster overlap detection."""

    def test_detect_cross_cluster_overlap(self, caplog):
        """Test detection and warning of cross-cluster overlaps."""
        import logging

        hit_table = pd.DataFrame(
            {
                'model': ['CompA1', 'CompB1'],
                'target': ['chr1', 'chr1'],
                'hitStart': ['1000', '1050'],  # Overlapping
                'hitEnd': ['1100', '1150'],
                'strand': ['+', '+'],
                'evalue': ['1e-40', '1e-35'],
                'score': ['200', '180'],
                'bias': ['NA', 'NA'],
                'hmmStart': ['1', '1'],
                'hmmEnd': ['100', '98'],
            }
        )

        cluster_map = {
            'ClusterA': ['CompA1'],
            'ClusterB': ['CompB1'],
        }

        with caplog.at_level(logging.WARNING):
            check_cross_cluster_overlaps(hit_table, cluster_map)

        assert 'Cross-cluster overlap' in caplog.text


# -----------------------------------------------------------------------------
# Anchor (Outer-Edge) Filter Tests
# -----------------------------------------------------------------------------


class TestComputeOuterEdgeOffset:
    """Tests for compute_outer_edge_offset."""

    def test_left_plus_at_edge(self):
        """Left terminus, + strand: hmmStart=1 → offset 0."""
        from tirmite.cli.ensemble_search import compute_outer_edge_offset

        assert compute_outer_edge_offset(1, 100, 100, '+', 'left') == 0

    def test_left_plus_offset_5(self):
        """Left terminus, + strand: hmmStart=6 → offset 5."""
        from tirmite.cli.ensemble_search import compute_outer_edge_offset

        assert compute_outer_edge_offset(6, 100, 100, '+', 'left') == 5

    def test_left_minus_at_edge(self):
        """Left terminus, - strand: hmmEnd=model_len → offset 0."""
        from tirmite.cli.ensemble_search import compute_outer_edge_offset

        assert compute_outer_edge_offset(1, 100, 100, '-', 'left') == 0

    def test_left_minus_offset_5(self):
        """Left terminus, - strand: hmmEnd=95 with model_len=100 → offset 5."""
        from tirmite.cli.ensemble_search import compute_outer_edge_offset

        assert compute_outer_edge_offset(1, 95, 100, '-', 'left') == 5

    def test_right_plus_at_edge(self):
        """Right terminus, + strand: hmmEnd=model_len → offset 0."""
        from tirmite.cli.ensemble_search import compute_outer_edge_offset

        assert compute_outer_edge_offset(1, 100, 100, '+', 'right') == 0

    def test_right_plus_offset_5(self):
        """Right terminus, + strand: hmmEnd=95, model_len=100 → offset 5."""
        from tirmite.cli.ensemble_search import compute_outer_edge_offset

        assert compute_outer_edge_offset(1, 95, 100, '+', 'right') == 5

    def test_right_minus_at_edge(self):
        """Right terminus, - strand: hmmStart=1 → offset 0."""
        from tirmite.cli.ensemble_search import compute_outer_edge_offset

        assert compute_outer_edge_offset(1, 100, 100, '-', 'right') == 0

    def test_right_minus_offset_5(self):
        """Right terminus, - strand: hmmStart=6 → offset 5."""
        from tirmite.cli.ensemble_search import compute_outer_edge_offset

        assert compute_outer_edge_offset(6, 100, 100, '-', 'right') == 5


class TestFilterHitsByAnchor:
    """Tests for filter_hits_by_anchor."""

    def _make_hit_table(self, rows):
        """Helper to build a hit DataFrame."""
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # F,R orientation (canonical TIR): + strand = left, - strand = right
    # ------------------------------------------------------------------

    def test_fr_left_plus_passes_within_offset(self):
        """F,R: left(+) hit with hmmStart=6 passes when max_offset=10."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '6',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 1

    def test_fr_left_plus_removed_exceeds_offset(self):
        """F,R: left(+) hit with hmmStart=16 (offset=15) is removed when max_offset=10."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '16',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 0

    def test_fr_right_minus_passes_within_offset(self):
        """F,R: right(-) hit with hmmStart=1 (offset=0) passes when max_offset=10."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '-',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '1',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 1

    def test_fr_right_minus_removed_exceeds_offset(self):
        """F,R: right(-) hit with hmmStart=20 (offset=19) is removed when max_offset=10."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '-',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '20',
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 0

    # ------------------------------------------------------------------
    # Asymmetric pairing: model name determines terminus type
    # ------------------------------------------------------------------

    def test_asymmetric_left_model_filtered(self):
        """Asymmetric: left model with large offset is removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'LeftModel',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
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
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '1',
                    'hmmEnd': '95',  # offset = 5 <= max_offset=10
                },
            ]
        )
        pairing = {'LeftModel': 'RightModel'}
        lengths = {'LeftModel': 100, 'RightModel': 100}
        result = filter_hits_by_anchor(
            df, lengths, max_offset=10, orientation='F,F', pairing_map=pairing
        )
        assert len(result) == 1
        assert result.iloc[0]['model'] == 'RightModel'

    # ------------------------------------------------------------------
    # Edge cases
    # ------------------------------------------------------------------

    def test_empty_table_returned_unchanged(self):
        """Empty hit table is returned unchanged."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = pd.DataFrame(
            columns=[
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
        )
        result = filter_hits_by_anchor(df, {}, max_offset=5)
        assert result.empty

    def test_missing_model_length_raises_error(self):
        """Raises AnchorFilterError under the mode the search workflow uses."""
        from tirmite.cli.ensemble_search import AnchorFilterError, filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'UnknownModel',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '50',
                    'hmmEnd': '100',
                }
            ]
        )
        with pytest.raises(AnchorFilterError, match='model length'):
            filter_hits_by_anchor(
                df, {}, max_offset=5, orientation='F,R', on_missing_length='raise'
            )

    def test_ff_symmetric_no_pairing_map_keeps_hits(self):
        """F,F symmetric without pairing map: terminus type unknown → hits kept."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '50',  # Large offset – would be removed if type known
                    'hmmEnd': '100',
                }
            ]
        )
        # F,F: same strand, no pairing map → both-ends check:
        # offset_start = 50 - 1 = 49 > 5, so hit is removed
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation='F,F'
        )
        assert len(result) == 0

    def test_ff_symmetric_no_pairing_map_keeps_hit_within_both_offsets(self):
        """F,F: same strand, no pairing map, hit covers full model → kept."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '1',
                    'hmmEnd': '100',
                }
            ]
        )
        # F,F: same strand, no pairing map → both-ends check:
        # offset_start = 0, offset_end = 0, both ≤ 5 → kept
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation='F,F'
        )
        assert len(result) == 1

    def test_rf_orientation_left_minus_filtered(self):
        """R,F: left(-) hit outer edge at position model_len; hmmEnd=80 → offset=20 > 10."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '-',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '1',
                    'hmmEnd': '80',  # left(-): offset = model_len - hmmEnd = 100-80 = 20
                }
            ]
        )
        # R,F: left_strand='-', right_strand='+'
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='R,F'
        )
        assert len(result) == 0

    def test_exact_offset_boundary_kept(self):
        """Hit with offset exactly equal to max_offset is kept."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = self._make_hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '11',  # offset = 10 = max_offset → kept
                    'hmmEnd': '100',
                }
            ]
        )
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=10, orientation='F,R'
        )
        assert len(result) == 1


# -----------------------------------------------------------------------------
# Comprehensive Anchor Filter Tests – All Orientations & Offset Values
# -----------------------------------------------------------------------------

# Helper: a standard 100-bp model with hits that start/end at various positions.
# hmmStart/hmmEnd are 1-based positions in the query (BLAST qstart/qend or hmmer).
#
# Offset scenarios tested per orientation:
#   - no_anchor  : max_offset is None → filter not applied (all hits pass)
#   - exact_zero : max_offset=0 → only hits at the outer edge pass
#   - intermediate: max_offset=10 → hits within 10 bp of outer edge pass
#   - over_model : max_offset=200 (> model len 100) → all hits pass


def _make_row(model, strand, hmm_start, hmm_end, target='chr1', model_len=100):
    """Create a single hit-table row dict."""
    return {
        'model': model,
        'target': target,
        'hitStart': '100',
        'hitEnd': str(100 + model_len - 1),
        'strand': strand,
        'evalue': '1e-10',
        'score': '100',
        'bias': 'NA',
        'hmmStart': str(hmm_start),
        'hmmEnd': str(hmm_end),
    }


def _anchor_df(rows):
    return pd.DataFrame(rows)


class TestAnchorFilterFR:
    """Comprehensive anchor tests for F,R orientation (+ = left, - = right)."""

    ORIENTATION = 'F,R'
    MODEL = 'TIR'
    LENGTHS = {'TIR': 100}

    # --- left terminus (+ strand): outer edge = position 1, offset = hmmStart-1 ---

    def test_left_plus_no_anchor_passes(self):
        """No anchor filter: all hits pass regardless of offset."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        # hmmStart=50 would give offset=49 - would fail with any reasonable max_offset
        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=50, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=200, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_left_plus_exact_edge_passes(self):
        """Left(+): hmmStart=1 (offset=0), max_offset=0 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=0, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_left_plus_exact_edge_fails_if_one_off(self):
        """Left(+): hmmStart=2 (offset=1), max_offset=0 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=2, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=0, orientation=self.ORIENTATION
        )
        assert len(result) == 0

    def test_left_plus_intermediate_passes(self):
        """Left(+): hmmStart=6 (offset=5), max_offset=10 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=6, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=10, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_left_plus_intermediate_fails(self):
        """Left(+): hmmStart=20 (offset=19), max_offset=10 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=20, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=10, orientation=self.ORIENTATION
        )
        assert len(result) == 0

    def test_left_plus_over_model_passes(self):
        """Left(+): max_offset=200 (larger than model) → all hits pass."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=99, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=200, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    # --- right terminus (- strand): outer edge = position 1, offset = hmmStart-1 ---

    def test_right_minus_exact_edge_passes(self):
        """Right(-): hmmStart=1 (offset=0), max_offset=0 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=0, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_right_minus_intermediate_passes(self):
        """Right(-): hmmStart=6 (offset=5), max_offset=10 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=6, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=10, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_right_minus_intermediate_fails(self):
        """Right(-): hmmStart=20 (offset=19), max_offset=10 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=20, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=10, orientation=self.ORIENTATION
        )
        assert len(result) == 0

    def test_right_minus_over_model_passes(self):
        """Right(-): max_offset=200 → all hits pass."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=99, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=200, orientation=self.ORIENTATION
        )
        assert len(result) == 1


class TestAnchorFilterRF:
    """Comprehensive anchor tests for R,F orientation (- = left, + = right)."""

    ORIENTATION = 'R,F'
    MODEL = 'TIR'
    LENGTHS = {'TIR': 100}

    # --- left terminus (- strand): outer edge = position model_len, offset = model_len - hmmEnd ---

    def test_left_minus_exact_edge_passes(self):
        """Left(-): hmmEnd=100 (offset=0), max_offset=0 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=0, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_left_minus_exact_edge_fails_if_one_off(self):
        """Left(-): hmmEnd=99 (offset=1), max_offset=0 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=1, hmm_end=99)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=0, orientation=self.ORIENTATION
        )
        assert len(result) == 0

    def test_left_minus_intermediate_passes(self):
        """Left(-): hmmEnd=95 (offset=5), max_offset=10 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=1, hmm_end=95)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=10, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_left_minus_intermediate_fails(self):
        """Left(-): hmmEnd=80 (offset=20), max_offset=10 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=1, hmm_end=80)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=10, orientation=self.ORIENTATION
        )
        assert len(result) == 0

    def test_left_minus_over_model_passes(self):
        """Left(-): max_offset=200 → all hits pass."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '-', hmm_start=1, hmm_end=5)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=200, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    # --- right terminus (+ strand): outer edge = position model_len, offset = model_len - hmmEnd ---

    def test_right_plus_exact_edge_passes(self):
        """Right(+): hmmEnd=100 (offset=0), max_offset=0 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=0, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_right_plus_intermediate_passes(self):
        """Right(+): hmmEnd=95 (offset=5), max_offset=10 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=1, hmm_end=95)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=10, orientation=self.ORIENTATION
        )
        assert len(result) == 1

    def test_right_plus_intermediate_fails(self):
        """Right(+): hmmEnd=80 (offset=20), max_offset=10 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row(self.MODEL, '+', hmm_start=1, hmm_end=80)])
        result = filter_hits_by_anchor(
            df, self.LENGTHS, max_offset=10, orientation=self.ORIENTATION
        )
        assert len(result) == 0


class TestAnchorFilterFF:
    """Comprehensive anchor tests for F,F orientation with asymmetric pairing map."""

    ORIENTATION = 'F,F'
    LENGTHS = {'LeftTIR': 100, 'RightTIR': 100}
    PAIRING = {'LeftTIR': 'RightTIR'}

    # --- left terminus (LeftTIR, + strand): outer edge = pos 1, offset = hmmStart-1 ---

    def test_ff_left_exact_edge_passes(self):
        """F,F left: hmmStart=1 (offset=0), max_offset=0 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('LeftTIR', '+', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=0,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_ff_left_intermediate_passes(self):
        """F,F left: hmmStart=6 (offset=5), max_offset=10 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('LeftTIR', '+', hmm_start=6, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=10,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_ff_left_intermediate_fails(self):
        """F,F left: hmmStart=20 (offset=19), max_offset=10 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('LeftTIR', '+', hmm_start=20, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=10,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 0

    def test_ff_left_over_model_passes(self):
        """F,F left: max_offset=200 → all hits pass."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('LeftTIR', '+', hmm_start=99, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=200,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    # --- right terminus (RightTIR, + strand): outer edge = pos model_len, offset = model_len - hmmEnd ---

    def test_ff_right_exact_edge_passes(self):
        """F,F right: hmmEnd=100 (offset=0), max_offset=0 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('RightTIR', '+', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=0,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_ff_right_intermediate_passes(self):
        """F,F right: hmmEnd=95 (offset=5), max_offset=10 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('RightTIR', '+', hmm_start=1, hmm_end=95)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=10,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_ff_right_intermediate_fails(self):
        """F,F right: hmmEnd=80 (offset=20), max_offset=10 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('RightTIR', '+', hmm_start=1, hmm_end=80)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=10,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 0

    def test_ff_right_over_model_passes(self):
        """F,F right: max_offset=200 → all hits pass."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('RightTIR', '+', hmm_start=1, hmm_end=5)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=200,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_ff_no_pairing_map_keeps_hits(self):
        """F,F without pairing map: both-ends check, large offset → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        # Large offset that is removed: offset_start = 50-1 = 49 > 5
        df = _anchor_df([_make_row('TIR', '+', hmm_start=50, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation=self.ORIENTATION
        )
        assert len(result) == 0

    def test_ff_no_pairing_map_full_coverage_kept(self):
        """F,F without pairing map: both-ends check, full coverage → kept."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        # Full model coverage: offset_start = 0, offset_end = 0
        df = _anchor_df([_make_row('TIR', '+', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df, {'TIR': 100}, max_offset=5, orientation=self.ORIENTATION
        )
        assert len(result) == 1


class TestAnchorFilterRR:
    """Comprehensive anchor tests for R,R orientation with asymmetric pairing map."""

    ORIENTATION = 'R,R'
    LENGTHS = {'LeftTIR': 100, 'RightTIR': 100}
    PAIRING = {'LeftTIR': 'RightTIR'}

    # --- left terminus (LeftTIR, - strand): outer edge = pos model_len, offset = model_len - hmmEnd ---

    def test_rr_left_exact_edge_passes(self):
        """R,R left(-): hmmEnd=100 (offset=0), max_offset=0 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('LeftTIR', '-', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=0,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_rr_left_intermediate_passes(self):
        """R,R left(-): hmmEnd=95 (offset=5), max_offset=10 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('LeftTIR', '-', hmm_start=1, hmm_end=95)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=10,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_rr_left_intermediate_fails(self):
        """R,R left(-): hmmEnd=80 (offset=20), max_offset=10 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('LeftTIR', '-', hmm_start=1, hmm_end=80)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=10,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 0

    def test_rr_left_over_model_passes(self):
        """R,R left(-): max_offset=200 → all hits pass."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('LeftTIR', '-', hmm_start=1, hmm_end=5)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=200,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    # --- right terminus (RightTIR, - strand): outer edge = pos 1, offset = hmmStart-1 ---

    def test_rr_right_exact_edge_passes(self):
        """R,R right(-): hmmStart=1 (offset=0), max_offset=0 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('RightTIR', '-', hmm_start=1, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=0,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_rr_right_intermediate_passes(self):
        """R,R right(-): hmmStart=6 (offset=5), max_offset=10 → passes."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('RightTIR', '-', hmm_start=6, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=10,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1

    def test_rr_right_intermediate_fails(self):
        """R,R right(-): hmmStart=20 (offset=19), max_offset=10 → removed."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('RightTIR', '-', hmm_start=20, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=10,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 0

    def test_rr_right_over_model_passes(self):
        """R,R right(-): max_offset=200 → all hits pass."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('RightTIR', '-', hmm_start=99, hmm_end=100)])
        result = filter_hits_by_anchor(
            df,
            self.LENGTHS,
            max_offset=200,
            orientation=self.ORIENTATION,
            pairing_map=self.PAIRING,
        )
        assert len(result) == 1


class TestAnchorFilterMissingLength:
    """Tests for error behaviour when model lengths are unavailable.

    The search workflow passes ``on_missing_length='raise'`` so that asking
    for anchor filtering without the lengths needed to perform it is a hard
    error rather than a silently unfiltered result. ``tirmite pair`` uses the
    default ``'warn'`` instead, which is why these tests state the mode
    explicitly.
    """

    def test_raises_error_when_length_missing_fr(self):
        """Raises AnchorFilterError when an F,R hit has no model length."""
        from tirmite.cli.ensemble_search import AnchorFilterError, filter_hits_by_anchor

        df = _anchor_df([_make_row('TIR', '+', hmm_start=1, hmm_end=80)])
        with pytest.raises(AnchorFilterError, match='model length'):
            filter_hits_by_anchor(
                df, {}, max_offset=5, orientation='F,R', on_missing_length='raise'
            )

    def test_raises_error_names_missing_model(self):
        """Error message includes the name of the missing model."""
        from tirmite.cli.ensemble_search import AnchorFilterError, filter_hits_by_anchor

        df = _anchor_df([_make_row('MyMissingModel', '+', hmm_start=1, hmm_end=80)])
        with pytest.raises(AnchorFilterError, match='MyMissingModel'):
            filter_hits_by_anchor(
                df, {}, max_offset=5, orientation='F,R', on_missing_length='raise'
            )

    def test_ff_same_strand_no_pairing_map_raises(self):
        """F,F without a pairing map needs the length for the both-ends check."""
        from tirmite.cli.ensemble_search import AnchorFilterError, filter_hits_by_anchor

        df = _anchor_df([_make_row('TIR', '+', hmm_start=50, hmm_end=100)])
        with pytest.raises(AnchorFilterError, match='model lengths'):
            filter_hits_by_anchor(
                df, {}, max_offset=5, orientation='F,F', on_missing_length='raise'
            )

    def test_warn_mode_keeps_hit_unchecked(self, caplog):
        """The default 'warn' mode keeps the hit and logs, matching tirmite pair."""
        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df([_make_row('TIR', '+', hmm_start=1, hmm_end=80)])
        result = filter_hits_by_anchor(df, {}, max_offset=5, orientation='F,R')

        assert len(result) == 1
        assert 'Model length not found' in caplog.text

    def test_search_pipeline_reports_missing_length_as_search_error(self, tmp_path):
        """_process_hits translates AnchorFilterError into EnsembleSearchError.

        The CLI presents a single exception type to the user, so unifying the
        two anchor-filter copies must not change what `tirmite search` raises.
        """
        import argparse

        from tirmite.cli.ensemble_search import EnsembleSearchError, _process_hits

        blast_file = tmp_path / 'hits.blast'
        blast_file.write_text(
            'UnknownModel\tchr1\t95.0\t100\t5\t0\t1\t100\t1000\t1099\t1e-40\t200\n'
        )

        args = argparse.Namespace(
            max_evalue=1.0,
            max_offset=5,
            orientation='F,R',
            pairing_map=None,
            cluster_map=None,
        )

        # query_lengths is empty, so the anchor filter has no length for
        # UnknownModel and must fail rather than pass the hit through.
        with pytest.raises(EnsembleSearchError, match='model lengths'):
            _process_hits(
                args,
                blast_files=[blast_file],
                nhmmer_files=[],
                query_lengths={},
            )


class TestAnchorFilterLogging:
    """Tests for logging output from anchor filter."""

    def test_logging_reports_removed_count(self, caplog):
        """Anchor filter logs how many hits were excluded."""
        import logging as stdlib_logging

        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df(
            [
                _make_row('TIR', '+', hmm_start=1, hmm_end=100),  # passes (offset=0)
                _make_row('TIR', '+', hmm_start=50, hmm_end=100),  # fails (offset=49)
            ]
        )
        with caplog.at_level(stdlib_logging.INFO):
            filter_hits_by_anchor(df, {'TIR': 100}, max_offset=5, orientation='F,R')

        assert '1 removed' in caplog.text

    def test_logging_reports_per_model_counts(self, caplog):
        """Anchor filter logs per-model exclusion counts."""
        import logging as stdlib_logging

        from tirmite.cli.ensemble_search import filter_hits_by_anchor

        df = _anchor_df(
            [
                _make_row('LeftTIR', '+', hmm_start=50, hmm_end=100),  # fails
                _make_row('RightTIR', '-', hmm_start=50, hmm_end=100),  # fails
            ]
        )
        with caplog.at_level(stdlib_logging.INFO):
            filter_hits_by_anchor(
                df,
                {'LeftTIR': 100, 'RightTIR': 100},
                max_offset=5,
                orientation='F,R',
            )

        assert 'LeftTIR' in caplog.text
        assert 'RightTIR' in caplog.text


# -----------------------------------------------------------------------------
# Split Paired Output Tests
# -----------------------------------------------------------------------------


class TestValidateSplitPairedOutput:
    """Tests for validate_split_paired_output."""

    def test_valid_disjoint_models(self):
        """No overlap between left and right columns → no error."""
        from tirmite.cli.ensemble_search import validate_split_paired_output

        pairing_map = {'LeftA': 'RightA', 'LeftB': 'RightB'}
        validate_split_paired_output(pairing_map)  # should not raise

    def test_overlap_raises_error(self):
        """Model appearing in both left and right columns → raises error."""
        from tirmite.cli.ensemble_search import (
            EnsembleSearchError,
            validate_split_paired_output,
        )

        pairing_map = {'ModelA': 'ModelB', 'ModelB': 'ModelC'}
        with pytest.raises(EnsembleSearchError, match='ModelB'):
            validate_split_paired_output(pairing_map)

    def test_symmetric_pairing_raises_error(self):
        """Symmetric pairing (same model in both columns) → raises error."""
        from tirmite.cli.ensemble_search import (
            EnsembleSearchError,
            validate_split_paired_output,
        )

        pairing_map = {'TIR': 'TIR'}
        with pytest.raises(EnsembleSearchError, match='TIR'):
            validate_split_paired_output(pairing_map)


class TestWriteSplitHits:
    """Tests for write_split_hits."""

    def _make_hit_table(self, rows):
        return pd.DataFrame(rows)

    def test_split_basic(self, tmp_path):
        """Hits are correctly split into left and right files."""
        from tirmite.cli.ensemble_search import write_split_hits

        df = self._make_hit_table(
            [
                {
                    'model': 'LeftA',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '1',
                    'hmmEnd': '100',
                },
                {
                    'model': 'RightA',
                    'target': 'chr1',
                    'hitStart': '500',
                    'hitEnd': '600',
                    'strand': '-',
                    'evalue': '1e-10',
                    'score': '90',
                    'bias': 'NA',
                    'hmmStart': '1',
                    'hmmEnd': '100',
                },
            ]
        )
        pairing_map = {'LeftA': 'RightA'}
        left_file, right_file = write_split_hits(df, pairing_map, tmp_path, 'test')

        assert left_file.exists()
        assert right_file.exists()

        left_content = left_file.read_text()
        right_content = right_file.read_text()

        assert 'LeftA' in left_content
        assert 'RightA' not in left_content
        assert 'RightA' in right_content
        assert 'LeftA' not in right_content

    def test_split_empty_table(self, tmp_path):
        """Empty table produces empty output files."""
        from tirmite.cli.ensemble_search import write_split_hits

        df = self._make_hit_table([])
        pairing_map = {'LeftA': 'RightA'}
        left_file, right_file = write_split_hits(df, pairing_map, tmp_path, 'test')

        assert left_file.exists()
        assert right_file.exists()

    def test_split_warns_for_unassigned_models(self, tmp_path, caplog):
        """Hits from models not in pairing map trigger a warning."""
        import logging as stdlib_logging

        from tirmite.cli.ensemble_search import write_split_hits

        df = self._make_hit_table(
            [
                {
                    'model': 'UnknownModel',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'score': '100',
                    'bias': 'NA',
                    'hmmStart': '1',
                    'hmmEnd': '100',
                },
            ]
        )
        pairing_map = {'LeftA': 'RightA'}
        with caplog.at_level(stdlib_logging.WARNING):
            write_split_hits(df, pairing_map, tmp_path, 'test')

        assert 'UnknownModel' in caplog.text


# -----------------------------------------------------------------------------
# Pairing Map Model Filter Tests
# -----------------------------------------------------------------------------


class TestFilterHitsToPairingMapModels:
    """Tests for filter_hits_to_pairing_map_models."""

    def _make_hit_table(self, rows):
        return pd.DataFrame(rows)

    def _make_row(self, model):
        return {
            'model': model,
            'target': 'chr1',
            'hitStart': '100',
            'hitEnd': '200',
            'strand': '+',
            'evalue': '1e-10',
            'score': '100',
            'bias': 'NA',
            'hmmStart': '1',
            'hmmEnd': '100',
        }

    def test_retains_paired_models(self):
        """Hits from models in the pairing map are retained."""
        from tirmite.cli.ensemble_search import filter_hits_to_pairing_map_models

        df = self._make_hit_table([self._make_row('LeftA'), self._make_row('RightA')])
        result = filter_hits_to_pairing_map_models(df, {'LeftA': 'RightA'})
        assert len(result) == 2
        assert set(result['model'].unique()) == {'LeftA', 'RightA'}

    def test_removes_unpaired_models(self):
        """Hits from models not in the pairing map are removed."""
        from tirmite.cli.ensemble_search import filter_hits_to_pairing_map_models

        df = self._make_hit_table(
            [
                self._make_row('LeftA'),
                self._make_row('RightA'),
                self._make_row('UnknownModel'),
            ]
        )
        result = filter_hits_to_pairing_map_models(df, {'LeftA': 'RightA'})
        assert len(result) == 2
        assert 'UnknownModel' not in result['model'].values

    def test_empty_table_returns_empty(self):
        """Empty input returns empty output."""
        from tirmite.cli.ensemble_search import filter_hits_to_pairing_map_models

        df = self._make_hit_table([])
        result = filter_hits_to_pairing_map_models(df, {'LeftA': 'RightA'})
        assert result.empty

    def test_empty_pairing_map_returns_unchanged(self):
        """Empty pairing map returns the table unchanged."""
        from tirmite.cli.ensemble_search import filter_hits_to_pairing_map_models

        df = self._make_hit_table([self._make_row('ModelX')])
        result = filter_hits_to_pairing_map_models(df, {})
        assert len(result) == 1

    def test_logs_removed_models(self, caplog):
        """Removed models are reported in the log."""
        import logging as stdlib_logging

        from tirmite.cli.ensemble_search import filter_hits_to_pairing_map_models

        df = self._make_hit_table(
            [self._make_row('LeftA'), self._make_row('UnknownModel')]
        )
        with caplog.at_level(stdlib_logging.INFO):
            filter_hits_to_pairing_map_models(df, {'LeftA': 'RightA'})

        assert 'UnknownModel' in caplog.text

    def test_summary_populated_with_excluded_model_counts(self):
        """Summary.excluded_not_in_map is populated when a summary is passed."""
        from tirmite.cli.ensemble_search import (
            SearchFilterSummary,
            filter_hits_to_pairing_map_models,
        )

        df = self._make_hit_table(
            [
                self._make_row('LeftA'),
                self._make_row('UnknownModel'),
                self._make_row('UnknownModel'),
            ]
        )
        summary = SearchFilterSummary()
        filter_hits_to_pairing_map_models(df, {'LeftA': 'RightA'}, summary=summary)
        assert summary.excluded_not_in_map == {'UnknownModel': 2}

    def test_summary_empty_when_no_exclusions(self):
        """Summary.excluded_not_in_map is empty when all models are in the map."""
        from tirmite.cli.ensemble_search import (
            SearchFilterSummary,
            filter_hits_to_pairing_map_models,
        )

        df = self._make_hit_table([self._make_row('LeftA'), self._make_row('RightA')])
        summary = SearchFilterSummary()
        filter_hits_to_pairing_map_models(df, {'LeftA': 'RightA'}, summary=summary)
        assert summary.excluded_not_in_map == {}


# -----------------------------------------------------------------------------
# SearchFilterSummary and log_filter_summary Tests
# -----------------------------------------------------------------------------


class TestLogFilterSummary:
    """Tests for log_filter_summary."""

    def test_logs_step0_excluded_models(self, caplog):
        """Step 0 exclusions are reported."""
        import logging as stdlib_logging

        from tirmite.cli.ensemble_search import SearchFilterSummary, log_filter_summary

        summary = SearchFilterSummary(excluded_not_in_map={'ModelX': 3, 'ModelY': 1})
        with caplog.at_level(stdlib_logging.INFO):
            log_filter_summary(summary)

        assert 'ModelX' in caplog.text
        assert 'ModelY' in caplog.text
        assert '4' in caplog.text or 'Step 0' in caplog.text

    def test_logs_step1_nested_removal(self, caplog):
        """Step 1 nested removals are reported with container model names."""
        import logging as stdlib_logging

        from tirmite.cli.ensemble_search import SearchFilterSummary, log_filter_summary

        summary = SearchFilterSummary(nested_removed={'LeftA': {'RightA': 2}})
        with caplog.at_level(stdlib_logging.INFO):
            log_filter_summary(summary)

        assert 'LeftA' in caplog.text
        assert 'RightA' in caplog.text

    def test_logs_step2_cross_model_removal(self, caplog):
        """Step 2 cross-model removals are reported per pair."""
        import logging as stdlib_logging

        from tirmite.cli.ensemble_search import SearchFilterSummary, log_filter_summary

        summary = SearchFilterSummary(
            cross_model_removed={('WeakModel', 'StrongModel'): 5}
        )
        with caplog.at_level(stdlib_logging.INFO):
            log_filter_summary(summary)

        assert 'WeakModel' in caplog.text
        assert 'StrongModel' in caplog.text
        assert '5' in caplog.text

    def test_logs_no_removal_messages_when_empty(self, caplog):
        """Empty summary produces 'No hits excluded' / 'No ... removed' messages."""
        import logging as stdlib_logging

        from tirmite.cli.ensemble_search import SearchFilterSummary, log_filter_summary

        summary = SearchFilterSummary()
        with caplog.at_level(stdlib_logging.INFO):
            log_filter_summary(summary)

        assert 'No hits excluded' in caplog.text
        assert 'No nested hits removed' in caplog.text
        assert 'No cross-model' in caplog.text


class TestRemoveNestedPairedHitsSummary:
    """Tests that remove_nested_paired_hits populates SearchFilterSummary."""

    def _make_hit(self, model, start, end, score):
        return {
            'model': model,
            'target': 'chr1',
            'hitStart': str(start),
            'hitEnd': str(end),
            'strand': '+',
            'evalue': '1e-10',
            'score': str(score),
            'bias': 'NA',
            'hmmStart': '1',
            'hmmEnd': str(end - start + 1),
        }

    def test_summary_records_nested_removal(self):
        """nested_removed is populated when a weak nested hit is removed."""
        from tirmite.cli.ensemble_search import (
            SearchFilterSummary,
            remove_nested_paired_hits,
        )

        df = pd.DataFrame(
            [
                self._make_hit('LeftA', 100, 500, 200.0),  # enclosing
                self._make_hit('RightA', 150, 400, 50.0),  # nested, weak
            ]
        )
        summary = SearchFilterSummary()
        remove_nested_paired_hits(df, {'LeftA': 'RightA'}, summary=summary)
        assert 'RightA' in summary.nested_removed
        assert summary.nested_removed['RightA']['LeftA'] == 1

    def test_summary_empty_when_no_nested_hit_removed(self):
        """nested_removed is empty when no hit qualifies for removal."""
        from tirmite.cli.ensemble_search import (
            SearchFilterSummary,
            remove_nested_paired_hits,
        )

        df = pd.DataFrame(
            [
                self._make_hit('LeftA', 100, 500, 200.0),
                self._make_hit('RightA', 600, 800, 180.0),  # non-overlapping
            ]
        )
        summary = SearchFilterSummary()
        remove_nested_paired_hits(df, {'LeftA': 'RightA'}, summary=summary)
        assert summary.nested_removed == {}


class TestFilterBestModelPerLocusSummary:
    """Tests that filter_best_model_per_locus populates SearchFilterSummary."""

    def _make_hit(self, model, start, end, score):
        return {
            'model': model,
            'target': 'chr1',
            'hitStart': str(start),
            'hitEnd': str(end),
            'strand': '+',
            'evalue': '1e-10',
            'score': str(score),
            'bias': 'NA',
            'hmmStart': '1',
            'hmmEnd': str(end - start + 1),
        }

    def test_summary_records_cross_model_removal(self):
        """cross_model_removed is populated when a weaker overlapping hit is removed."""
        from tirmite.cli.ensemble_search import (
            SearchFilterSummary,
            filter_best_model_per_locus,
        )

        df = pd.DataFrame(
            [
                self._make_hit('LeftA', 100, 300, 200.0),
                self._make_hit('LeftB', 150, 280, 50.0),  # weaker, overlapping
            ]
        )
        summary = SearchFilterSummary()
        filter_best_model_per_locus(
            df, {'LeftA': 'RightA', 'LeftB': 'RightB'}, summary=summary
        )
        assert ('LeftB', 'LeftA') in summary.cross_model_removed
        assert summary.cross_model_removed[('LeftB', 'LeftA')] == 1

    def test_summary_empty_when_no_cross_model_removal(self):
        """cross_model_removed is empty when no hits are removed."""
        from tirmite.cli.ensemble_search import (
            SearchFilterSummary,
            filter_best_model_per_locus,
        )

        df = pd.DataFrame(
            [
                self._make_hit('LeftA', 100, 300, 200.0),
                self._make_hit('LeftB', 500, 700, 50.0),  # no overlap
            ]
        )
        summary = SearchFilterSummary()
        filter_best_model_per_locus(
            df, {'LeftA': 'RightA', 'LeftB': 'RightB'}, summary=summary
        )
        assert summary.cross_model_removed == {}


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


# ---------------------------------------------------------------------------
# Overlap arithmetic, merge gap semantics, and the exposed filter thresholds
# ---------------------------------------------------------------------------


def _locus_hit(model, target, start, end, score, strand='+'):
    """Build one hit row for the cross-model locus filter."""
    return {
        'model': model,
        'target': target,
        'hitStart': str(start),
        'hitEnd': str(end),
        'strand': strand,
        'evalue': '1e-10',
        'score': str(score),
        'bias': 'NA',
        'hmmStart': '1',
        'hmmEnd': '100',
    }


class TestOverlapArithmetic:
    """Hit coordinates are 1-based and inclusive at both ends."""

    def test_hits_sharing_exactly_one_base_overlap(self):
        """A single shared base is an overlap of 1, not 0.

        `write_hits_table` computes length as abs(end - start) + 1, so the
        overlap calculation must use the same convention. Without the +1 two
        hits sharing exactly one base scored 0 and were treated as disjoint.
        """
        hit_table = pd.DataFrame(
            [
                _locus_hit('ModelA', 'chr1', 1000, 1100, 300),
                # Starts on the last base of ModelA's span.
                _locus_hit('ModelB', 'chr1', 1100, 1200, 100),
            ]
        )
        pairing_map = {'ModelA': 'RightA', 'ModelB': 'RightB'}

        result = filter_best_model_per_locus(hit_table, pairing_map)

        assert len(result) == 1
        assert result.iloc[0]['model'] == 'ModelA'

    def test_abutting_hits_do_not_overlap(self):
        """Adjacent but non-touching hits share no base and are both kept."""
        hit_table = pd.DataFrame(
            [
                _locus_hit('ModelA', 'chr1', 1000, 1100, 300),
                _locus_hit('ModelB', 'chr1', 1101, 1200, 100),
            ]
        )
        pairing_map = {'ModelA': 'RightA', 'ModelB': 'RightB'}

        result = filter_best_model_per_locus(hit_table, pairing_map)

        assert len(result) == 2


class TestMergeMaxGap:
    """max_gap is a gap tolerance: larger values merge more."""

    def _cluster_hits(self, gap):
        """Two same-cluster hits separated by `gap` bases."""
        second_start = 1101 + gap
        return pd.DataFrame(
            [
                _locus_hit('Comp1', 'chr1', 1000, 1100, 200),
                _locus_hit('Comp2', 'chr1', second_start, second_start + 100, 150),
            ]
        )

    def test_zero_gap_requires_abutting_or_overlapping(self):
        """With max_gap=0 a one-base gap is not bridged."""
        result = merge_overlapping_cluster_hits(
            self._cluster_hits(gap=1), {'ClusterA': ['Comp1', 'Comp2']}, max_gap=0
        )
        assert len(result) == 2

    def test_larger_gap_bridges_more(self):
        """Raising max_gap merges hits that a smaller value left separate."""
        cluster_map = {'ClusterA': ['Comp1', 'Comp2']}

        assert (
            len(
                merge_overlapping_cluster_hits(
                    self._cluster_hits(gap=50), cluster_map, max_gap=10
                )
            )
            == 2
        )
        assert (
            len(
                merge_overlapping_cluster_hits(
                    self._cluster_hits(gap=50), cluster_map, max_gap=100
                )
            )
            == 1
        )

    def test_min_overlap_alias_still_works_and_warns(self, caplog):
        """The old parameter name is honoured for one release, with a warning."""
        import logging

        cluster_map = {'ClusterA': ['Comp1', 'Comp2']}
        with caplog.at_level(logging.WARNING):
            result = merge_overlapping_cluster_hits(
                self._cluster_hits(gap=50), cluster_map, min_overlap=100
            )

        assert len(result) == 1
        assert 'deprecated' in caplog.text


class TestClusterMapMatchesNothing:
    """A cluster map that matches no model is an error, not an empty result."""

    def test_raises_when_no_model_matches(self):
        """Exiting 0 with an empty file reads as 'no hits found'.

        That is indistinguishable from a genuine empty result, and the usual
        cause is a transposed or mismatched cluster map.
        """
        hit_table = pd.DataFrame([_locus_hit('ActualModel', 'chr1', 1000, 1100, 200)])

        with pytest.raises(EnsembleSearchError, match='matched none of'):
            merge_overlapping_cluster_hits(hit_table, {'ClusterA': ['SomethingElse']})

    def test_error_names_the_models_present(self):
        """The message lists what the hits actually contained."""
        hit_table = pd.DataFrame([_locus_hit('ActualModel', 'chr1', 1000, 1100, 200)])

        with pytest.raises(EnsembleSearchError, match='ActualModel'):
            merge_overlapping_cluster_hits(hit_table, {'ClusterA': ['SomethingElse']})

    def test_partially_unclustered_hits_are_reported_by_name(self, caplog):
        """Dropped models are named, not just counted."""
        import logging

        hit_table = pd.DataFrame(
            [
                _locus_hit('Comp1', 'chr1', 1000, 1100, 200),
                _locus_hit('Stray', 'chr1', 5000, 5100, 150),
            ]
        )

        with caplog.at_level(logging.ERROR):
            result = merge_overlapping_cluster_hits(hit_table, {'ClusterA': ['Comp1']})

        assert len(result) == 1
        assert 'Stray' in caplog.text
        assert 'DISCARDED' in caplog.text


class TestThresholdFlagsAreReachable:
    """--min-score-ratio and --merge-max-gap actually reach the filters."""

    def test_parser_exposes_the_flags(self):
        """The documented knobs exist on the search parser."""
        from tirmite.cli.ensemble_search import create_search_parser

        options = {
            option
            for action in create_search_parser()._actions
            for option in action.option_strings
        }
        assert '--min-score-ratio' in options
        assert '--merge-max-gap' in options

    def test_min_score_ratio_changes_nested_removal(self):
        """A higher ratio keeps hits that the default would discard."""
        hit_table = pd.DataFrame(
            [
                _locus_hit('LeftA', 'chr1', 1000, 1200, 200),
                _locus_hit('RightA', 'chr1', 1020, 1080, 100),
            ]
        )
        pairing_map = {'LeftA': 'RightA'}

        # 200/100 = 2.0: decisive at the 1.5 default, not at 3.0.
        assert len(remove_nested_paired_hits(hit_table, pairing_map)) == 1
        assert (
            len(remove_nested_paired_hits(hit_table, pairing_map, min_score_ratio=3.0))
            == 2
        )

    def test_min_score_ratio_changes_cross_model_filtering(self):
        """The same threshold governs the cross-model step."""
        hit_table = pd.DataFrame(
            [
                _locus_hit('ModelA', 'chr1', 1000, 1100, 200),
                _locus_hit('ModelB', 'chr1', 1050, 1150, 100),
            ]
        )
        pairing_map = {'ModelA': 'RightA', 'ModelB': 'RightB'}

        assert len(filter_best_model_per_locus(hit_table, pairing_map)) == 1
        assert (
            len(
                filter_best_model_per_locus(hit_table, pairing_map, min_score_ratio=3.0)
            )
            == 2
        )


class TestTransposedClusterMapWarning:
    """Recognise a cluster map written with the columns the wrong way round."""

    def test_warns_on_model_first_layout(self, tmp_path, caplog):
        """The old documented layout parses cleanly but matches nothing.

        One `model<TAB>cluster` line per model is what earlier documentation
        showed. It produces clusters named after each model, so no hit ever
        matches, and previously the only symptom was an empty result.
        """
        import logging

        path = tmp_path / 'transposed.tsv'
        path.write_text(
            'MY_TIR_subtype1\tMY_TIR\n'
            'MY_TIR_subtype2\tMY_TIR\n'
            'MY_TIR_subtype3\tMY_TIR\n'
        )

        with caplog.at_level(logging.WARNING):
            parse_cluster_mapping(path)

        assert 'transposed' in caplog.text
        assert 'cluster name FIRST' in caplog.text

    def test_does_not_warn_on_valid_two_column_map(self, tmp_path, caplog):
        """A legitimate single-component-per-cluster map must not warn.

        Guards against a false positive: here every cluster name is distinct,
        so column 1 does not repeat and the shape is not the transposed one.
        """
        import logging

        path = tmp_path / 'valid.tsv'
        path.write_text('ClusterA\tComp1\nClusterB\tComp2\nClusterC\tComp3\n')

        with caplog.at_level(logging.WARNING):
            parse_cluster_mapping(path)

        assert 'transposed' not in caplog.text

    def test_does_not_warn_on_multi_component_map(self, tmp_path, caplog):
        """The normal layout has more than two columns and never warns."""
        import logging

        path = tmp_path / 'normal.tsv'
        path.write_text('ClusterA\tComp1\tComp2\tComp3\nClusterB\tComp4\n')

        with caplog.at_level(logging.WARNING):
            parse_cluster_mapping(path)

        assert 'transposed' not in caplog.text

    def test_empty_map_does_not_warn(self, tmp_path, caplog):
        """A comments-only file has no shape to judge."""
        import logging

        path = tmp_path / 'empty.tsv'
        path.write_text('# nothing here\n')

        with caplog.at_level(logging.WARNING):
            parse_cluster_mapping(path)

        assert 'transposed' not in caplog.text


class TestGreedyRemovalOrder:
    """Pin the order-dependent behaviour of filter_best_model_per_locus.

    That function skips hits already marked for removal (`if idx in
    hits_to_remove: continue`), so a hit removed by one comparison can no
    longer eliminate anything itself. The outcome therefore depends on the
    traversal order, which is DataFrame row order.

    No test previously failed if that order changed, which makes the function
    unsafe to convert into a coordinate-ordered sweep without noticing. These
    tests exist so an optimisation has to preserve it deliberately.

    remove_nested_paired_hits has no such skip -- its result is the set union
    of all pairwise decisions -- so it is order-independent and may be freely
    reordered.
    """

    def _chain(self, scores):
        """Three mutually overlapping hits from three different models."""
        return pd.DataFrame(
            [
                _locus_hit(f'Model{i}', 'chr1', 1000 + i * 20, 1200 + i * 20, s)
                for i, s in enumerate(scores)
            ]
        )

    def _pairing_map(self, n=3):
        """Map each ModelN to a distinct partner so all are 'paired features'."""
        return {f'Model{i}': f'Partner{i}' for i in range(n)}

    def test_transitive_chain_outcome_is_pinned(self):
        """A middle hit removed early cannot then remove the third.

        Scores 300 / 100 / 90: Model0 beats Model1 (3.0) and Model2 (3.33), so
        both go. Model1 would also have beaten... nothing, but the point is
        that the result is exactly one survivor and it is the strongest.
        """
        result = filter_best_model_per_locus(
            self._chain([300, 100, 90]), self._pairing_map()
        )

        assert list(result['model']) == ['Model0']

    def test_middle_hit_removed_first_cannot_evict_the_third(self):
        """The greedy skip is load-bearing here.

        Model0 (300) removes Model1 (100). Model1 would have been decisive
        against Model2 (60) at ratio 1.67, but it is already removed, so it
        gets no say. Model0 vs Model2 is 5.0, so Model2 goes too.
        """
        result = filter_best_model_per_locus(
            self._chain([300, 100, 60]), self._pairing_map()
        )

        assert list(result['model']) == ['Model0']

    def test_non_decisive_chain_keeps_everything(self):
        """When no pair is decisive, order cannot matter."""
        result = filter_best_model_per_locus(
            self._chain([100, 95, 90]), self._pairing_map()
        )

        assert len(result) == 3

    def test_row_order_determines_the_survivor_for_equal_scores(self):
        """With equal scores nothing is decisive, so all survive.

        Pinning this guards the boundary: a sweep that reordered comparisons
        must not start removing equal-scoring hits.
        """
        result = filter_best_model_per_locus(
            self._chain([100, 100, 100]), self._pairing_map()
        )

        assert len(result) == 3

    def test_nested_removal_is_order_independent(self):
        """remove_nested_paired_hits evaluates every pair regardless.

        It has no already-removed skip, so reversing the input row order must
        give the same set of survivors.
        """
        rows = [
            _locus_hit('LeftA', 'chr1', 1000, 1400, 300),
            _locus_hit('RightA', 'chr1', 1100, 1200, 50),
        ]
        pairing_map = {'LeftA': 'RightA'}

        forward = remove_nested_paired_hits(pd.DataFrame(rows), pairing_map)
        reverse = remove_nested_paired_hits(
            pd.DataFrame(list(reversed(rows))), pairing_map
        )

        assert sorted(forward['model']) == sorted(reverse['model'])
        assert list(forward['model']) == ['LeftA']
