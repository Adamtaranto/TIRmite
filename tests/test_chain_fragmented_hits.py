"""
Test suite for chain_fragmented_hits function.

This test suite validates that the chaining logic correctly applies all criteria:
1. Hits must belong to the same query
2. Hits are sequential and non-overlapping in query
3. Hits are sequential and non-overlapping in target
4. Hits are all in the same orientation relative to the query
5. Hits are separated by less than max_gap bases in the target sequence
"""

from tirmite.cli.hmm_build import BlastHit, chain_fragmented_hits


class TestChainFragmentedHits:
    """Test cases for chain_fragmented_hits function."""

    def test_empty_hits_list(self):
        """Test that empty list returns empty list."""
        result = chain_fragmented_hits([])
        assert result == []

    def test_single_hit(self):
        """Test that single hit returns single chain."""
        hit = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit])
        assert len(result) == 1
        assert len(result[0]) == 1
        assert result[0][0] == hit

    def test_valid_chain_same_query_forward_strand(self):
        """Test that valid fragments from same query on forward strand are chained."""
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=51,
            query_end=100,
            subject_start=170,  # 20bp gap in target
            subject_end=220,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 1, "Should form a single chain"
        assert len(result[0]) == 2, "Chain should contain both hits"

    def test_valid_chain_same_query_reverse_strand(self):
        """Test that valid fragments from same query on reverse strand are chained."""
        # On reverse strand, subject_start > subject_end
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=220,
            subject_end=170,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=51,
            query_end=100,
            subject_start=150,
            subject_end=100,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 1, "Should form a single chain on reverse strand"
        assert len(result[0]) == 2, "Chain should contain both hits"

    def test_different_query_not_chained(self):
        """Test that hits from different queries are not chained."""
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query2',  # Different query
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=170,
            subject_end=220,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 2, "Should create two separate chains for different queries"

    def test_different_subject_not_chained(self):
        """Test that hits on different subjects/chromosomes are not chained."""
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr2',  # Different chromosome
            query_start=51,
            query_end=100,
            subject_start=170,
            subject_end=220,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 2, "Should create two separate chains for different subjects"

    def test_different_strand_not_chained(self):
        """Test that hits on different strands are not chained."""
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,  # Forward strand
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=51,
            query_end=100,
            subject_start=220,
            subject_end=170,  # Reverse strand
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 2, "Should create two separate chains for different strands"

    def test_overlapping_query_not_chained(self):
        """Test that hits with overlapping query regions are not chained."""
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=40,  # Overlaps with hit1 query (ends at 50)
            query_end=90,
            subject_start=170,
            subject_end=220,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 2, "Should create two separate chains for overlapping query regions"

    def test_overlapping_target_forward_not_chained(self):
        """Test that hits with overlapping target regions on forward strand are not chained."""
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=51,
            query_end=100,
            subject_start=140,  # Overlaps with hit1 target (ends at 150)
            subject_end=190,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 2, "Should create two separate chains for overlapping target regions"

    def test_overlapping_target_reverse_not_chained(self):
        """Test that hits with overlapping target regions on reverse strand are not chained."""
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=220,
            subject_end=170,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=51,
            query_end=100,
            subject_start=180,  # Overlaps with hit1 target (170-220)
            subject_end=130,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 2, "Should create two separate chains for overlapping target regions on reverse strand"

    def test_gap_too_large_not_chained(self):
        """Test that hits with gap larger than max_gap are not chained."""
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=51,
            query_end=100,
            subject_start=700,  # Gap of 550bp (700 - 150 = 550)
            subject_end=750,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 2, "Should create two separate chains when gap exceeds max_gap"

    def test_full_length_hits_not_chained(self):
        """Test that neighboring full-length hits are not chained together.

        This is the specific bug mentioned in the problem statement:
        Two complete, full-length hits from the same query should NOT be chained
        just because they are neighbors in the target sequence.
        """
        # Two full-length hits of the same 100bp query at different locations
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=100,  # Full length query
            subject_start=1000,
            subject_end=1100,
            length=100,
            identity=95.0,
            query_len=100,
            subject_len=10000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,  # Also starts at 1 - not sequential in query!
            query_end=100,  # Full length query
            subject_start=1200,  # Close by in target
            subject_end=1300,
            length=100,
            identity=95.0,
            query_len=100,
            subject_len=10000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 2, "Full-length hits should not be chained (not sequential in query)"
        assert len(result[0]) == 1, "Each chain should contain only one hit"
        assert len(result[1]) == 1, "Each chain should contain only one hit"

    def test_multiple_valid_chains(self):
        """Test that multiple valid chains are correctly formed."""
        # Chain 1: query1, chr1, positions 1-50 and 51-100
        hit1 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=1,
            query_end=50,
            subject_start=100,
            subject_end=150,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit2 = BlastHit(
            query_id='query1',
            subject_id='chr1',
            query_start=51,
            query_end=100,
            subject_start=170,
            subject_end=220,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        # Chain 2: query2, chr2, positions 1-50 and 51-100
        hit3 = BlastHit(
            query_id='query2',
            subject_id='chr2',
            query_start=1,
            query_end=50,
            subject_start=500,
            subject_end=550,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        hit4 = BlastHit(
            query_id='query2',
            subject_id='chr2',
            query_start=51,
            query_end=100,
            subject_start=570,
            subject_end=620,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=1000,
        )
        result = chain_fragmented_hits([hit1, hit2, hit3, hit4], max_gap=500)
        assert len(result) == 2, "Should form two separate chains"
        assert all(len(chain) == 2 for chain in result), "Each chain should have 2 hits"

    def test_example_from_problem_statement(self):
        """Test the specific example given in the problem statement.

        Query is 100bp long, generates 2 non-overlapping partial hits where
        target contains a 20bp insertion in the middle.
        Query coords: hit1: 1-50, hit2: 51-100
        Target coords: hit1: 1-50, hit2: 70-120
        Gap in target is 20bp (70-50=20).
        """
        hit1 = BlastHit(
            query_id='seed_seq',
            subject_id='genome1',
            query_start=1,
            query_end=50,
            subject_start=1,
            subject_end=50,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=5000,
        )
        hit2 = BlastHit(
            query_id='seed_seq',
            subject_id='genome1',
            query_start=51,
            query_end=100,
            subject_start=70,
            subject_end=120,
            length=50,
            identity=95.0,
            query_len=100,
            subject_len=5000,
        )
        result = chain_fragmented_hits([hit1, hit2], max_gap=500)
        assert len(result) == 1, "Should form a single chain"
        assert len(result[0]) == 2, "Chain should contain both hits"

        # Verify the order is correct
        chain = result[0]
        assert chain[0].query_start == 1
        assert chain[0].query_end == 50
        assert chain[1].query_start == 51
        assert chain[1].query_end == 100
