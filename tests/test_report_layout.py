"""Tests for tirmite.report.layout — annotation track stacking."""

import random

import pytest

from tirmite.report.layout import assign_rows, min_gap_for_contig, stack_contig


class FakeHit:
    """Minimal stand-in for a HitRecord: stacking only reads three fields."""

    def __init__(self, uid, start, end):
        self.uid = uid
        self.start = start
        self.end = end


class TestAssignRows:
    def test_empty(self):
        rows, overflow = assign_rows([])
        assert rows == {}
        assert overflow == set()

    def test_disjoint_intervals_share_one_row(self):
        intervals = [(1, 10, 'a'), (20, 30, 'b'), (40, 50, 'c')]
        rows, overflow = assign_rows(intervals)
        assert rows == {'a': 0, 'b': 0, 'c': 0}
        assert overflow == set()

    def test_mutually_overlapping_intervals_get_distinct_rows(self):
        intervals = [(1, 100, 'a'), (2, 100, 'b'), (3, 100, 'c')]
        rows, _ = assign_rows(intervals)
        assert sorted(rows.values()) == [0, 1, 2]

    def test_row_is_reused_once_free(self):
        # b overlaps a, but c starts after a ends and so reclaims row 0.
        intervals = [(1, 10, 'a'), (5, 15, 'b'), (11, 20, 'c')]
        rows, _ = assign_rows(intervals)
        assert rows['a'] == 0
        assert rows['b'] == 1
        assert rows['c'] == 0

    def test_min_gap_keeps_abutting_intervals_apart(self):
        intervals = [(1, 10, 'a'), (12, 20, 'b')]
        assert assign_rows(intervals)[0] == {'a': 0, 'b': 0}
        # With a gap requirement of 5 they are visually touching.
        rows, _ = assign_rows(intervals, min_gap=5)
        assert rows['a'] != rows['b']

    def test_deterministic_under_input_shuffling(self):
        intervals = [
            (1, 50, 'a'),
            (10, 60, 'b'),
            (20, 30, 'c'),
            (70, 80, 'd'),
            (75, 90, 'e'),
        ]
        reference, _ = assign_rows(intervals)
        rng = random.Random(0)
        for _ in range(20):
            shuffled = intervals[:]
            rng.shuffle(shuffled)
            assert assign_rows(shuffled)[0] == reference

    def test_max_rows_caps_height_and_flags_overflow(self):
        intervals = [(1, 100, str(i)) for i in range(6)]
        rows, overflow = assign_rows(intervals, max_rows=3)
        assert max(rows.values()) == 2
        assert len(overflow) == 3
        # Everything that overflowed is parked on the last row.
        assert all(rows[key] == 2 for key in overflow)

    def test_max_rows_not_reached_gives_no_overflow(self):
        intervals = [(1, 100, 'a'), (2, 100, 'b')]
        _, overflow = assign_rows(intervals, max_rows=10)
        assert overflow == set()

    def test_row_count_equals_max_overlap_depth(self):
        # Greedy colouring of intervals is optimal, so the row count is
        # exactly the deepest pile-up.
        intervals = [(1, 100, 'a'), (10, 20, 'b'), (30, 40, 'c'), (35, 45, 'd')]
        rows, _ = assign_rows(intervals)
        assert max(rows.values()) + 1 == 3


class TestMinGapForContig:
    @pytest.mark.parametrize(
        'length,expected',
        [(0, 1), (1, 1), (2000, 1), (200_000, 100), (30_000_000, 15000)],
    )
    def test_scales_with_contig_length(self, length, expected):
        assert min_gap_for_contig(length) == expected


class TestStackContig:
    def test_pair_members_share_a_row(self):
        hits = [FakeHit(0, 100, 200), FakeHit(1, 5000, 5100)]
        rows, _ = stack_contig(hits, [(0, 1)], contig_length=100_000)
        assert rows[0] == rows[1]

    def test_hit_inside_an_element_gets_its_own_row(self):
        # An unpaired hit within the element's span must not be drawn over the
        # element's connecting line.
        hits = [
            FakeHit(0, 100, 200),
            FakeHit(1, 5000, 5100),
            FakeHit(2, 2000, 2100),
        ]
        rows, _ = stack_contig(hits, [(0, 1)], contig_length=100_000)
        assert rows[0] == rows[1]
        assert rows[2] != rows[0]

    def test_hit_outside_an_element_reuses_the_row(self):
        hits = [
            FakeHit(0, 100, 200),
            FakeHit(1, 5000, 5100),
            FakeHit(2, 50_000, 50_100),
        ]
        rows, _ = stack_contig(hits, [(0, 1)], contig_length=100_000)
        assert rows[2] == rows[0]

    def test_nested_elements_stack(self):
        hits = [
            FakeHit(0, 100, 200),
            FakeHit(1, 9000, 9100),
            FakeHit(2, 1000, 1100),
            FakeHit(3, 8000, 8100),
        ]
        rows, _ = stack_contig(hits, [(0, 1), (2, 3)], contig_length=100_000)
        assert rows[0] == rows[1]
        assert rows[2] == rows[3]
        assert rows[0] != rows[2]

    def test_pair_with_a_missing_member_falls_back_to_single_stacking(self):
        hits = [FakeHit(0, 100, 200)]
        rows, overflow = stack_contig(hits, [(0, 99)], contig_length=100_000)
        assert rows == {0: 0}
        assert overflow == set()

    def test_no_hits(self):
        assert stack_contig([], [], contig_length=1000) == ({}, set())

    def test_max_rows_flags_both_pair_members(self):
        hits = []
        pairs = []
        for i in range(4):
            hits.append(FakeHit(2 * i, 100, 200))
            hits.append(FakeHit(2 * i + 1, 5000, 5100))
            pairs.append((2 * i, 2 * i + 1))
        rows, overflow = stack_contig(hits, pairs, contig_length=100_000, max_rows=2)
        assert max(rows.values()) == 1
        # Overflow is reported for whole pairs, never half of one.
        for left, right in pairs:
            assert (left in overflow) == (right in overflow)

    def test_deterministic_under_input_shuffling(self):
        hits = [FakeHit(i, i * 1000 + 1, i * 1000 + 900) for i in range(10)]
        pairs = [(0, 3), (4, 7)]
        reference, _ = stack_contig(hits, pairs, contig_length=100_000)
        rng = random.Random(1)
        for _ in range(10):
            shuffled = hits[:]
            rng.shuffle(shuffled)
            assert stack_contig(shuffled, pairs, contig_length=100_000)[0] == reference
