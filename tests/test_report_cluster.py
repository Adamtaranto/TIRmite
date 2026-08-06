"""Tests for tirmite.report.cluster — UPGMA for ordering the heatmap axes."""

import random

import pytest

from tirmite.report.cluster import overlap_distances, upgma


class TestOverlapDistances:
    def test_no_overlap_is_maximum_distance(self):
        matrix = overlap_distances(['a', 'b'], {}, {'a': 10, 'b': 10})
        assert matrix == [[0.0, 1.0], [1.0, 0.0]]

    def test_complete_sharing_is_zero_distance(self):
        # Both models hit 10 loci and share all 10.
        matrix = overlap_distances(['a', 'b'], {('a', 'b'): 10}, {'a': 10, 'b': 10})
        assert matrix[0][1] == pytest.approx(0.0)

    def test_distance_is_symmetric_with_a_zero_diagonal(self):
        matrix = overlap_distances(
            ['a', 'b', 'c'], {('a', 'b'): 3}, {'a': 10, 'b': 10, 'c': 10}
        )
        assert matrix[0][1] == matrix[1][0]
        assert all(matrix[i][i] == 0.0 for i in range(3))

    def test_hit_counts_scale_the_similarity(self):
        # The same raw count means much more for a model with few hits.
        few = overlap_distances(['a', 'b'], {('a', 'b'): 3}, {'a': 4, 'b': 4})
        many = overlap_distances(['a', 'b'], {('a', 'b'): 3}, {'a': 400, 'b': 400})
        assert few[0][1] < many[0][1]

    def test_without_hit_counts_the_largest_count_sets_the_scale(self):
        matrix = overlap_distances(['a', 'b', 'c'], {('a', 'b'): 10, ('a', 'c'): 5})
        assert matrix[0][1] == pytest.approx(0.0)
        assert matrix[0][2] == pytest.approx(0.5)

    def test_self_pairs_are_ignored(self):
        matrix = overlap_distances(['a', 'b'], {('a', 'a'): 99}, {'a': 99, 'b': 1})
        assert matrix[0][0] == 0.0
        assert matrix[0][1] == 1.0

    def test_unknown_models_are_skipped(self):
        matrix = overlap_distances(['a'], {('a', 'ghost'): 5}, {'a': 5})
        assert matrix == [[0.0]]


class TestUpgma:
    def test_no_labels(self):
        tree = upgma([], [])
        assert tree.order == []
        assert tree.merges == []

    def test_single_label_has_no_merges(self):
        tree = upgma(['only'], [[0.0]])
        assert tree.order == ['only']
        assert tree.merges == []
        assert tree.height == 0.0

    def test_two_labels_make_one_merge(self):
        tree = upgma(['a', 'b'], [[0.0, 0.4], [0.4, 0.0]])
        assert sorted(tree.order) == ['a', 'b']
        assert len(tree.merges) == 1
        assert tree.merges[0].height == pytest.approx(0.2)

    def test_closest_pair_joins_first(self):
        labels = ['a', 'b', 'c']
        # a and c are closest.
        d = [
            [0.0, 0.9, 0.1],
            [0.9, 0.0, 0.8],
            [0.1, 0.8, 0.0],
        ]
        tree = upgma(labels, d)
        assert tree.merges[0].height == pytest.approx(0.05)
        assert abs(tree.order.index('a') - tree.order.index('c')) == 1

    def test_groups_stay_contiguous_in_the_leaf_order(self):
        labels = ['a1', 'a2', 'a3', 'b1', 'b2']
        overlaps = {
            ('a1', 'a2'): 20,
            ('a1', 'a3'): 18,
            ('a2', 'a3'): 19,
            ('b1', 'b2'): 22,
            ('a1', 'b1'): 1,
            ('a3', 'b2'): 1,
        }
        hits = dict.fromkeys(labels, 25)
        tree = upgma(labels, overlap_distances(labels, overlaps, hits))

        a_positions = sorted(tree.order.index(m) for m in ('a1', 'a2', 'a3'))
        b_positions = sorted(tree.order.index(m) for m in ('b1', 'b2'))
        # Each group occupies a run of consecutive positions.
        assert a_positions == list(range(a_positions[0], a_positions[0] + 3))
        assert b_positions == list(range(b_positions[0], b_positions[0] + 2))

    def test_merge_heights_never_decrease(self):
        labels = ['a', 'b', 'c', 'd']
        d = [
            [0.0, 0.1, 0.8, 0.9],
            [0.1, 0.0, 0.7, 0.9],
            [0.8, 0.7, 0.0, 0.2],
            [0.9, 0.9, 0.2, 0.0],
        ]
        tree = upgma(labels, d)
        for merge in tree.merges:
            # A join must sit at or above both the subtrees it connects, or
            # the drawn tree doubles back on itself.
            assert merge.height >= merge.left_height
            assert merge.height >= merge.right_height

    def test_merge_positions_lie_between_their_children(self):
        labels = ['a', 'b', 'c', 'd']
        d = [
            [0.0, 0.1, 0.8, 0.9],
            [0.1, 0.0, 0.7, 0.9],
            [0.8, 0.7, 0.0, 0.2],
            [0.9, 0.9, 0.2, 0.0],
        ]
        tree = upgma(labels, d)
        for merge in tree.merges:
            assert merge.left_x <= merge.x <= merge.right_x
            # Positions are leaf indices, so they must lie on the axis.
            assert 0 <= merge.left_x <= len(labels) - 1
            assert 0 <= merge.right_x <= len(labels) - 1

    def test_every_label_appears_exactly_once(self):
        labels = [f'm{i}' for i in range(9)]
        rng = random.Random(3)
        d = [[0.0] * 9 for _ in range(9)]
        for i in range(9):
            for j in range(i + 1, 9):
                value = rng.random()
                d[i][j] = d[j][i] = value
        tree = upgma(labels, d)
        assert sorted(tree.order) == sorted(labels)
        assert len(tree.merges) == len(labels) - 1

    def test_result_is_deterministic(self):
        labels = ['a', 'b', 'c', 'd']
        d = [
            [0.0, 0.5, 0.5, 0.9],
            [0.5, 0.0, 0.5, 0.9],
            [0.5, 0.5, 0.0, 0.9],
            [0.9, 0.9, 0.9, 0.0],
        ]
        # Ties everywhere: without a deterministic tiebreak, dict iteration
        # order would let two runs draw different trees.
        first = upgma(labels, d)
        for _ in range(10):
            assert upgma(labels, d).order == first.order

    def test_identical_models_join_at_zero(self):
        labels = ['a', 'b', 'c']
        d = [
            [0.0, 0.0, 0.9],
            [0.0, 0.0, 0.9],
            [0.9, 0.9, 0.0],
        ]
        tree = upgma(labels, d)
        assert tree.merges[0].height == pytest.approx(0.0)
        assert tree.height == pytest.approx(0.45)
