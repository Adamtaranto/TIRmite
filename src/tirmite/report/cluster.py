"""
Average-linkage hierarchical clustering, for ordering a heatmap's axes.

This is UPGMA on a small dense distance matrix, written out rather than taken
from SciPy. The matrix is one row per query model -- tens of entries, not
thousands -- so the naive form is instant, and adding SciPy to the dependency
set for one dendrogram would be a poor trade.

Clustering builds a tree; positions for drawing are then read off the tree in
one pass, so leaf order and merge coordinates cannot disagree.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

__all__ = ['Dendrogram', 'Merge', 'overlap_distances', 'upgma']


@dataclass
class _Leaf:
    """A single label at the bottom of the tree."""

    label: str
    height: float = 0.0


@dataclass
class _Node:
    """A join of two subtrees at a given height."""

    left: Union['_Leaf', '_Node']
    right: Union['_Leaf', '_Node']
    height: float


@dataclass(frozen=True)
class Merge:
    """
    One join in the tree, positioned for drawing.

    Attributes
    ----------
    left_x, right_x : float
        Horizontal positions of the two joined subtrees, in leaf-index units.
    left_height, right_height : float
        Heights the two subtrees were already at.
    height : float
        Height of this join.
    x : float
        Horizontal position of the joined node, midway between its children.
    """

    left_x: float
    right_x: float
    left_height: float
    right_height: float
    height: float
    x: float


@dataclass
class Dendrogram:
    """
    A clustering result ready to plot.

    Attributes
    ----------
    order : list of str
        Leaf labels, left to right.
    merges : list of Merge
        Joins, deepest first, positioned against `order`.
    """

    order: List[str] = field(default_factory=list)
    merges: List[Merge] = field(default_factory=list)

    @property
    def height(self) -> float:
        """
        Return the height of the root join.

        Returns
        -------
        float
            The tallest join, or 0.0 when there is nothing to join.
        """
        return max((merge.height for merge in self.merges), default=0.0)


def overlap_distances(
    models: Sequence[str],
    overlaps: Dict[Tuple[str, str], int],
    hits_per_model: Optional[Dict[str, int]] = None,
) -> List[List[float]]:
    """
    Turn pairwise overlap counts into a distance matrix.

    Parameters
    ----------
    models : sequence of str
        Model names, defining the matrix order.
    overlaps : dict
        Mapping of ``(a, b)`` to the number of loci the two models share, with
        the pair in either order. Self-pairs are ignored.
    hits_per_model : dict, optional
        Total hits per model. When available the distance is a Dice
        dissimilarity, which accounts for how much each model *could* have
        shared; without it, counts are scaled by the largest observed.

    Returns
    -------
    list of list of float
        Symmetric distances in [0, 1], zero on the diagonal.

    Notes
    -----
    Two models that never share a locus sit at distance 1, so they join last.
    Scaling by hit counts matters: a model with three hits that shares all
    three is more similar to its partner than one with three hundred that
    shares three, and a raw count cannot tell those apart.
    """
    size = len(models)
    index = {model: i for i, model in enumerate(models)}
    matrix = [[1.0] * size for _ in range(size)]
    for i in range(size):
        matrix[i][i] = 0.0

    peak = max((count for (a, b), count in overlaps.items() if a != b), default=0)

    for (a, b), count in overlaps.items():
        if a == b or a not in index or b not in index:
            continue
        if hits_per_model:
            total = hits_per_model.get(a, 0) + hits_per_model.get(b, 0)
            similarity = (2 * count / total) if total else 0.0
        else:
            similarity = (count / peak) if peak else 0.0
        distance = min(1.0, max(0.0, 1.0 - similarity))
        i, j = index[a], index[b]
        matrix[i][j] = distance
        matrix[j][i] = distance

    return matrix


def upgma(labels: Sequence[str], distances: List[List[float]]) -> Dendrogram:
    """
    Cluster labels by average linkage and return a drawable tree.

    Parameters
    ----------
    labels : sequence of str
        Leaf labels.
    distances : list of list of float
        Symmetric distance matrix, in the same order as `labels`.

    Returns
    -------
    Dendrogram
        Leaf order and positioned merges. With fewer than two labels the tree
        is the labels themselves and no merges.

    Notes
    -----
    Ties are broken on label order, so the same input always produces the same
    tree. Without that, dictionary iteration order would let two runs over the
    same data draw different figures.
    """
    if len(labels) < 2:
        return Dendrogram(order=list(labels))

    trees: Dict[int, Union[_Leaf, _Node]] = {
        i: _Leaf(label) for i, label in enumerate(labels)
    }
    sizes: Dict[int, int] = dict.fromkeys(range(len(labels)), 1)
    first_label: Dict[int, str] = dict(enumerate(labels))

    dist: Dict[Tuple[int, int], float] = {}
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            dist[(i, j)] = distances[i][j]

    next_id = len(labels)

    while len(trees) > 1:
        best = min(
            (
                (value, first_label[a], first_label[b], a, b)
                for (a, b), value in dist.items()
                if a in trees and b in trees
            ),
            default=None,
        )
        if best is None:
            break
        value, _la, _lb, a, b = best

        height = value / 2.0
        # Average linkage over an inconsistent matrix can propose a join no
        # taller than one of its children, which would draw a tree doubling
        # back on itself.
        height = max(height, trees[a].height, trees[b].height)
        node = _Node(left=trees[a], right=trees[b], height=height)

        for other in list(trees):
            if other in (a, b):
                continue
            d_a = dist.get(_key(a, other))
            d_b = dist.get(_key(b, other))
            if d_a is None or d_b is None:
                continue
            dist[_key(next_id, other)] = (d_a * sizes[a] + d_b * sizes[b]) / (
                sizes[a] + sizes[b]
            )

        for other in list(trees):
            dist.pop(_key(a, other), None)
            dist.pop(_key(b, other), None)

        sizes[next_id] = sizes[a] + sizes[b]
        first_label[next_id] = min(first_label[a], first_label[b])
        del trees[a]
        del trees[b]
        trees[next_id] = node
        next_id += 1

    root = next(iter(trees.values()))
    order: List[str] = []
    _collect_leaves(root, order)
    positions = {label: float(i) for i, label in enumerate(order)}

    merges: List[Merge] = []
    _position(root, positions, merges)
    return Dendrogram(order=order, merges=merges)


def _key(a: int, b: int) -> Tuple[int, int]:
    """
    Return an order-independent key for a pair of cluster ids.

    Parameters
    ----------
    a, b : int
        Cluster identifiers.

    Returns
    -------
    tuple of (int, int)
        The pair, smaller first.
    """
    return (a, b) if a < b else (b, a)


def _collect_leaves(node: Union[_Leaf, _Node], out: List[str]) -> None:
    """
    Append a subtree's leaf labels in left-to-right order.

    Parameters
    ----------
    node : _Leaf or _Node
        Subtree root.
    out : list of str
        Accumulator, extended in place.

    Returns
    -------
    None
        Updates `out`.
    """
    if isinstance(node, _Leaf):
        out.append(node.label)
        return
    _collect_leaves(node.left, out)
    _collect_leaves(node.right, out)


def _position(
    node: Union[_Leaf, _Node],
    positions: Dict[str, float],
    merges: List[Merge],
) -> Tuple[float, float]:
    """
    Compute a subtree's x position, recording a merge for each join.

    Parameters
    ----------
    node : _Leaf or _Node
        Subtree root.
    positions : dict
        Leaf label to x position.
    merges : list of Merge
        Accumulator, appended to in post-order so children precede parents.

    Returns
    -------
    tuple of (float, float)
        The subtree's ``(x, height)``.
    """
    if isinstance(node, _Leaf):
        return positions[node.label], node.height

    left_x, left_height = _position(node.left, positions, merges)
    right_x, right_height = _position(node.right, positions, merges)
    if left_x > right_x:
        left_x, right_x = right_x, left_x
        left_height, right_height = right_height, left_height

    x = (left_x + right_x) / 2.0
    merges.append(
        Merge(
            left_x=left_x,
            right_x=right_x,
            left_height=left_height,
            right_height=right_height,
            height=node.height,
            x=x,
        )
    )
    return x, node.height
