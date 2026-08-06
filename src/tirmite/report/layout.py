"""
Interval stacking for the report's annotation tracks.

Annotations that overlap must not obscure one another, so each is assigned a
row within its contig's track. This is the classic interval-graph colouring
problem, solved greedily with a sweep: the greedy solution over intervals is
optimal, so the number of rows produced is exactly the maximum depth of overlap.

The stacking runs in Python at report-build time rather than in the browser
because the assignment must be **stable across zoom**. Recomputing rows for the
visible window would make glyphs jump vertically while the user scrolls, which
defeats the point of stacking them at all. Doing it here also means the track
height is known when the HTML is written, and it is unit-testable with plain
pytest -- the repo has no JavaScript test harness.
"""

import heapq
import logging
from typing import (
    Any,
    Dict,
    Hashable,
    List,
    Optional,
    Protocol,
    Sequence,
    Set,
    Tuple,
    cast,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['Stackable', 'assign_rows', 'min_gap_for_contig', 'stack_contig']


class Stackable(Protocol):
    """
    The three fields :func:`stack_contig` needs from a hit.

    Structural rather than nominal so that stacking stays decoupled from
    :class:`tirmite.report.model.HitRecord` and remains testable with a stub.
    """

    # Declared read-only so that frozen dataclasses satisfy the protocol.
    @property
    def uid(self) -> int:
        """Stable identifier for the hit."""

    @property
    def start(self) -> int:
        """Inclusive lower genomic coordinate."""

    @property
    def end(self) -> int:
        """Inclusive higher genomic coordinate."""


# Roughly one screen pixel at whole-contig zoom. Two annotations closer than
# this are visually touching even though their coordinates do not overlap, so
# they are stacked apart.
_PIXELS_ACROSS_TRACK = 2000


def min_gap_for_contig(contig_length: int) -> int:
    """
    Return the separation below which two annotations visually collide.

    Parameters
    ----------
    contig_length : int
        Length of the contig in base pairs.

    Returns
    -------
    int
        Minimum gap in base pairs. At least 1, so that abutting features are
        always separated.
    """
    return max(1, contig_length // _PIXELS_ACROSS_TRACK)


def assign_rows(
    intervals: Sequence[Tuple[int, int, Hashable]],
    *,
    min_gap: int = 0,
    max_rows: Optional[int] = None,
) -> Tuple[Dict[Hashable, int], Set[Hashable]]:
    """
    Assign each interval to a row so that no two in a row overlap.

    Parameters
    ----------
    intervals : sequence of (int, int, hashable)
        ``(start, end, key)`` triples, inclusive coordinates. Keys must be
        unique; a repeated key overwrites its earlier assignment.
    min_gap : int, default 0
        Additional separation required between two intervals sharing a row.
        Use :func:`min_gap_for_contig` to keep visually touching features
        apart.
    max_rows : int, optional
        Cap on the number of rows. Intervals that would need a row beyond the
        cap are placed on the last row and reported as overflow.

    Returns
    -------
    rows : dict
        Mapping of interval key to zero-based row index.
    overflow : set
        Keys that exceeded `max_rows` and were forced onto the last row. Empty
        when `max_rows` is None.

    Notes
    -----
    Intervals are sorted by ``(start, end, key)`` and the lowest-indexed free
    row is always reused, so the result depends only on the intervals
    themselves and not on the order they were supplied in. Tests rely on that
    determinism, and so does a stable rendering across runs.
    """
    rows: Dict[Hashable, int] = {}
    overflow: Set[Hashable] = set()
    if not intervals:
        return rows, overflow

    # Sorting by key as the final tiebreak makes the output independent of
    # input order. Keys are heterogeneous in principle, so compare their reprs
    # rather than the keys themselves.
    ordered = sorted(intervals, key=lambda iv: (iv[0], iv[1], repr(iv[2])))

    # Rows currently occupied, as a heap of (row_end, row_index). Freed rows
    # move to `free`, a heap of row indices, so the lowest free row is reused.
    occupied: List[Tuple[int, int]] = []
    free: List[int] = []
    n_rows = 0

    for start, end, key in ordered:
        # Release every row whose occupant finished before this interval needs
        # to begin. The heap is ordered by row_end, so one pass suffices.
        while occupied and occupied[0][0] < start - min_gap:
            _, freed_row = heapq.heappop(occupied)
            heapq.heappush(free, freed_row)

        if free:
            row = heapq.heappop(free)
        elif max_rows is not None and n_rows >= max_rows:
            # Out of rows: park the interval on the last one. It will overlap
            # something, which the UI flags rather than hides.
            row = max_rows - 1
            overflow.add(key)
            rows[key] = row
            continue
        else:
            row = n_rows
            n_rows += 1

        rows[key] = row
        heapq.heappush(occupied, (end, row))

    return rows, overflow


def stack_contig(
    hits: Sequence[Stackable],
    pairs: Sequence[Tuple[int, int]],
    contig_length: int,
    *,
    max_rows: int = 30,
) -> Tuple[Dict[int, int], Set[int]]:
    """
    Assign track rows to every hit on one contig.

    Parameters
    ----------
    hits : sequence
        Hit objects carrying ``uid``, ``start`` and ``end`` attributes. All
        must lie on the same contig.
    pairs : sequence of (int, int)
        ``(left_uid, right_uid)`` for each terminus pair on this contig.
    contig_length : int
        Contig length in base pairs, used to derive the visual collision gap.
    max_rows : int, default 30
        Cap on track height.

    Returns
    -------
    rows : dict
        Mapping of hit uid to row index.
    overflow : set
        Uids forced onto the last row because the cap was reached.

    Notes
    -----
    A pair is stacked as a **single composite interval** spanning both of its
    termini, and both members inherit that row. That is what keeps the
    connecting line between a pair short and horizontal; stacking the termini
    independently would routinely split a pair across rows and turn every link
    into a long diagonal. Unpaired hits are stacked as their own intervals in
    the same pass, so they fill whatever rows are free -- including inside the
    span of an element, where they correctly land on a different row rather
    than being drawn over the element's link line.

    Rows are assigned across all pairing groups at once, so a contig gets one
    coherent track rather than one per group.
    """
    by_uid = {h.uid: h for h in hits}
    min_gap = min_gap_for_contig(contig_length)

    intervals: List[Tuple[int, int, Hashable]] = []
    paired_uids: Set[int] = set()

    for left_uid, right_uid in pairs:
        left = by_uid.get(left_uid)
        right = by_uid.get(right_uid)
        if left is None or right is None:
            # A pair whose members are not both present on this contig cannot
            # be drawn as one interval; its termini stack individually below.
            continue
        intervals.append(
            (
                min(left.start, right.start),
                max(left.end, right.end),
                ('pair', left_uid, right_uid),
            )
        )
        paired_uids.add(left_uid)
        paired_uids.add(right_uid)

    for hit in hits:
        if hit.uid not in paired_uids:
            intervals.append((hit.start, hit.end, ('hit', hit.uid)))

    assigned, overflow_keys = assign_rows(intervals, min_gap=min_gap, max_rows=max_rows)

    rows: Dict[int, int] = {}
    overflow: Set[int] = set()
    for key, row in assigned.items():
        # Keys are the tuples built above: ('pair', left, right) or ('hit', uid).
        parts = cast(Tuple[Any, ...], key)
        uids = parts[1:]
        for uid in uids:
            rows[uid] = row
            if key in overflow_keys:
                overflow.add(uid)

    return rows, overflow
