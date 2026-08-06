"""Benchmarks for HTML report construction.

Report building runs once at the end of a pairing run, but it touches every hit
the run retained, so its cost scales with the whole result set rather than with
anything the user chose. Two stages are worth watching:

1. **Row stacking** (``stack_contig``) is the only super-linear-looking step:
   it sorts and sweeps every interval on a contig. It is O(n log n), and these
   benchmarks exist so a future change that makes it quadratic -- an inner scan
   over open rows, say -- shows up immediately rather than as a report that
   quietly takes minutes on a real genome.

2. **Data assembly** (``PairReportAccumulator``) walks the hit table once per
   hit and then again per pairing group. It is the larger absolute cost by an
   order of magnitude, and it is where added per-hit work would land.

Measured at 2k, 20k and 100k hits (all on one contig, every element
overlapping its neighbour, which is the worst case for stacking)::

    stack_contig        0.002 s    0.021 s    0.166 s
    build_report_data   0.066 s    0.617 s    3.52 s
    to_json             0.013 s    0.253 s    0.666 s

Both scale about linearly, so nothing here needs attention at realistic sizes:
a genome producing 100k terminus hits costs about four seconds of report
building on top of the pairing run that found them.

The JSON payload is the thing to watch. 100k hits serialise to 18 MB, and more
than half of that is the ``elements`` array, which is still a list of objects
while the hits are columnar. If report size ever becomes the complaint,
columnarising the elements is the change to make -- not shrinking the hit
columns.

Rendering and figure drawing are deliberately not benchmarked here: they pull
in jinja2 and matplotlib, whose own import and draw times would dominate and
drift with those dependencies rather than with TIRmite.

Setup is excluded from every measurement via ``benchmark.pedantic(setup=...)``.
The accumulator is stateful -- ``add_group`` mutates it and ``finalise``
consumes it -- so each round needs a fresh one, and building the hit table is
expensive enough to swamp the result if it were timed too.
"""

from conftest import make_pairing_hit_rows, rows_to_hit_table
import pytest

from tirmite.core.pairing import PairingConfig, table2dict
from tirmite.report.collect import PairReportAccumulator
from tirmite.report.layout import stack_contig

# Each element contributes two hits, so these are 2k, 20k and 100k hits.
ELEMENT_COUNTS = [1_000, 10_000, 50_000]

MODEL_LENGTHS = {'LeftTIR': 100, 'RightTIR': 100}


class FakeHit:
    """The three fields stacking reads, without building a full HitRecord."""

    __slots__ = ('uid', 'start', 'end')

    def __init__(self, uid, start, end):
        self.uid = uid
        self.start = start
        self.end = end


def make_stacking_input(n_elements: int):
    """
    Build hits and pair memberships for one contig.

    Parameters
    ----------
    n_elements : int
        Number of elements; each contributes two hits and one pair.

    Returns
    -------
    tuple
        ``(hits, pairs, contig_length)`` ready for :func:`stack_contig`.

    Notes
    -----
    Every element overlaps its neighbour, which is the expensive case: rows
    are never freed early, so the sweep carries the maximum number of open
    rows throughout.
    """
    hits = []
    pairs = []
    for i in range(n_elements):
        base = 1000 + i * 500
        left_uid = 2 * i
        right_uid = 2 * i + 1
        hits.append(FakeHit(left_uid, base, base + 99))
        hits.append(FakeHit(right_uid, base + 2000, base + 2099))
        pairs.append((left_uid, right_uid))
    return hits, pairs, 1000 + n_elements * 500 + 3000


def make_accumulator_input(n_elements: int):
    """
    Build the inputs one pairing group's worth of report data needs.

    Parameters
    ----------
    n_elements : int
        Number of elements to describe.

    Returns
    -------
    tuple
        ``(hit_table, hit_index, paired, config)``.
    """
    rows = make_pairing_hit_rows(n_elements)
    hit_table = rows_to_hit_table(rows)
    _, hit_index = table2dict(hit_table)
    paired = {'LeftTIR': [{2 * i, 2 * i + 1} for i in range(n_elements)]}
    config = PairingConfig(
        orientation='F,R', left_model='LeftTIR', right_model='RightTIR'
    )
    return hit_table, hit_index, paired, config


@pytest.mark.parametrize('n_elements', ELEMENT_COUNTS)
def test_stack_contig(benchmark, n_elements):
    """Measure row assignment for one densely populated contig."""
    state = {}

    def setup():
        state['args'] = make_stacking_input(n_elements)
        return (), {}

    def run():
        hits, pairs, length = state['args']
        return stack_contig(hits, pairs, length, max_rows=30)

    benchmark.pedantic(run, setup=setup, rounds=3, iterations=1)


@pytest.mark.parametrize('n_elements', ELEMENT_COUNTS)
def test_build_report_data(benchmark, n_elements):
    """Measure assembling a full report's data from a pairing result."""
    state = {}

    def setup():
        hit_table, hit_index, paired, config = make_accumulator_input(n_elements)
        accumulator = PairReportAccumulator(
            hit_table=hit_table,
            model_lengths=MODEL_LENGTHS,
            # A constant length keeps sequence-source I/O out of the timing;
            # this measures the assembly, not pyfaidx.
            contig_length=lambda name: 100_000_000,
        )
        state['acc'] = accumulator
        state['args'] = (hit_index, paired, config)
        return (), {}

    def run():
        hit_index, paired, config = state['args']
        accumulator = state['acc']
        accumulator.add_group(
            left_feature='LeftTIR',
            right_feature='RightTIR',
            config=config,
            paired=paired,
            hit_index=hit_index,
        )
        # Alignment panels are off: they read sequence and may shell out to
        # MAFFT, neither of which belongs in a data-assembly measurement.
        return accumulator.finalise(embed_sequences=False, msa_mode='off')

    benchmark.pedantic(run, setup=setup, rounds=3, iterations=1)


@pytest.mark.parametrize('n_elements', [10_000])
def test_serialise_report_data(benchmark, n_elements):
    """Measure JSON serialisation of an assembled report."""
    state = {}

    def setup():
        hit_table, hit_index, paired, config = make_accumulator_input(n_elements)
        accumulator = PairReportAccumulator(
            hit_table=hit_table,
            model_lengths=MODEL_LENGTHS,
            contig_length=lambda name: 100_000_000,
        )
        accumulator.add_group(
            left_feature='LeftTIR',
            right_feature='RightTIR',
            config=config,
            paired=paired,
            hit_index=hit_index,
        )
        state['data'] = accumulator.finalise(embed_sequences=False, msa_mode='off')
        return (), {}

    def run():
        return state['data'].to_json()

    benchmark.pedantic(run, setup=setup, rounds=3, iterations=1)
