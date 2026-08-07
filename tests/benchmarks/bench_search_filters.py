"""Benchmarks for the `tirmite search` hit filters.

``remove_nested_paired_hits`` and ``filter_best_model_per_locus`` are explicit
O(n^2) nested loops over each target group, using ``.loc`` per row. Both were
recently changed to group by target alone rather than by (target, strand),
which roughly doubles the size of each group and therefore quadruples the
pair comparisons -- these benchmarks exist to keep that cost visible.
"""

from conftest import make_hit_table
import pytest

from tirmite.cli.ensemble_search import (
    check_cross_cluster_overlaps,
    filter_best_model_per_locus,
    filter_hits_by_evalue,
    merge_overlapping_cluster_hits,
    remove_nested_paired_hits,
)
from tirmite.core.hit_filters import filter_hits_by_anchor
from tirmite.core.parsers import sort_hit_table

# Sizes chosen from measurement, not guesswork. remove_nested_paired_hits costs
# ~39 ms at 100 hits, ~250 ms at 250 and ~900 ms at 500; at 5000 hits a single
# call takes 93 SECONDS. The quadratic growth is the point of these benchmarks,
# but the top size has to stay runnable in CI, so the sweep stops at 500.
HIT_COUNTS = [100, 250, 500]

# The cheap, vectorised filters scale linearly and can afford a larger table.
LINEAR_HIT_COUNTS = [1000, 5000]

# Clusters matching the model names make_hit_table generates.
CLUSTER_MAP = {'ClusterA': ['Model0', 'Model1'], 'ClusterB': ['Model2', 'Model3']}
PAIRING_MAP = {'Model0': 'Model1', 'Model2': 'Model3'}
MODEL_LENGTHS = {f'Model{i}': 250 for i in range(4)}


@pytest.mark.parametrize('n_hits', HIT_COUNTS)
def test_remove_nested_paired_hits(benchmark, n_hits):
    """Quadratic scan for hits nested inside their pairing partner's."""
    hit_table = make_hit_table(n_hits)

    result = benchmark(remove_nested_paired_hits, hit_table, PAIRING_MAP)

    assert result is not None


@pytest.mark.parametrize('n_hits', HIT_COUNTS)
def test_filter_best_model_per_locus(benchmark, n_hits):
    """Quadratic scan for overlapping hits from competing models."""
    hit_table = make_hit_table(n_hits)

    result = benchmark(filter_best_model_per_locus, hit_table, PAIRING_MAP)

    assert result is not None


@pytest.mark.parametrize('n_hits', LINEAR_HIT_COUNTS)
def test_merge_overlapping_cluster_hits(benchmark, n_hits):
    """Linear sweep merging same-cluster hits after sorting."""
    hit_table = make_hit_table(n_hits)

    result = benchmark(merge_overlapping_cluster_hits, hit_table, CLUSTER_MAP)

    assert result is not None


@pytest.mark.parametrize('n_hits', HIT_COUNTS)
def test_check_cross_cluster_overlaps(benchmark, n_hits):
    """Warning-only quadratic scan across cluster boundaries."""
    hit_table = make_hit_table(n_hits)

    benchmark(check_cross_cluster_overlaps, hit_table, CLUSTER_MAP)


@pytest.mark.parametrize('n_hits', LINEAR_HIT_COUNTS)
def test_filter_hits_by_anchor(benchmark, n_hits):
    """Row-wise anchor filter; iterrows over the whole table."""
    hit_table = make_hit_table(n_hits)

    result = benchmark(
        filter_hits_by_anchor,
        hit_table,
        MODEL_LENGTHS,
        50,
        'F,R',
        PAIRING_MAP,
    )

    assert result is not None


@pytest.mark.parametrize('n_hits', LINEAR_HIT_COUNTS)
def test_filter_hits_by_evalue(benchmark, n_hits):
    """Vectorised threshold filter; the cheap baseline for comparison."""
    hit_table = make_hit_table(n_hits)

    result = benchmark(filter_hits_by_evalue, hit_table, 1e-5)

    assert result is not None


@pytest.mark.parametrize('n_hits', LINEAR_HIT_COUNTS)
def test_sort_hit_table(benchmark, n_hits):
    """Numeric-coordinate sort applied by every importer."""
    hit_table = make_hit_table(n_hits)

    result = benchmark(sort_hit_table, hit_table)

    assert result is not None


def test_full_filter_chain(benchmark, small_hit_table):
    """All three pairing-map filter steps in the order the pipeline runs them.

    Individually fast steps can still add up, and this catches a regression
    that only shows when they compose.
    """

    def run():
        table = filter_hits_by_evalue(small_hit_table, 1e-3)
        table = remove_nested_paired_hits(table, PAIRING_MAP)
        return filter_best_model_per_locus(table, PAIRING_MAP)

    result = benchmark(run)

    assert result is not None
