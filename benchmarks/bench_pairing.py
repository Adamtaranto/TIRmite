"""Benchmarks for the terminus pairing engine.

Two things measurement established here, both contrary to the obvious guess:

1. **Candidate discovery dominates, not iteration.** ``parseHitsGeneral``
   (which calls ``_find_candidates``) costs ~0.6 s for 250 asymmetric
   elements, against ~0.017 s for the ``iterateGetPairsAsymmetric`` pass that
   follows. The iterate functions loop until stability and look expensive;
   they are not. Optimisation effort for the asymmetric path belongs in
   ``_find_candidates``.

2. **The symmetric path is orders of magnitude slower than the asymmetric
   one**, and getting worse with size:

   =========  =========================  =============================
   Elements   ``iterateGetPairsCustom``  ``iterateGetPairsAsymmetric``
   =========  =========================  =============================
   50         0.96 s                     0.0016 s
   100        12.4 s                     0.0040 s
   250        427 s                      0.0168 s
   =========  =========================  =============================

   That is a 600x gap at 50 elements widening to 25,000x at 250, and the
   symmetric growth is far worse than quadratic. ``F,F`` and ``R,R`` are the
   documented LTR use case, so this is the configuration users are pointed at
   for LTR elements. The symmetric sweep is therefore capped well below the
   asymmetric one to keep this suite runnable.

Setup is excluded from every measurement via ``benchmark.pedantic(setup=...)``:
the ``iterateGetPairs*`` functions mutate the hit index they are given, so each
round needs a fresh one, and building it is expensive enough to swamp the
result if it were timed too.
"""

from conftest import make_pairing_hit_rows, rows_to_hit_table
import pytest

from tirmite.core.pairing import (
    PairingConfig,
    getPairsAsymmetric,
    getPairsSymmetric,
    iterateGetPairsAsymmetric,
    iterateGetPairsCustom,
    parseHitsGeneral,
    table2dict,
)

# Asymmetric sweep. parseHitsGeneral is the cost here: ~0.02 s at 50 elements,
# ~0.08 s at 100, ~0.6 s at 250.
ELEMENT_COUNTS = [50, 100, 250]

# Symmetric sweep, deliberately smaller. iterateGetPairsCustom takes ~12 s at
# 100 elements and ~427 s at 250, so anything above 50 is unusable in CI.
SYMMETRIC_ELEMENT_COUNTS = [25, 50]


def _asymmetric_inputs(n_elements, orientation='F,R'):
    """
    Build a fresh hit index and config for an asymmetric pairing run.

    Parameters
    ----------
    n_elements : int
        Number of synthetic elements to generate.
    orientation : str, default 'F,R'
        Orientation string for the pairing configuration.

    Returns
    -------
    tuple
        ``(hitIndex, config)``, with the index freshly built so the caller may
        mutate it.
    """
    rows = make_pairing_hit_rows(n_elements)
    hits_dict, hit_index = table2dict(rows_to_hit_table(rows))
    config = PairingConfig(
        orientation=orientation, left_model='LeftTIR', right_model='RightTIR'
    )
    hit_index = parseHitsGeneral(hitsDict=hits_dict, hitIndex=hit_index, config=config)
    return hit_index, config


def _symmetric_inputs(n_elements, orientation='F,F'):
    """
    Build a fresh hit index and config for a symmetric pairing run.

    Parameters
    ----------
    n_elements : int
        Number of synthetic elements to generate.
    orientation : str, default 'F,F'
        Orientation string for the pairing configuration.

    Returns
    -------
    tuple
        ``(hitIndex, config)``, with the index freshly built so the caller may
        mutate it.
    """
    rows = make_pairing_hit_rows(n_elements, left_model='LTR', right_model='LTR')
    # Same-strand symmetric: one model supplies both termini, both on '+'.
    for row in rows:
        row['strand'] = '+'
    hits_dict, hit_index = table2dict(rows_to_hit_table(rows))
    config = PairingConfig(orientation=orientation, single_model='LTR')
    hit_index = parseHitsGeneral(hitsDict=hits_dict, hitIndex=hit_index, config=config)
    return hit_index, config


@pytest.mark.parametrize('n_elements', ELEMENT_COUNTS)
def test_table2dict(benchmark, n_elements):
    """Index a hit table by model and hit id."""
    hit_table = rows_to_hit_table(make_pairing_hit_rows(n_elements))

    result = benchmark(table2dict, hit_table)

    assert result is not None


@pytest.mark.parametrize('n_elements', ELEMENT_COUNTS)
def test_parse_hits_general_asymmetric(benchmark, n_elements):
    """Candidate discovery: the asymmetric hot path.

    This is where ``_find_candidates`` runs, and it is where nearly all of the
    asymmetric pairing cost lives.
    """
    rows = make_pairing_hit_rows(n_elements)
    hit_table = rows_to_hit_table(rows)
    hits_dict, _ = table2dict(hit_table)
    config = PairingConfig(
        orientation='F,R', left_model='LeftTIR', right_model='RightTIR'
    )

    # parseHitsGeneral mutates hitIndex, so each round is given a fresh one.
    # setup() must supply every argument: pedantic forbids mixing a returning
    # setup with declared args/kwargs.
    result = benchmark.pedantic(
        parseHitsGeneral,
        setup=lambda: (
            (),
            {
                'hitsDict': hits_dict,
                'hitIndex': table2dict(hit_table)[1],
                'config': config,
            },
        ),
    )

    assert result is not None


@pytest.mark.parametrize('n_elements', SYMMETRIC_ELEMENT_COUNTS)
def test_parse_hits_general_symmetric(benchmark, n_elements):
    """Candidate discovery for the single-model symmetric configuration."""
    rows = make_pairing_hit_rows(n_elements, left_model='LTR', right_model='LTR')
    for row in rows:
        row['strand'] = '+'
    hit_table = rows_to_hit_table(rows)
    hits_dict, _ = table2dict(hit_table)
    config = PairingConfig(orientation='F,F', single_model='LTR')

    result = benchmark.pedantic(
        parseHitsGeneral,
        setup=lambda: (
            (),
            {
                'hitsDict': hits_dict,
                'hitIndex': table2dict(hit_table)[1],
                'config': config,
            },
        ),
    )

    assert result is not None


@pytest.mark.parametrize('n_elements', ELEMENT_COUNTS)
def test_get_pairs_asymmetric(benchmark, n_elements):
    """A single reciprocal-best-match pass over distinct left/right models."""
    _, config = _asymmetric_inputs(n_elements)

    result = benchmark.pedantic(
        getPairsAsymmetric,
        setup=lambda: (
            (),
            {'hitIndex': _asymmetric_inputs(n_elements)[0], 'config': config},
        ),
    )

    assert result is not None


@pytest.mark.parametrize('n_elements', SYMMETRIC_ELEMENT_COUNTS)
def test_get_pairs_symmetric(benchmark, n_elements):
    """A single pass with one model acting as both termini."""
    _, config = _symmetric_inputs(n_elements)

    result = benchmark.pedantic(
        getPairsSymmetric,
        setup=lambda: (
            (),
            {
                'hitIndex': _symmetric_inputs(n_elements)[0],
                'model_name': 'LTR',
                'config': config,
            },
        ),
    )

    assert result is not None


@pytest.mark.parametrize('n_elements', ELEMENT_COUNTS)
def test_iterate_get_pairs_asymmetric(benchmark, n_elements):
    """The full iterate-to-stability asymmetric workflow.

    Cheap relative to the candidate discovery that precedes it; kept so that a
    change moving work into the iteration loop is still visible.
    """
    _, config = _asymmetric_inputs(n_elements)

    result = benchmark.pedantic(
        iterateGetPairsAsymmetric,
        setup=lambda: (
            (),
            {
                'hitIndex': _asymmetric_inputs(n_elements)[0],
                'config': config,
                'stableReps': 5,
            },
        ),
    )

    assert result is not None


@pytest.mark.parametrize('n_elements', SYMMETRIC_ELEMENT_COUNTS)
def test_iterate_get_pairs_custom(benchmark, n_elements):
    """The full iterate-to-stability symmetric workflow.

    The slowest operation in TIRmite by a wide margin. See the module
    docstring: ~0.96 s at 50 elements against ~0.0016 s for the equivalent
    asymmetric run, growing to ~427 s at 250 elements. The sweep is capped at
    50 for that reason.
    """
    _, config = _symmetric_inputs(n_elements)

    result = benchmark.pedantic(
        iterateGetPairsCustom,
        setup=lambda: (
            (),
            {
                'hitIndex': _symmetric_inputs(n_elements)[0],
                'config': config,
                'stableReps': 5,
            },
        ),
    )

    assert result is not None


@pytest.mark.parametrize('orientation', ['F,R', 'R,F'])
def test_asymmetric_orientation_variants(benchmark, orientation):
    """Both asymmetric orientations at a fixed size."""
    n_elements = 100
    _, config = _asymmetric_inputs(n_elements, orientation=orientation)

    result = benchmark.pedantic(
        iterateGetPairsAsymmetric,
        setup=lambda: (
            (),
            {
                'hitIndex': _asymmetric_inputs(n_elements, orientation=orientation)[0],
                'config': config,
                'stableReps': 5,
            },
        ),
    )

    assert result is not None


@pytest.mark.parametrize('orientation', ['F,F', 'R,R'])
def test_symmetric_orientation_variants(benchmark, orientation):
    """Both same-strand orientations at a fixed, deliberately small size."""
    n_elements = 25
    _, config = _symmetric_inputs(n_elements, orientation=orientation)

    result = benchmark.pedantic(
        iterateGetPairsCustom,
        setup=lambda: (
            (),
            {
                'hitIndex': _symmetric_inputs(n_elements, orientation=orientation)[0],
                'config': config,
                'stableReps': 5,
            },
        ),
    )

    assert result is not None
