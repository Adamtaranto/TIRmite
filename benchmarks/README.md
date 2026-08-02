# TIRmite benchmarks

Performance benchmarks run with [pytest-codspeed](https://github.com/CodSpeedHQ/pytest-codspeed).

## Running them

Benchmarks live outside pytest's `testpaths`, so a plain `pytest` never runs
them. Name the directory explicitly:

```bash
# Correctness check: every benchmark runs once.
pytest benchmarks/ --no-cov

# Measured run, with the CodSpeed reporting table.
pytest benchmarks/ --codspeed --no-cov
```

`--no-cov` matters. Coverage instrumentation would otherwise be measured along
with the code under test.

All fixture data is generated synthetically from a fixed seed, so a run
measures code changes rather than input changes. The repository ships no test
data and none is needed here.

## What is measured

| File | Covers |
|---|---|
| `bench_pairing.py` | `tirmite.core.pairing`: hit indexing, candidate discovery, and both pairing engines |
| `bench_search_filters.py` | `tirmite search` hit filters, including the quadratic ones |
| `bench_extraction.py` | Sequence extraction paths and nhmmer parsing |

Setup is excluded from measurement via `benchmark.pedantic(setup=...)` wherever
the function under test mutates its input, which most of the pairing engine
does.

## Findings

Two results from the initial run are worth recording, because both contradict
the obvious guess about where time goes.

### Candidate discovery dominates pairing, not iteration

`parseHitsGeneral` (which calls `_find_candidates`) costs roughly 0.6 s for 250
asymmetric elements, against roughly 0.017 s for the `iterateGetPairsAsymmetric`
pass that follows it — about 35x. The iterate functions loop until stability and
look like the expensive part; they are not. Optimisation effort on the
asymmetric path belongs in `_find_candidates`.

Most of that 0.6 s turned out to be logging. `_check_distance` runs once per
(hit, candidate) pair and built five f-strings on every call, passed as
arguments and therefore evaluated whether or not DEBUG was enabled: 941,261
`logger.debug` calls per invocation at 250 elements. Guarding the hot call sites
with `logger.isEnabledFor` brought `parseHitsGeneral` from 0.6 s to 0.136 s at
250 elements, and from ~18 s to 7.9 s at 2000. It remains quadratic by design.

### The `tirmite search` filters were dominated by pandas row access

`remove_nested_paired_hits` and `filter_best_model_per_locus` did a
`DataFrame.loc[label]` lookup per inner-loop iteration — ~67 µs each, which is
89% and 79% of their runtime, against nanoseconds for the integer comparisons
they exist to perform. Pulling each group's columns out as numpy arrays once:

| Hits | `remove_nested_paired_hits` | `filter_best_model_per_locus` |
|---|---|---|
| 500 | 0.90 s → 0.015 s | 0.63 s → 0.028 s |
| 2000 | 14.8 s → 0.15 s | 11.1 s → 0.49 s |
| 5000 | 93 s → 1.07 s | 81 s → 3.17 s |

Both remain O(n²); only the constant changed, and neither the algorithm nor the
comparison order did. `filter_best_model_per_locus` is deliberately *not*
converted to a sorted sweep: it skips hits already marked for removal, so the
result depends on traversal order, and a coordinate-ordered sweep would silently
change which hits survive.

`check_cross_cluster_overlaps` had the same defect but is genuinely sweepable —
its group is already sorted by start, so the inner loop can break once a
candidate begins after the current hit ends. 0.33 s → 0.004 s at 500 hits, and
0.011 s at 5000: effectively linear.

### The symmetric pairing path was orders of magnitude slower — now fixed

These benchmarks were written to measure the problem and did their job: they
located it precisely enough to find the cause, which turned out to be a
correctness bug rather than an algorithmic one.

`_find_candidates` appended into a candidate list shared by both direction
searches and then re-sorted the *whole* list by the current direction. Under
`F,F` and `R,R` both searches run for every hit, so the other direction's
entries received a negative sort key — most negative for the farthest hit — and
`candidates[0]` became the globally farthest hit rather than the nearest. That
inverted the reciprocal-best-match test, so one pair formed per iteration and
convergence took M/2 rounds instead of 1.

| Elements | `iterateGetPairsCustom` before | after |
|---|---|---|
| 50 | 0.96 s | 0.0003 s |
| 100 | 12.4 s | 0.0005 s |
| 250 | 427 s | 0.0013 s |

It was also a user-visible correctness bug: with the CLI default of
`--stable-reps 0` the run stopped after that single pair, so
`tirmite pair --orientation F,F` returned one pair regardless of input size.

The symmetric sweep in `bench_pairing.py` is still capped at 50 elements. It no
longer needs to be, and raising it would give better coverage of the fixed
path.

## CI

`.github/workflows/codspeed.yml` runs the suite on pushes to `main` and on pull
requests. It needs a `CODSPEED_TOKEN` repository secret, obtained by registering
the repository at [codspeed.io](https://codspeed.io). Until that secret exists
the workflow will fail at the upload step; the benchmarks themselves still run
locally without any account.
