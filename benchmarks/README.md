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

### The symmetric pairing path is orders of magnitude slower

| Elements | `iterateGetPairsCustom` (F,F) | `iterateGetPairsAsymmetric` (F,R) | Ratio |
|---|---|---|---|
| 50 | 0.96 s | 0.0016 s | ~600x |
| 100 | 12.4 s | 0.0040 s | ~3,100x |
| 250 | 427 s | 0.0168 s | ~25,000x |

Growth is far worse than quadratic: a 2.5x increase in input (100 → 250) costs
35x more time. `F,F` and `R,R` are the documented LTR configuration, so this is
the path users are directed to for LTR elements, and it is effectively unusable
beyond a few hundred elements. A real LTR search on a whole genome would produce
thousands of hits.

This has not been fixed. The symmetric sweep is capped at 50 elements so the
suite stays runnable, and these benchmarks will show the improvement when
someone takes it on.

## CI

`.github/workflows/codspeed.yml` runs the suite on pushes to `main` and on pull
requests. It needs a `CODSPEED_TOKEN` repository secret, obtained by registering
the repository at [codspeed.io](https://codspeed.io). Until that secret exists
the workflow will fail at the upload step; the benchmarks themselves still run
locally without any account.
