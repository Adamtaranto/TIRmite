# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Common Changelog](https://common-changelog.org/).

## [Unreleased]

This release fixes several bugs that made whole subcommands unusable or silently
produced wrong results, and reorganises the package internals. Three of the
fixes change what TIRmite outputs for the same input; those are called out
below.

### Fixed

#### Subcommands that could not work

- **`tirmite seed` could never build an HMM.** `build_hmm_from_alignment_pyhmmer` validated its input with `SeqIO.parse(file, 'stockholm')`, but every seed path feeds it the FASTA that MAFFT writes, so the subcommand failed with `Did not find STOCKHOLM header` on every run that reached HMM building. The alignment format is now detected. Only the `--update` path, which aligns with `hmmalign`, produced genuine Stockholm.
- **`tirmite pair --orientation F,F` or `R,R` returned exactly one pair, regardless of input size.** Candidate ranking re-sorted the whole accumulated candidate list by the current search direction, and under the same-strand orientations both directions run for every hit, so the other direction's entries received a negative sort key and `candidates[0]` became the *farthest* hit rather than the nearest. That inverted the reciprocal-best-match test, so one pair formed per iteration and the run stopped there under the default `--stable-reps 0`. `F,F`/`R,R` is the documented LTR configuration. Every element now pairs, and the same fix removed a severe slowdown: 250 elements took 427 s and now take 0.0013 s.
- **`tirmite seed --right-seed` raised `NameError` whenever seed similarities were found.** Three stale references to a renamed loop variable in the comparison-report writer. That block is now a module-level function, so it cannot shadow anything, and it is tested.
- **`tirmite seed` crashed in its alignment fallback path.** The recovery path taken when `SequenceFile` cannot read the alignment built the MSA with `DigitalMSA(alphabet, sequences)`. `sequences` is keyword-only, so positionally it landed in `name`, which expects bytes, and the fallback raised `TypeError: Expected bytes, got list` every time it ran. All arguments are now passed by keyword, which also settles pyhmmer's warning that `alphabet` stops accepting a positional argument after v0.12.0.

#### Runs that reported success while producing nothing

- **`tirmite validate` reported success while validating nothing.** A missing MAFFT or `blastdbcmd`, a standard 12-column BLAST file passed to `--blast-results`, or a mistyped results path each produced a full summary of `0.0` errors and exit 0. Required tools are now pre-flighted, unparseable BLAST input is an error naming the required `-outfmt`, and `check_tsd_gaps` returns `None` for an inconclusive comparison instead of `0`, so it can no longer be averaged in as agreement. This is the same failure mode fixed for `compare_tsds` in 1.5.0, recurring in a different place.
- **`tirmite seed` reported success when no hit could be extracted.** Per-hit extraction failures are logged and skipped, so a BLAST database built without `-parse_seqids` yielded an empty list; the seed records were then appended, and a multi-record seed file passed the "at least 2 sequences" guard, producing an "HMM" built from the seeds aligned to each other with no genomic evidence at all.

#### Wrong results

- **`tirmite validate --tsd-in-model` was accepted and ignored, and inverted the sign of the reported error.** For an in-model TSD, over-declaring the length *inserts* element bases rather than trimming flank, so the query is longer, MAFFT gaps the target, and the tool reported "too short" for a TSD that was too long — pushing the user to lengthen it further. The junction was also assumed to sit at the query midpoint, which is wrong by half the TSD length for out-of-model sites.
- **`tirmite search --max-offset` silently discarded every valid reverse-inserted asymmetric hit.** The anchor filter existed as two independent copies, one in `tirmite pair` and one in `tirmite search`. The reverse-insertion strand-swap correction shipped in 1.5.0 reached only the `pair` copy, so the `search` copy measured the offset against each hit's *inner* model edge whenever an element was inserted in reverse. Both subcommands now share one implementation in `tirmite.core.hit_filters`.
- **`tirmite search --max-offset` mishandled symmetric models named in a pairing map.** A row naming the same feature on both sides (`SymTIR<TAB>SymTIR`) describes a symmetric element, whose model has no fixed terminus role. The `search` copy labelled such a model `left` and immediately overwrote it with `right`, so every hit was tested against one model end instead of both, and hits covering only the model's second half were wrongly retained.
- **`tirmite search` compared nested and cross-model hits only within a single strand.** Under the canonical `F,R` orientation the two termini of one element lie on opposite strands, so nested-hit removal never fired for the case it exists to handle, and a weaker cross-model hit could survive merely by sitting on the other strand.
- **Nested hits were removed using an inverted score ratio.** A nested hit was discarded unless it scored at least 1.5x *better* than the hit containing it, so equal-scoring and even better-scoring nested hits were deleted. The test now reads `enclosing / nested >= min_score_ratio`, matching the cross-model step.
- **A cluster-level pairing map did not reach the anchor filter.** The anchor filter runs before cluster merging renames hit models, so a pairing map naming clusters matched nothing. For the same-strand orientations `F,F` and `R,R` that fell back to requiring *both* model ends to be anchored, silently discarding valid termini. `F,R` runs were unaffected.
- **Minus-strand hits were never de-overlapped.** `resolve_overlapping_hits` compared raw BLAST coordinates, which are inverted for minus-strand HSPs, so the computed overlap was always clamped to zero. Every reverse-oriented element contributed redundant, near-identical copies to the HMM training set.
- **Left and right HMMs could overwrite each other.** Model names were truncated to 20 characters, so `<name>_left` and `<name>_right` collided for any name of 15 characters or more; the right model overwrote the left and both returned paths pointed at the same file.
- **A discarded left hit could still evict a right hit** in `resolve_asymmetric_conflicts`, which iterated the original rather than the surviving list.
- **With `--genome-list`, hits were extracted from the wrong genome.** Only the first genome was indexed, so hits on a same-named contig in another assembly were extracted from the first genome at the second's coordinates. All genomes are now indexed, and a contig name appearing in more than one is an error.
- **Hit coordinates were sorted lexicographically**, so `'1000'` sorted before `'999'`. Row order decides which of two mutually-eliminating hits a greedy filter keeps.
- **Overlap was understated by one base.** Hit coordinates are 1-based and inclusive, so two hits sharing a single base overlap by 1, not 0.
- **`tirmite seed` misreported the number of sequences in an alignment.** `len(msa)` is the alignment width in columns, not the sequence count, so "Building HMM for X from N sequences" reported the alignment length.
- **`detect_input_format` used the wrong nhmmer column positions** (strand at index 9 and e-value at 15, which are `env to` and the description), so a genuine nhmmer file was reported as `unknown` while `import_nhmmer` parsed it correctly.

#### Documentation

- **The `--cluster-map` documentation showed the columns in the wrong order.** `docs/tutorials/tirmite-search.md` presented the file as one `model<TAB>cluster` line per model, while the parser, the CLI help and the tests all expect `cluster<TAB>component1<TAB>component2...`. Following the tutorial produced a map in which every *model* name became a cluster, so nothing matched and the run silently returned no clustered hits.

### Changed

#### Breaking

- **The minimum supported Python is now 3.10** (was 3.9). `requires-python`, the ruff and mypy targets and the CI matrix all move together. Users on 3.9 must upgrade or pin an earlier release.
- **Usage errors now exit with status 2** rather than 1, matching argparse's convention and distinguishing "you invoked this wrongly" from "the run failed". This affects a bare `tirmite`, a bare subcommand, an unknown subcommand, and invalid or insufficient arguments to `tirmite search`. Scripts checking for a specific non-zero exit code may need updating.
- **`tirmite validate` exits 1 when it could not validate a single target site.** It previously reported success and exited 0 whether it validated 500 sites or none.
- **`tirmite search` exits 1 when the cluster map matches none of the models present in the hits**, instead of writing an empty output file and exiting 0.
- **The `validate` summary replaces `num_valid_hits` with `num_filtered_hits` and `num_compared_hits`.** The old column counted hits that passed BLAST filtering, ignoring every downstream extraction and alignment failure, so it overstated the evidence behind a verdict. `predicted_tsd_error` is now `NA` rather than `0.0` when nothing was measured.

#### Behaviour

- Running a subcommand with no arguments now prints **that subcommand's** help. Previously only a bare `tirmite` was special-cased: `tirmite search` parsed successfully, created its output directory, initialised logging, and only then reported that it had no inputs.
- `tirmite search` validates its arguments **before** creating the output directory or configuring logging, so an invocation that cannot run leaves nothing on disk.
- `tirmite validate` reads `flank_len`, `tsd_len` and `tsd_in_model` from the target-site FASTA header written by `tirmite pair`, rather than requiring the user to re-supply matching values. A `--tsd-in-model` that disagrees with the header warns and defers to the header. The `--tsd-length-map` lookup now tries both key orders, which a reverse-inserted element needs.
- Model names are capped at 60 characters rather than 20, with a digest appended when truncation is required, so distinct models cannot collide.
- `merge_overlapping_cluster_hits`'s `min_overlap` parameter is renamed `max_gap`. It was always a gap tolerance — a larger value merges *more* — which is the opposite of what the name and docstring said. The old name is accepted for one release.

#### Performance

- `tirmite pair` symmetric pairing: 427 s to 0.0013 s at 250 elements (see the F,F fix above).
- `remove_nested_paired_hits`: 93 s to 1.07 s, and `filter_best_model_per_locus` 81 s to 3.17 s, at 5,000 hits. Both were dominated by per-row `pandas.DataFrame.loc` access, not by their O(n²) structure, which is unchanged.
- `check_cross_cluster_overlaps`: 0.33 s to 0.004 s at 500 hits, by breaking out of the inner loop once hits can no longer overlap.
- `parseHitsGeneral`: 0.6 s to 0.136 s at 250 elements, by not building debug strings that are discarded when DEBUG logging is off.
- `add_blast_header` streams the BLAST table through a temporary file instead of reading it into memory. Peak allocation was ~2x the file size (21 MB for a 10.6 MB table) and is now constant. It also no longer truncates the original before rewriting, so an interrupted run can no longer destroy the search results.
- `tirmite validate` and `tirmite seed` build one BLAST database source per run instead of one per hit. `BlastDBSource` caches contig lengths because each miss costs a `blastdbcmd` subprocess, and rebuilding it discarded that cache every time: 40 subprocess spawns for 20 hits on one contig, now 21. `tirmite seed` also indexed the genome once per target-ID check *and* again for extraction — three times over for an asymmetric run — and now indexes it once, lazily.

#### Internals

- `_write_pair_summary` in `tirmite.cli.hmm_pair` no longer both counts and formats. The counting is `tirmite.report.stats.pair_summary_stats` and the formatting is `format_pair_summary`; the CLI function writes the file and returns the `PairSummary`. The HTML report consumes the same object, so the two outputs cannot report different numbers for one run. The `*_summary.txt` file is byte-identical, which a golden-string test enforces.
- `make_source` in `tirmite.cli.hmm_pair` is built whether or not the filtered hit table is empty, so the report can ask it for sequence lengths. Its arguments and behaviour are unchanged.
- `check_cross_cluster_overlaps` returns the overlaps it finds instead of only logging them. Its log output is unchanged, including the ten-warning cap.
- `MultiFastaSource` in `tirmite.utils.extract`, a sequence source spanning several indexed genomes and dispatching on sequence name. A run over `--genome-list` can now read every hit for its report rather than only those in the first genome.
- `SearchFilterSummary` gains per-stage hit counts, the cross-cluster overlaps and the per-component hit counts taken before cluster merging. `log_filter_summary` does not read them and its output is unchanged. `_process_hits` accepts one to populate; without it the pipeline behaves exactly as before.
- The parts of report building that do not depend on pairing — contig lengths, contig-end clipping, row stacking, per-contig and per-model summaries — moved to `tirmite.report.collect.ReportAccumulatorBase`, shared by the pairing and search accumulators.
- The core algorithms moved out of the 5,196-line `tirmite/tirmitetools.py` into focused modules under `tirmite/core/`: `parsers`, `filters`, `termini`, `extraction`, `flanks`, `tsd`, `pairing`, `output` and `hit_filters`. `tirmite.tirmitetools` remains as a re-export shim, so both `from tirmite.tirmitetools import X` and `import tirmite.tirmitetools as tirmite` keep working unchanged. New code should import from the specific `tirmite.core` module.
- Logging uses a named logger per module (`logging.getLogger(__name__)`) instead of the root-logger convenience functions, across ~650 call sites, so records carry the emitting module name and per-module verbosity is possible.
- `init_logging` configures the `tirmite` logger rather than the root logger, and no longer clears existing root handlers — importing TIRmite and calling a `main()` function previously replaced the handlers of whatever application was hosting it. Library modules attach a `NullHandler`, so importing TIRmite emits nothing until logging is configured, and propagation stays enabled so host handlers still receive TIRmite's records. Repeated calls no longer duplicate console output.
- `compute_outer_edge_offset` and `filter_hits_by_anchor` now live in `tirmite.core.hit_filters` and are re-exported from `tirmite.cli.hmm_pair` and `tirmite.cli.ensemble_search`, so existing imports keep working. The filter gained an `on_missing_length` parameter: `tirmite search` passes `'raise'`, `tirmite pair` keeps the previous `'warn'`.
- Argparse validators, `cleanID`, `extract_model_name_from_path` and the MAFFT wrappers each existed as two or three diverged copies. They are now single implementations in `tirmite.cli._argtypes`, `tirmite.utils.utils` and `tirmite.runners.mafft`. The validators previously reported different error messages for the same bad input depending on the subcommand.
- `tirmite.cli.validate.run_blastn` is renamed `_run_validation_blastn`, removing a name collision with `tirmite.runners.runBlastn.run_blastn`. The two are not interchangeable: the runners version applies `-word_size 4 -perc_identity 60`, which are wrong for validation.
- `extract_hit_sequence` takes a sequence source rather than a database path.

### Added

- **`tirmite pair --report` writes a self-contained HTML report.** One file per run, with no external stylesheets, scripts, fonts or images, so it opens from disk and still works years later out of an archive. It contains:
  - Genome-browser annotation tracks for every sequence carrying a hit, with scroll-to-zoom, drag-to-pan and keyboard navigation. Each annotation shows strand as an arrow, pairing group as colour, and a thin line linking the two termini of an element. Clicking a paired terminus opens its element sequence with copy-to-clipboard and reverse-complement buttons.
  - Stacked alignments of every hit to each terminus model, aligned with MAFFT when it is on `PATH` and placed by model coordinates otherwise.
  - Distribution plots — element lengths, pairing outcome, model coverage against the `--mincov` setting, hits per sequence, and hit significance — as print-ready inline SVG.
  - Sortable statistics tables covering pairing groups, models, sequences and the filters applied.

  Two distinctions the report is careful about. A hit that only partially matched its model is drawn with a torn edge; one that is short because the contig ended gets a contig-end cap instead, because those are different findings. In the alignment panels the same split appears as grey (a model position the hit did not match, though the genome continues) versus blank (sequence that does not exist). The report also reports model coverage and the span-based coverage `--mincov` filters on as separate, labelled quantities — they differ, and the latter can exceed 1.

  Report generation runs after every other output is on disk and its failures are logged rather than raised, so a broken visualisation cannot cost a completed analysis. Options: `--report-out`, `--report-title`, `--no-report-sequences`, `--report-max-seq-mb`, `--report-msa`, `--report-msa-max-rows`, `--report-max-hits`, `--report-max-rows`.
- **`tirmite search --report` writes a self-contained HTML report**, reusing the same annotation tracks, alignment panels and table machinery as the pairing report. It answers the questions a search run raises rather than a pairing run's:
  - **Per-query statistics** — hits, sequences, model coverage and contig-end truncation for every query that matched.
  - **Cross-matches** — every locus where two groups both hit, as a group-by-group matrix and a table of individual loci. `check_cross_cluster_overlaps` logs only its first ten warnings so a bad cluster map cannot bury a run; it now also returns the full list, which is what the report shows.
  - **Clustering** — what each cluster was assembled from, how many hits merging collapsed, and which model won each contested locus from the nested-hit and cross-model filter steps.
  - **Filtering** — hits remaining after each pipeline stage with the change at each step, so a run that lost everything shows where. Previously reconstructable only by subtracting log lines.

  Query alignment panels are available with `--report-msa`, off by default because they re-read sequence for every hit. Supplying `--genome` is worth it either way: the report indexes it for true sequence lengths on the track axes, which are otherwise estimated from the hits.

  Cluster merging renames every hit to its cluster while a lengths file is keyed by component, which would leave model coverage and its jagged edges unavailable for every clustered run. A cluster's length is now resolved from its components when they all agree; when they differ there is no single denominator — a merged hit's alignment coordinates come from whichever component scored best, and that component is not recoverable from the merged row — so those clusters are left without one and the report says so. Options: `--report-out`, `--report-title`, `--report-msa` (off by default here, since panels re-read sequence a search run has no other reason to touch), `--report-msa-max-rows`, `--report-max-hits`, `--report-max-rows`.
- `jinja2` and `matplotlib` as runtime dependencies, for the report. Both are imported lazily, so a run that builds no report does not pay for them.
- `--min-score-ratio` and `--merge-max-gap` for `tirmite search`. The score ratio was documented as configurable but was hardcoded at 1.5 and unreachable from the CLI.
- `--logfile` and a full `--loglevel` choice list for `tirmite legacy` and `tirmite seed`. The other three subcommands already had them; these two could not write a log file at all.
- A progress bar over the per-genome loop in `tirmite search`, shown only when searching more than one genome and only when stderr is a terminal.
- `parse_cluster_mapping` warns when a cluster map looks like it has its columns transposed (every line exactly two columns, with column 2 repeating). That layout parses cleanly but matches no hit, and was what earlier documentation showed.
- A benchmark suite under `benchmarks/`, run with `pytest benchmarks/ --codspeed --no-cov` and wired to CodSpeed in CI. It is excluded from the normal test run.
- mypy and numpydoc run in CI as blocking checks. mypy was configured but had never actually run; its 9-error backlog is cleared, and the four globally-disabled error codes (`arg-type`, `attr-defined`, `union-attr`, `index`) are now enabled everywhere and suppressed only in the three CLI modules that still need them, restoring real type checking across 27 of the 30 source modules.
- mypy now checks only TIRmite's own code (`no_site_packages`), with the third-party imports declared in one place. It previously consumed whatever inline annotations happened to be installed, so its result depended on the environment rather than on the source: biopython, pyhmmer and rich ship `py.typed` while pandas and pyfaidx do not, and which is true decides whether an import raises `import-untyped` or `import-not-found`. Since CI runs Python 3.10 and local development runs 3.12, resolving different dependency versions, a routine upgrade could flip an error code and fail a blocking build for reasons unrelated to the change under test. The run now produces identical results with the dependencies installed or absent. The cost is that calls into biopython, pyhmmer and rich are no longer checked against their signatures. `warn_unused_ignores` is enabled, which this makes meaningful, and it retired 51 of the 61 `# type: ignore` comments — every remaining one is a genuine code-level suppression rather than an environment-specific import workaround.
- A PEP 561 `py.typed` marker, so downstream type checkers honour TIRmite's inline annotations.

### Removed

- 15 functions that were never called from anywhere in the package or its test suite, totalling ~630 lines. Six had carried a `DeprecationWarning` since 1.4; the other nine were silently dead. None were part of a documented API.
  - `tirmite.cli.hmm_build`: `build_hmm_from_alignment` (superseded by `build_hmm_from_alignment_pyhmmer`, which runs in-process rather than shelling out to `hmmbuild`).
  - `tirmite.runners.runBlastn`: `run_blast_batch`, `makeBlast`, `run_blast`, and the legacy `Error` class.
  - `tirmite.runners.wrapping`: `write_script_file`, `run_script_file`, `syscall`, `run_cmd`, `_write_script`, `decode`, and the legacy `Error` class. This removed ~45% of the module.
  - `tirmite.utils.utils`: `dochecks`, `isfile`, `manageTemp`, `checkUniqueID`.
  - `tirmite.utils.extract`: `get_contig_length` (a thin wrapper over `SequenceSource.contig_length`, never called).

  Both `Error` classes were raised only from functions in this list, so nothing that survives can raise them.

## [1.5.0] - 2026-08-01

### Added

- Unified sequence extraction API in `tirmite/utils/extract.py`. `FastaSource` and `BlastDBSource` share one primitive with a single documented contract: 1-based inclusive plus-strand coordinates in, uppercase plus-strand sequence out, clamped and length-verified. All extraction sites route through it.
- N-padding of regions that extend past a contig boundary, so flanks, in-model TSDs and `--padlen` windows are always the requested width. Padded records carry a `padded:<left>,<right>` field in their FASTA description. Regions that fall entirely outside a contig are still skipped.
- `--no-pad-flanks` for `tirmite pair`: truncate boundary-overrunning regions instead of padding them, restoring the previous behaviour.
- `--extend-hits-to-model` for `tirmite pair`: emit hits that only partially cover the model at full model length, extended outward by the uncovered model positions.
- `element_orientation` field in reconstructed target site headers, recording whether the element is inserted forward or reverse.
- `tests/test_extraction_equivalence.py`: a differential suite that builds a FASTA and a BLAST database from the same file and asserts byte-identical results across a coordinate matrix, both strands, soft-masking, and out-of-bounds windows.
- Test coverage for symmetric same-strand pairing, reverse-inserted asymmetric elements, and `--maxdist` semantics, none of which were previously exercised.
- BLAST+ is installed in CI for the Python 3.14 job, so the extraction equivalence tests actually run rather than skipping silently.

### Fixed

- **Extraction from a BLAST database and from an indexed FASTA returned different sequences.** Each of ~10 call sites reimplemented coordinate conversion, clamping and strand handling in its own branch, and the two had drifted:
  - With `--padlen` on a minus-strand hit, lowercase pad markers were computed in forward-genomic space but applied to a sequence `blastdbcmd` had already reverse-complemented, so they marked the wrong end and mislabelled which bases were the hit.
  - Elements and left TIRs were reverse-complemented from the BLAST path but returned plus-strand from the FASTA path. Plus strand is now canonical, matching the GFF coordinates.
  - `blastdbcmd` silently returns the entire sequence, with a zero exit status, when a range starts past the contig end. A returned-length check now rejects this rather than substituting a whole contig for the intended region.
  - `blastdbcmd` discards soft-masking while `pyfaidx` preserves it; both are now uppercased.
  - Record IDs, descriptions and failure modes are now identical between backends.
- **A TSD truncated at a contig end was reported as a perfect match.** The length mismatch caused the comparison to be skipped and written as `tsd_hamming=0`, indistinguishable from a verified duplication. Unverifiable TSDs now report `NA`, and padded positions are excluded from the comparison.
- **Symmetric `F,F` pairing never produced any pairs**, so the documented LTR use case returned nothing. Candidate discovery only ever searched in one direction when both termini are on the same strand, so reciprocity could never be established.
- **Symmetric `R,R` pairing paired hits with themselves.** The reference hit was not excluded from its own candidate list, and on the minus strand the 5'/3' swap makes the self-distance positive, so a lone hit produced a "pair" of one.
- **`--max-offset` measured the wrong model edge for reverse-inserted asymmetric elements**, silently discarding valid hits before pairing.
- **Flanks for unpaired asymmetric hits were taken from inside the element** when the element was inserted in reverse.
- A pairing-map row naming the same feature in both columns was always treated as a right terminus and never received the both-ends anchor test.
- Negative model-coverage offsets from an incorrect `--lengths-file` or mismatched HMM shifted the flank window into the hit, extracting element sequence and labelling it as flanking. Offsets are now clamped and the inconsistency is reported.
- 1-based nhmmer coordinates were passed into a 0-based `pyfaidx` slice in `tirmite seed`, dropping the first base of every extracted hit.
- `--orientation` was not upper-cased when building the pairing configuration, so `f,r` silently set both strands to `-`. Malformed values produced zero pairs with no error; they are now rejected.
- Inverted hit windows produced a nonsensical arrangement of flanks rather than being skipped.
- Corrected usage examples throughout the tutorials. Every command was checked against the CLI: `tirmite search` was documented with an older interface (`--hmm-file`, `--blast-query`, `--nhmmer-file`, `--maxeval`, `--mincov`, `--query-len`), `tirmite seed` was shown with the removed `--max-gap` and with `--blast-file` instead of `--blast-hits`, the HMM-update example omitted `--update`, several `--insertion-site` examples omitted the required `--flanks`/`--flanks-paired`, and asymmetric output filenames were given as `_LEFT.hmm`/`_RIGHT.hmm` rather than `_left.hmm`/`_right.hmm`.

### Changes

- **Documentation now builds with [Zensical](https://zensical.org/)** instead of Material for MkDocs, which has entered maintenance mode. `mkdocs.yml` is replaced by `zensical.toml`; the Markdown sources are unchanged and the site keeps the same appearance. The deploy workflow publishes via the Pages artifact rather than a `gh-pages` branch, so the repository's Pages source must be set to "GitHub Actions".
- **`--maxdist` now measures the gap between the facing inner edges of the two terminus hits**, i.e. the length of the element interior. The pairing path previously measured to the partner's far edge, adding the whole length of that terminus, and disagreed with the legacy path. Runs that leave `--maxdist` unset are unaffected; an explicit value now admits a slightly narrower set of pairs, short by roughly one terminus length. Overlapping termini are also now rejected.
- Extracted sequences are always uppercase. Soft-masking in the input genome is no longer carried through, so the lowercase convention marking `--padlen` flanks is unambiguous.
- For asymmetric elements, output files and `_L`/`_R` suffixes follow the terminus role rather than genomic order, so each model's termini stay in that model's file regardless of insertion direction. Sequences and coordinates are unchanged.
- `reconstruct_target_site` returns `None` rather than `0` for a TSD that could not be compared.
- Coverage is measured only on the Python 3.14 CI job; the other matrix jobs run without instrumentation.
- Default Python version bumped to 3.14; `isort` now runs via `ruff` rather than as a separate step.

---

## [1.4.0] - 2026-04-20

### Added

- `tirmite pair` fully implemented outer edge flank extraction and target site reconstruction accounting for hit offset from query end and accounting for TSD or DR.
- `tirmite pair` now writes output to sub-dirs per query pair, and writes summary reports.
- Experimental module `tirmite validate` to compare reconstructed insertion sites to a blast db of genomes to find natural empty sites and check if DR length prediction was correct.
- `tirmite search` uses pairing map when identifying cross-matches between paired models, and also filters nested matches to other models, also filters on hit proximity to outer edge of query.
- `--split-paired-output` option for `tirmite-search`: write left and right model hits to separate output files (`<prefix>_left_hits.tab` and `<prefix>_right_hits.tab`) based on the pairing map. Requires `--pairing-map`.
- `filter_hits_to_pairing_map_models` function: retain only hits from models listed in the pairing map, discarding hits from unrecognised models before downstream filtering steps.
- `SearchFilterSummary` dataclass and `log_filter_summary` function: accumulate and report structured hit-filtering statistics across all pairing map pipeline steps (model exclusion, nested hit removal, cross-model overlap removal).

### Fixed

- Fixed max-offset anchor filter for same-strand symmetric (F,F or R,R) model pairs when no pairing map is provided.
- Fixed max-offset anchor filter for asymmetrical model pairs in `tirmite-search`.
- Fixed writing flanks for external element edges only.

### Changes

- Breaking changes in several cmd line args. Standardised to kebab-case format.

---

## [1.3.0] - 2026-01-14

### Added

- Support for asymmetrical termini where right and left ends of elements are conserved but distinct from one another
- New entry points for modular workflow:
  - `tirmite-build`: Build initial HMM from seed sequences
  - `tirmite-pair`: Run pairing algorithm on nhmmer search results
  - `tirmite-legacy`: Run full workflow (original behavior)
- Numpydoc-style docstrings for all functions and modules
- Type hints with mypy validation throughout codebase
- pyhmmer integration for HMM searches (replacing direct nhmmer calls)
- rich library for enhanced logging output
- New dependencies: pyfaidx, pyhmmer, rich, mafft
- Threading support for BLAST operations
- Log file output option
- GitHub Actions for CFF validation and automatic citation file updating

### Changed

- Refactored CLI structure into modular entry points (cli.py, hmm_build.py, hmm_pair.py, legacy.py)
- Reorganized codebase with new directory structure:
  - `src/tirmite/cli/`: Command-line interface modules
  - `src/tirmite/runners/`: Wrapper modules for external tools
  - `src/tirmite/utils/`: Utility functions and logging
- Manual model orientation setting for pairing of asymmetric models:
  - LTR: symmetric F,F
  - TIR: symmetric F,R
  - Starship: asymmetric F,R
- Renamed `cores` parameter to `threads` for consistency
- Updated pre-commit hooks configuration
- Updated GitHub Actions workflows for pytest and ruff linting
- Improved error handling and logging throughout codebase
- Write nhmmer results to .out file
- Write contig descriptions in output files

### Removed

- Removed `getTimestring` utility function
- Removed unused `parse_seqids` function
- Removed `env_osx64.yml` environment file
- Removed `CODE_OF_CONDUCT.md`

### Fixed

- Resolved mypy type errors in CLI modules related to Optional[Namespace] false positives
- Fixed reverse strand asymmetric pairing logic
- Fixed import and formatting issues identified by ruff
- Fixed error chaining in exception handling

### Security

- License changed from MIT to GPL-3.0-or-later

---

## [1.0.0] - 2017-10-02

Initial stable release with core TIRmite functionality.
