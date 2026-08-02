# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Common Changelog](https://common-changelog.org/).

## [Unreleased]

### Changed

- Logging now uses a named logger per module (`logging.getLogger(__name__)`) instead of the root-logger convenience functions, across all ~650 call sites. Log records carry the emitting module name, and per-module verbosity can be controlled.
- `init_logging` configures the `tirmite` logger rather than the root logger, and no longer clears existing root handlers. Importing TIRmite and calling a `main()` function previously replaced the handlers of whatever application was hosting it. Library modules (`tirmite.core`, `tirmite.utils`, `tirmite.runners`) attach a `NullHandler`, so importing TIRmite emits nothing until logging is configured. Propagation stays enabled, so a host application's root handlers still receive TIRmite's records.
- `init_logging` is now safe to call more than once; repeated calls no longer duplicate console output.

### Added

- `--logfile` and a full `--loglevel` choice list for `tirmite legacy` and `tirmite seed`. The other three subcommands already had them; these two could not write a log file at all.

### Fixed

- **`tirmite search --max-offset` silently discarded every valid reverse-inserted asymmetric hit.** The anchor filter existed as two independent copies, one in `tirmite pair` and one in `tirmite search`. The reverse-insertion strand-swap correction shipped in 1.5.0 was applied only to the `pair` copy. Without it, the `search` copy measured the offset against each hit's *inner* model edge whenever an element was inserted in reverse, so genuine termini were filtered out. Both subcommands now share one implementation in `tirmite.core.hit_filters`.
- **`tirmite search --max-offset` mishandled symmetric models named in a pairing map.** A pairing map row naming the same feature on both sides (`SymTIR<TAB>SymTIR`) describes a symmetric element whose model has no fixed terminus role. The `search` copy labelled such a model `left` and then immediately overwrote it with `right`, so every hit was tested against one model end instead of both. Hits covering only the model's second half were wrongly retained.

### Changed

- `compute_outer_edge_offset` and `filter_hits_by_anchor` now live in `tirmite.core.hit_filters` and are re-exported from `tirmite.cli.hmm_pair` and `tirmite.cli.ensemble_search`, so existing imports keep working. The filter gained an `on_missing_length` parameter: `tirmite search` passes `'raise'` (asking for anchor filtering without the required model lengths is an error) while `tirmite pair` keeps the previous `'warn'` behaviour. Calling `filter_hits_by_anchor` directly now defaults to `'warn'`; the `search` workflow's behaviour is unchanged, as `_process_hits` translates the new `AnchorFilterError` into `EnsembleSearchError`.

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
