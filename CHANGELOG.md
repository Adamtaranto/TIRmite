# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Common Changelog](https://common-changelog.org/).

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
