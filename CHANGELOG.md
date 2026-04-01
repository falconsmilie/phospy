# Changelog

All notable changes to this project are documented here.

## [1.1.1]

Current package metadata version reflected in this repository snapshot.

### Documentation and release-surface alignment

- documented that file-based total and phospho inputs are read as tab-delimited text tables, while `predMat` is read as CSV with a phosphosite index column
- documented the pipeline output bundle more completely, including `run_manifest.json`
- updated version-specific wording in the parity and roadmap docs so it matches the current supported surface
- aligned the licensing note in `NOTICE.md` with the package metadata

## [1.0.0] - 2026-03-26

First supported PhosPy release.

### Scope

PhosPy 1.0.0 covers:

- core preprocessing from total and phospho inputs to corrected phosphosite matrices
- downstream kinase analysis from `predMat`
- a native `KinaseWorkflow` with seam-level validation against committed references
- a small supported root-level public API
- CLI support for the core preprocessing plus `predMat` path

### Added

- a clear 1.0.0 scope statement in the README
- a documented supported public API
- validation and parity guides
- a public roadmap for the next likely expansion areas
- a runnable native workflow example at `examples/native_workflow_demo.py`
- an end-to-end smoke test for the documented example pipeline workflow

### Changed

- tightened the README quick-start around the bundled example data
- clarified `KinaseWorkflow` parity wording across the docs
- documented the release gate as `pre-commit`, the non-parity suite, and the parity suite
- refreshed the packaged example outputs to match current CLI behaviour
- simplified package metadata around supported 1.0.0 dependencies and extras

### Known Limitations

- PhosPy 1.0.0 is a selective subset, not a full PhosR replacement.
- Parity claims are limited to committed fixture-backed seams.
- `KinaseWorkflow` is native first; `svm_mode="r_parity"` narrows learner-seam comparison but does not imply full numerical equivalence.
- The CLI does not yet expose the native kinase workflow.
