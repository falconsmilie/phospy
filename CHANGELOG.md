# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - 2026-03-26

This is the first supported PhosPy release.

### Scope

PhosPy 1.0.0 is an unofficial Python implementation of selected PhosR-style workflows for phosphoproteomics.

The supported scope for 1.0.0 is:

- core preprocessing from total and phospho inputs to corrected phosphosite matrices
- downstream kinase analysis from `predMat`
- a native `KinaseWorkflow` with seam-level validation against committed references
- a small supported root-level public API
- CLI support for the core preprocessing plus `predMat` path

### Added

- an explicit 1.0.0 scope statement in the README
- a documented supported public API for 1.0.0
- a dedicated validation and parity guide
- a public roadmap describing the most likely next expansion areas after 1.0.0
- a runnable native workflow example script under `examples/native_workflow_demo.py`
- an end-to-end smoke test covering the documented example pipeline workflow

### Changed

- tightened README quickstart examples around the bundled example data
- clarified `KinaseWorkflow` parity wording in the project docs
- documented the release gate as `pre-commit`, the non-parity suite, and the parity suite
- refreshed the packaged example outputs so they match current CLI behaviour
- simplified package metadata to focus on supported 1.0.0 dependencies and extras

### Known Limitations

- PhosPy 1.0.0 is a selective subset, not a full PhosR replacement.
- Parity claims are limited to committed fixture-backed seams.
- `KinaseWorkflow` is native first; `svm_mode="r_parity"` narrows learner-seam comparison but does not imply full numerical equivalence.
- The CLI does not yet expose the native kinase workflow.
