# Changelog

All notable changes to this project are documented here.

## [1.2.3] - 2026-04-09

### Added

- a first-class public workflow for **PredMat generation**, giving users one supported path from phosphosite inputs and sequence data to scoring outputs, prediction outputs, and PredMat results.
- a dedicated public workflow for **Signalome construction** from scoring and prediction outputs.
- a stable **Signalome result contract** with structured access to signalome outputs and downstream export-friendly data.
- **signalome map-ready output generation** for downstream visualisation workflows.
- **kinase-network output generation** derived from signalome results.
- a user-facing **Signalome workflow guide and example** covering the supported path from scoring and prediction to signalomes, map data, and network outputs.
- a central **parity contract matrix** in `docs/parity.md` covering fixture families, protected seams, protected metric classes, and supported modes.
- a reproducible benchmark harness for comparing `svm_mode="default"` and `svm_mode="r_parity"`.
- an ADR recording the public support decision for `r_parity` and clarifying the intended role of both public presets.
- parity-sensitive release review guidance so mode changes are reviewed against explicit benchmark and regression expectations.

### Changed

- Refactored the `validation/` package to group validators by **validation type** rather than by validated process or workflow.
- Reorganised validation into clearer type-based modules, improving discoverability and reducing process-shaped validator sprawl.
- Promoted reusable validation rules into more obvious public validation functions and reduced private helper clutter across the package.
- Clarified the public prediction mode story across the documentation:
  - `default` remains the recommended standard mode.
  - `r_parity` remains a supported parity-oriented preset.
- Made release thresholds for `default` and `r_parity` explicit in the parity documentation instead of leaving them mainly embedded in test assertions.
- Centralised parity threshold constants in regression tests to make the release bar easier to audit.
- Hardened documented workflow smoke coverage so the public PredMat and Signalome examples run in both supported prediction modes where relevant.
- Improved consistency between documentation, benchmark workflow, release guidance, and protected parity test contracts.

### Fixed

- Restored the mode-comparison benchmark harness to the repository so documented benchmark commands point to a real in-repo tool.
- Fixed benchmark output path handling so generated reports can be written outside the repository root.
- Fixed README workflow examples and mode guidance so the documented public paths align with the current API and release expectations.
- Fixed documentation gaps around PredMat and Signalome usage by adding clearer recommended workflow paths.
- Reduced duplication and ambiguity in validation logic by consolidating rules under the new validation module structure.

### Testing

- Full regression suite passes: `384 passed`.
- Benchmark threshold checks pass for both `default` and `r_parity`.
- Public workflow smoke tests pass for documented PredMat and Signalome flows.
- Validation tests were updated and continue to pass against the reorganised validation package.
- End-to-end coverage now includes the supported public PredMat and Signalome workflows.

### Notes

- Benchmark harness code is versioned in the repository.
- Generated benchmark reports remain local review artifacts and are not committed by default.
- This release significantly improves the package’s public workflow surface, downstream analysis support, validation maintainability, and parity release confidence without changing the recommended default prediction mode.

## [1.2.1] - 2026-04-04

Current package metadata version reflected in this repository snapshot.

### Documentation

- simplified the README and docs pages so the main workflows are easier to find
- aligned API docs with the current public classes, methods, CLI options, and output files
- clarified where `predMat` validation happens in `PhosRPipeline`
- tightened validation and parity docs to focus on user-facing behaviour

## [1.0.0] - 2026-03-26

First supported PhosPy release.

### Scope

PhosPy 1.0.0 covers:

- core preprocessing from total and phospho inputs to corrected phosphosite matrices
- downstream kinase analysis from `predMat`
- a native `KinaseWorkflow` with seam-level validation against committed references
- a small supported root-level public API
- CLI support for the core preprocessing plus `predMat` path
