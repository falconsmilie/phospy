# Changelog

All notable changes to this project are documented here.

## [Unreleased]

### Changed

- No unreleased entries yet.

## [1.5.0] - 2026-04-22

### Release Framing

- Presents PhosPy as the current focused product shape: one analysis-ready
  dataset boundary and one supported workflow chain (`build -> kinase ->
  optional signalome`).
- Emphasizes the shipped rewrite contract and science/reference boundaries as
  they exist now, rather than continuity with older public surfaces.
- Publishes dedicated release-facing notes for this version:
  `docs/release_notes/1.5.0.md`.

### Changed

- Narrowed `phospy.prediction` package-level default exports to the stable prediction types and moved low-level mechanics out of the default surface.
- Removed package-level compatibility re-exports from `phospy.prediction`; advanced helpers now require concrete-module imports.
- Repaired core kinase ranking comparison surfaces to explicit like-for-like contracts and restored closure-grade ranking parity gates for the L6 supported lane.
- Added explicit ranking surface-contract assertions and centralized L6 ranking gate thresholds in `tests/support/l6_prediction_parity_thresholds.py`.
- Synchronized parity governance and release-facing docs to classify the core kinase scoring/prediction lane as parity-gated for ranking on the repaired surface.

### Scope Boundaries

- Scientific claims for this release are seam-level and fixture-backed; this is
  not a blanket legacy-equivalence claim.
- Bundled runtime references remain rat-only in this release; human and mouse
  lanes require caller-supplied `ReferenceBundle`.

## [1.4.0] - 2026-04-15

### Added

- Added `ReferenceBundle` and `ReferenceProvider` contracts for bundled reference handling.
- Added an `AnalysisReadyPhosphoDataset` boundary plus a supported adapter path for analysis-ready inputs.
- Added a first bundled reference provider for the rat L6 native lane.
- Added a shared `KinaseWorkflow` common-path API.
- Added a package skeleton and root migration map to support the package reshaping work.
- Added regression tests for thin API orchestration and domain delegation.

### Changed

- Reworked public workflows into an `api` package and separated simple and advanced workflow lanes.
- Split request validation into a dedicated `requests` package and reorganised the broader `validation` package.
- Centralised preprocessing under a dedicated package and moved dataset models, prediction execution, kinase activity
  analysis, signalome analysis, and reference handling into focused subpackages.
- Split workflows into focused public modules and narrowed package-level exports to supported seams, reducing broad root
  and package re-exports.
- Flattened request plumbing, simplified trusted orchestration inputs, and renamed trusted input bundles for clearer
  orchestration boundaries.
- Unified dataset file loading behind a single validated loader path and collapsed duplicate loader normalisation and
  helper paths.
- Replaced predictor introspection with an explicit execution contract and consolidated kinase workflow execution behind
  a shared internal path.
- Made scientific preprocessing policies explicit and configurable.
- Reduced frozen-style and false-immutability construction in favour of plainer model construction.
- Optimised prediction hot paths, precomputed per-site positions during prediction execution, reduced redundant
  dataframe copying, and vectorised signalome network and expansion-path processing.
- Reviewed and refreshed the documentation to match the refactored package surface.

### Fixed

- Rejected self-comparisons and reverse-duplicate comparison definitions earlier in validation.
- Tightened low-level activity scoring validation and routed public activity scoring helpers through request validation.
- Tightened `predMat` validation and overlap defaults.
- Hardened output publisher recovery-path validation.
- Removed misleading compatibility and encoding helpers that no longer matched the supported execution path.

### Documentation

- Added ADRs for the domain refactor and the high-level adapter workflow for kinase inference.
- Finalised architecture documentation and clarified retained root exports.
- Audited remaining imports and removed obsolete compatibility layers from the documented package surface.

### Testing

- Added regression coverage for thin API orchestration and domain delegation.
- Validation-oriented changes were accompanied by broader request, scoring, comparison, and recovery-path hardening
  across the refactor.

## [1.2.3] - 2026-04-09

### Added

- a first-class public workflow for **PredMat generation**, giving users one supported path from phosphosite inputs and
  sequence data to scoring outputs, prediction outputs, and PredMat results.
- a dedicated public workflow for **Signalome construction** from scoring and prediction outputs.
- a stable **Signalome result contract** with structured access to signalome outputs and downstream export-friendly
  data.
- **signalome map-ready output generation** for downstream visualisation workflows.
- **kinase-network output generation** derived from signalome results.
- a user-facing **Signalome workflow guide and example** covering the supported path from scoring and prediction to
  signalomes, map data, and network outputs.
- a central **parity contract matrix** in `docs/parity.md` covering fixture families, protected seams, protected metric
  classes, and supported modes.
- a reproducible benchmark harness for comparing `svm_mode="default"` and `svm_mode="r_parity"`.
- an ADR recording the public support decision for `r_parity` and clarifying the intended role of both public presets.
- parity-sensitive release review guidance so mode changes are reviewed against explicit benchmark and regression
  expectations.

### Changed

- Refactored the `validation/` package to group validators by **validation type** rather than by validated process or
  workflow.
- Reorganised validation into clearer type-based modules, improving discoverability and reducing process-shaped
  validator sprawl.
- Promoted reusable validation rules into more obvious public validation functions and reduced private helper clutter
  across the package.
- Clarified the public prediction mode story across the documentation:
    - `default` remains the recommended standard mode.
    - `r_parity` remains a supported parity-oriented preset.
- Made release thresholds for `default` and `r_parity` explicit in the parity documentation instead of leaving them
  mainly embedded in test assertions.
- Centralised parity threshold constants in regression tests to make the release bar easier to audit.
- Hardened documented workflow smoke coverage so the public PredMat and Signalome examples run in both supported
  prediction modes where relevant.
- Improved consistency between documentation, benchmark workflow, release guidance, and protected parity test contracts.

### Fixed

- Restored the mode-comparison benchmark harness to the repository so documented benchmark commands point to a real
  in-repo tool.
- Fixed benchmark output path handling so generated reports can be written outside the repository root.
- Fixed README workflow examples and mode guidance so the documented public paths align with the current API and release
  expectations.
- Fixed documentation gaps around PredMat and Signalome usage by adding clearer recommended workflow paths.
- Reduced duplication and ambiguity in validation logic by consolidating rules under the new validation module
  structure.

### Testing

- Full regression suite passes: `384 passed`.
- Benchmark threshold checks pass for both `default` and `r_parity`.
- Public workflow smoke tests pass for documented PredMat and Signalome flows.
- Validation tests were updated and continue to pass against the reorganised validation package.
- End-to-end coverage now includes the supported public PredMat and Signalome workflows.

### Notes

- Benchmark harness code is versioned in the repository.
- Generated benchmark reports remain local review artifacts and are not committed by default.
- This release significantly improves the package’s public workflow surface, downstream analysis support, validation
  maintainability, and parity release confidence without changing the recommended default prediction mode.

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
