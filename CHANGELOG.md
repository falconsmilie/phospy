# Changelog

All notable changes to this project are documented here.

## [Unreleased]

### Changed

- Consolidated user-facing docs into a flatter `docs/` layout and refreshed examples against the current public API.

## [1.5.0] - 2026-04-22

### Release Overview

- Clarified the supported product shape: build an analysis-ready dataset, run kinase scoring/prediction, and optionally run signalome analysis.
- Kept the public API centred on `run(request)` workflow entrypoints.
- Made `phospy.api` the main import surface for requests, configs, results, enums, references, and public exceptions.

### Scientific Scope

- Bundled runtime references are rat-only in this release.
- `ReferencePreset.AUTO` resolves bundled data for rat datasets.
- Human and mouse workflows require a caller-supplied `ReferenceBundle`.
- Kinase scoring and prediction keep fixture-backed parity checks on explicit comparison surfaces.
- Activity output is documented as thresholded substrate-mean activity and weighted activity, not full KSEA enrichment.
- Signalome requires explicit `site_metadata.protein_id` for every interpreted site.

### Documentation

- Release notes: `docs/release-notes-1.5.0.md`.
- Main user docs: `docs/quickstart.md`, `docs/api.md`, `docs/cli.md`, `docs/validation.md`, and `docs/output_bundles.md`.
- Scientific and maintainer docs: `docs/scientific-coverage.md`, `docs/parity.md`, `docs/performance.md`, and `docs/maintenance.md`.

## [1.4.0] - 2026-04-15

### Added

- Added `ReferenceBundle` contracts for bundled and caller-supplied reference handling.
- Added the strict `AnalysisReadyPhosphoDataset` boundary.
- Added a shared `KinaseWorkflow` path and broader workflow validation coverage.

### Changed

- Reorganised public workflows, dataset construction, validation, preprocessing, reference handling, and result models into clearer package areas.
- Made scientific preprocessing policies explicit and configurable.
- Improved prediction and signalome hot paths through dataframe-copy reduction and vectorised processing where appropriate.

### Fixed

- Tightened comparison validation, activity scoring validation, prediction-matrix validation, and output publisher recovery-path validation.

## [1.2.3] - 2026-04-09

### Added

- Added public workflow coverage for prediction matrix generation and signalome construction.
- Added structured signalome result outputs for downstream analysis.
- Added parity documentation and benchmark tooling for supported comparison surfaces.

### Changed

- Reorganised validation by validation type.
- Clarified prediction mode behaviour and release thresholds.
- Added smoke coverage for documented public workflow examples.

### Fixed

- Repaired benchmark output path handling.
- Updated README and workflow examples to match the supported API at the time.

## [1.2.1] - 2026-04-04

### Documentation

- Simplified README and docs pages so the main workflows were easier to find.
- Aligned API docs with current public classes, methods, CLI options, and output files.

## [1.0.0] - 2026-03-26

### Scope

- First supported PhosPy release.
- Covered core preprocessing, kinase analysis, a small public API, and CLI support for the initial workflow lane.
