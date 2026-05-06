# Changelog

All notable changes to this project are documented here.

## [1.5.2] - 2026-05-06

### Added

- explicit KSEA z-score activity scoring, shared threshold-membership policy handling, and condition-specific
  `activity_substrate_counts` reporting.
- FASTA-backed site-sequence resolution components with configurable conflict policies and durable preprocessing/kinase
  provenance reporting.
- opt-in missing-data preprocessing policies for MinProb and KNN imputation, plus explicit forbid-path diagnostics and
  preprocessing readiness reporting.
- structured identifier normalisation provenance and conflict diagnostics across dataset ingestion and reference-table
  boundaries.
- expanded signalome clustering components (candidate scoring/selection, module selection, tree building, scale guards,
  and backend diagnostics schemas) with explicit policy records.
- docs and governance additions including ADR-0016/ADR-0017, contributor guidance, workflow contracts, and consolidated
  release-notes pages.

### Changed

- split preprocessing configuration and processing-state responsibilities into focused modules/packages, with an
  authoritative stage registry and stricter diagnostics parsing.
- refactored kinase and signalome workflow orchestration into dedicated runner/result/provenance collaborators for
  clearer ownership boundaries.
- promoted high-impact scientific/workflow behaviour toggles to explicit enum-backed policy models with stricter public
  validation boundaries.
- expanded strict typing and CI quality gates (Pyright coverage, benchmark/performance contracts, and broader boundary
  contract tests).
- refreshed docs, CLI docs, examples, and MkDocs structure/styling to match the current public API surfaces.

### Fixed

- deterministic provenance hashing with typed label/index handling, explicit structure hashing, and composite
  stage-hash compatibility support.
- MinProb preprocessing stability by column identity plus preserved stage ordering and row-median provenance persistence.
- duplicate/conflict handling for site identifiers and reference accession normalisation, including explicit
  post-normalisation conflict/duplicate reporting.
- stricter scientific matrix guard behaviour (forbid-policy enforcement, bool-frame rejection, and fail-fast invalid
  preprocessing metadata handling).

## [1.5.1] - 2026-04-29

### Added

- typed run provenance with deterministic table fingerprints, persisted stage-table fingerprints, and bundle
  round-trip support.
- preprocessing provenance surfaces for duplicate/conflict reporting, comparison sidecar stats, row-audit
  filtering, and site/protein context tables.
- opt-in preprocessing lanes for log2, median centering, and quantile processing with explicit provenance
  reporting.
- lightweight internal table-schema wrappers and enforced public output schemas for signalome assignments,
  modules, and network tables.
- a pluggable signalome clustering backend protocol with `exact_python` and `scipy_hierarchical` engines,
  backend-aware reporting, and parity/contract coverage.
- strict motif sequence validation, explicit total-protein identity mapping validation, and quantitative-meaning
  export audits across user-facing payloads.
- performance-contract docs, benchmark scripts, and CI regression coverage for preprocessing and signalome
  clustering.
- pinned CI constraints and a Pyright gate for core scientific/API modules.

### Changed

- Split public API config models into focused modules and replaced overloaded clustering/backend configuration
  terminology with explicit engine/policy naming.
- Replaced `TransformationState` with explicit intensity-scale and processing-state models throughout preprocessing and
  bundle serialization.
- Refactored signalome workflow and clustering internals into focused components with centralised dispatch/diagnostics
  helpers and stricter guard ordering.
- Tightened API/workflow boundary behaviour to fail fast on invalid states (including invalid `module_count`) instead of
  silently repairing or clamping inputs.
- Consolidated user-facing docs into a flatter `docs/` layout and refreshed examples against the current public API.

### Fixed

- total-protein correction scale semantics and enforced strict versioned diagnostics state during bundle load/save
  paths.
- signalome backend/tree guard contracts, including explicit module-count provenance and failure ordering.
- score-preconditioning site-index retention and preserved undefined signalome correlation states in outputs.

## [1.5.0] - 2026-04-22

### Release Overview

- Clarified the supported product shape: build an analysis-ready dataset, run kinase scoring/prediction, and optionally
  run signalome analysis.
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
- Main user docs: `docs/quickstart.md`, `docs/api.md`, `docs/cli.md`, `docs/validation.md`, and
  `docs/output_bundles.md`.
- Scientific and maintainer docs: `docs/scientific-coverage.md`, `docs/parity.md`, `docs/performance.md`, and
  `docs/maintenance.md`.

## [1.4.0] - 2026-04-15

### Added

- Added `ReferenceBundle` contracts for bundled and caller-supplied reference handling.
- Added the strict `AnalysisReadyPhosphoDataset` boundary.
- Added a shared `KinaseWorkflow` path and broader workflow validation coverage.

### Changed

- Reorganised public workflows, dataset construction, validation, preprocessing, reference handling, and result models
  into clearer package areas.
- Made scientific preprocessing policies explicit and configurable.
- Improved prediction and signalome hot paths through dataframe-copy reduction and vectorised processing where
  appropriate.

### Fixed

- Tightened comparison validation, activity scoring validation, prediction-matrix validation, and output publisher
  recovery-path validation.

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
