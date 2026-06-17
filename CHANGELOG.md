# Changelog

All notable changes to this project are documented here.

## [1.6.0] - 17.06.2026

### Added

- protein-scoped phosphosite `site_key` row identity, derived `display_id` metadata, and builder support for deriving
  auditable site identity from explicit protein context.
- dedicated dataset identity validators for encoded site keys, display-label handling, and protein-scoped site metadata
  coherence.
- local reference-bundle builder support, stricter reference manifests/provenance, approved human and mouse reference
  bundle handling, and explicit rat bundled-reference licence/provenance metadata.
- Kinase Library reference schema loading, motif scoring, reference display-ambiguity policy handling, and KinaseWorkflow
  integration.
- ssGSEA-style kinase substrate-enrichment activity scoring, method-specific activity diagnostics, and method-neutral
  activity result matrices.
- MaxQuant phosphosite and FragPipe/PTMProphet importers on top of a shared phosphosite importer foundation, including
  importer validation, fixtures, documentation, and examples.
- fixed-effect differential design contracts covering paired policy defaults, fixed block terms, design-matrix
  construction, validation, execution, and provenance.
- linear batch-residualisation preprocessing, including configuration, metadata resolution, validation, reports, and
  integration into dataset preprocessing.
- protein-aware preparation preprocessing, including site mapping, sample-alignment diagnostics, result models, and
  executable workflow documentation.
- native offline enrichment analysis with gene/PTM set collections, GMT/CSV/TSV readers, ORA execution, p-value
  correction, request/result contracts, validation, and public workflow docs.
- imputation observation metadata, deterministic KNN imputation regression coverage, and imputation-aware differential
  safety/support contracts.
- ADR-0023 through ADR-0028 plus expanded importer, reference-bundle, enrichment, workflow, scientific-scope, and testing
  documentation.
- regression coverage for `site_key` identity across dataset construction, kinase, differential, signalome, bundle
  reconstruction, public output-table boundaries, importers, enrichment, fixed-effect differential designs, batch
  correction, and protein-aware preparation.

### Changed

- analysis-ready datasets and workflows now operate on encoded `site_key` indexes while preserving `display_id` in
  site-level outputs for human-readable reporting.
- split dataset convention normalisation and site-sequence resolution into focused collaborators for metadata
  normalisation, total-matrix handling, conflict resolution, diagnostics, and request building.
- tightened kinase, differential, and signalome workflow contracts around explicit `site_key` identity and metadata
  coherence.
- moved technical-replicate aggregation out of the differential validator and into dedicated workflow aggregation logic.
- treated `protein_id` as optional source/workflow metadata at the base dataset boundary, with completeness enforced only
  by workflows that require protein grouping such as signalome.
- split differential computation, interpretation, public result construction, and provenance responsibilities into more
  explicit workflow contracts.
- expanded public request/result validation boundaries, duplicate sample-metadata rejection, duplicate site-key
  resolution policies, dataframe/payload snapshot helper constraints, and strict output-table identity checks.
- preferred `KinaseActivityResult.activity_matrix` and method-neutral `substrate_count_matrix` for activity outputs while
  keeping legacy accessors during the deprecation window.
- tightened bundled-reference runtime policy, PhosR-inspired scope wording, unsupported-support-claim guards, and
  semi-public `phospy.science` import policy documentation.
- reorganised workflow API documentation into dedicated kinase, differential, signalome, and enrichment pages and
  refreshed README, quickstart, examples, MkDocs navigation, maintenance docs, and testing audit assets.
- stabilised CI performance contracts and Pyright formatting/typing coverage.
- expanded strict typing across dataset, validation, prediction, differential, kinase, and signalome boundaries.

### Deprecated

- `load_enrichment_sets_gmt`, `load_enrichment_sets_table`, `load_enrichment_sets_csv`, and
  `load_enrichment_sets_tsv` now emit `DeprecationWarning`; migrate to the matching `read_enrichment_sets_*`
  functions.
- `KinaseActivityResult.activity_scores` and `KinaseActivityResult.weighted_activity` now emit
  `DeprecationWarning`; migrate to `KinaseActivityResult.activity_matrix`.
- bare motif sequence strings passed to `build_motif_library_from_sequences` now emit `DeprecationWarning`; migrate to
  `ExplicitMotifSequence` values or structured mapping entries with stable `reference_id`, optional `site_id`,
  `kinase`, and `sequence` fields.
- `PeptideToSiteAggregationConfig(strategy="compat_best_p_value")` now emits `DeprecationWarning`; keep it only for
  historical minimum-p-value reproduction and migrate new analyses to `fixed_effect_meta`, `random_effect_meta`, or
  `stouffer_z` based on the intended uncertainty model.
- `TechnicalReplicateResolver` now emits `DeprecationWarning` when constructed or run; migrate low-level callers to
  `TechnicalReplicateAggregationPlanner` plus `TechnicalReplicateAggregator`, or use `DifferentialAnalysisWorkflow` for
  normal workflow execution.

### Removed

- BREAKING: removed legacy provenance and saved-bundle compatibility paths. PhosPy now supports only the current
  provenance and bundle schemas. Bundles or provenance payloads generated by older development versions may need to be
  regenerated.
- removed internal constructed-site-ID fallbacks, site-matrix compatibility wrappers, and unused signalome compatibility
  primitive helpers.

### Fixed

- rejected unsafe display-ID fallbacks, mismatched explicit builder `site_key` values, and semantic `site_key`/metadata
  incoherence at dataset and public result boundaries.
- resolved duplicate phosphosite handling by protein-scoped `site_key` identity instead of display labels, allowing
  repeated human-readable `display_id` values under unique rows.
- preserved `site_key` and `display_id` through workflow outputs and validated encoded site keys in signalome assignment
  tables.
- accepted NumPy integer phosphosite positions in site-key building and tightened strict `position`/`site_position`
  metadata validation.
- hardened differential result invariants and numerical edge cases, including public result identity coherence.
- pinned BH multiple-testing semantics to finite p-values and added explicit infinite-value coverage.
- fixed public API `__all__` compatibility exports and made request validation boundaries explicit.
- hardened MaxQuant and FragPipe/PTMProphet importers with realistic edge-case fixtures and dataset handoff regression
  coverage.
- preserved protein-aware preparation boundaries before modelling support and constrained unsupported PhosR-style
  SPS/RUV-III, ComBat, mixed-effect, and external enrichment support claims.
- tightened bundled rat reference provenance/licence metadata, Kinase Library workflow requirements, and motif-scoring
  support-boundary documentation.
- moved differential imputation summary logic behind a dataset-domain API and added direct dataset invariant tests for
  imputation observation-mask alignment.
- regenerated testing audit documentation and stabilised performance-contract expectations.

## [1.5.2] - 2026-05-22

### Added

- explicit KSEA z-score activity scoring, shared threshold-membership policy handling, and condition-specific
  `activity_substrate_counts` reporting.
- native limma-style moderated differential analysis as a first-class workflow, including robust eBayes trend
  moderation, explicit quantitative meaning/provenance, and technical replicate aggregation policies.
- peptide-evidence and multi-site ambiguity models with policy-driven peptide-to-site aggregation integrated into
  dataset construction.
- FASTA-backed site-sequence resolution components with configurable conflict policies and durable preprocessing/kinase
  provenance reporting.
- stricter phosphosite identity/localisation contracts, sequence provenance, and workflow scientific-eligibility
  reporting surfaces.
- opt-in missing-data preprocessing policies for MinProb and KNN imputation, plus explicit forbid-path diagnostics and
  preprocessing readiness reporting.
- structured identifier normalisation provenance and conflict diagnostics across dataset ingestion and reference-table
  boundaries.
- schema-aware table readers with strict metadata/numeric parsing and explicit exact-vs-tolerance table hash metadata.
- expanded signalome clustering components (candidate scoring/selection, module selection, tree building, scale guards,
  and backend diagnostics schemas) with explicit policy records.
- expanded scientific and governance docs, including ADR-0016 through ADR-0022, testing audit assets, workflow
  contracts, and a PhosR compatibility/scope matrix.

### Changed

- reorganised domain implementation under `phospy.science` and moved internal contract ownership from `phospy.api` into
  dedicated `phospy.contracts` modules.
- split preprocessing configuration and processing-state responsibilities into focused modules/packages, with an
  authoritative stage registry and stricter diagnostics parsing.
- split dataset-builder and preprocessing orchestration responsibilities into focused collaborators with stricter
  `site_sequence` and sample-metadata contract enforcement.
- refactored kinase and signalome workflow orchestration into dedicated runner/result/provenance collaborators for
  clearer ownership boundaries.
- promoted high-impact scientific/workflow behaviour toggles to explicit enum-backed policy models with stricter public
  validation boundaries.
- expanded strict typing and CI quality gates (Pyright coverage, realistic performance/data-scale benchmark contracts,
  and broader boundary/parity regression suites).
- refreshed docs, examples, and MkDocs structure/styling to match the current public API and scientific-scope claims.

### Removed

- removed the legacy `phospy` console-script CLI entry point and retired obsolete CLI workflow docs/tests.

### Fixed

- deterministic provenance hashing with typed label/index handling, explicit structure hashing, and composite
  stage-hash compatibility support.
- MinProb preprocessing stability by column identity plus preserved stage ordering and row-median provenance persistence.
- duplicate/conflict handling for site identifiers and reference accession normalisation, including explicit
  post-normalisation conflict/duplicate reporting.
- stricter scientific matrix guard behaviour (forbid-policy enforcement, bool-frame rejection, and fail-fast invalid
  preprocessing metadata handling).
- fixed centred site-sequence validation and phosphosite identity collision handling, and removed duplicate
  analysis-ready validation paths.
- enforced established log2 intensity scale before differential `logFC` emission and prevented unaudited
  intensity-scale establishment.

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

- Release notes: `docs/release-notes.md`.
- Main user docs: `docs/quickstart.md`, `docs/api/guide.md`, `docs/cli.md`, `docs/validation.md`, and
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
