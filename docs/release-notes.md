# PhosPy Release Notes

## Version 1.5.2 (2026-05-22)

## Release Overview

PhosPy keeps the supported public shape clear: build an
`AnalysisReadyPhosphoDataset`, run `KinaseWorkflow`, and optionally run
`SignalomeWorkflow` when protein identifiers are available.

This release focuses on stricter preprocessing/scientific boundaries, expanded
signalome and kinase workflow components, first-class differential analysis,
and stronger provenance and typing guarantees.

## Added

- Explicit KSEA z-score activity scoring, shared threshold-membership policy
  handling, and condition-specific `activity_substrate_counts` reporting.
- Native limma-style moderated differential analysis as a first-class workflow,
  including robust eBayes trend moderation, explicit quantitative
  meaning/provenance, and technical replicate aggregation policies.
- Peptide-evidence and multi-site ambiguity models with policy-driven
  peptide-to-site aggregation integrated into dataset construction.
- FASTA-backed site-sequence resolution components with configurable conflict
  policies and durable preprocessing/kinase provenance reporting.
- Stricter phosphosite identity/localisation contracts, sequence provenance, and
  workflow scientific-eligibility reporting surfaces.
- Opt-in missing-data preprocessing policies for MinProb and KNN imputation,
  plus explicit forbid-path diagnostics and preprocessing readiness reporting.
- Structured identifier normalisation provenance and conflict diagnostics across
  dataset ingestion and reference-table boundaries.
- Schema-aware table readers with strict metadata/numeric parsing and explicit
  exact-vs-tolerance table hash metadata.
- Expanded signalome clustering components (candidate scoring/selection, module
  selection, tree building, scale guards, and backend diagnostics schemas) with
  explicit policy records.
- Expanded scientific and governance docs, including ADR-0016 through ADR-0022,
  testing audit assets, workflow contracts, and a PhosR
  compatibility/scope matrix.

## Changed

- Reorganised domain implementation under `phospy.science` and moved internal
  contract ownership from `phospy.api` into dedicated `phospy.contracts`
  modules.
- Split preprocessing configuration and processing-state responsibilities into
  focused modules/packages, with an authoritative stage registry and stricter
  diagnostics parsing.
- Split dataset-builder and preprocessing orchestration responsibilities into
  focused collaborators with stricter `site_sequence` and sample-metadata
  contract enforcement.
- Refactored kinase and signalome workflow orchestration into dedicated
  runner/result/provenance collaborators for clearer ownership boundaries.
- Promoted high-impact scientific/workflow behaviour toggles to explicit
  enum-backed policy models with stricter public validation boundaries.
- Expanded strict typing and CI quality gates (Pyright coverage,
  realistic performance/data-scale benchmark contracts, and broader
  boundary/parity regression suites).
- Refreshed docs, examples, and MkDocs structure/styling to match the current
  public API and scientific-scope claims.

## Removed

- Removed the legacy `phospy` console-script CLI entry point and retired
  obsolete CLI workflow docs/tests.

## Fixed

- Deterministic provenance hashing with typed label/index handling, explicit
  structure hashing, and composite stage-hash compatibility support.
- MinProb preprocessing stability by column identity plus preserved stage
  ordering and row-median provenance persistence.
- Duplicate/conflict handling for site identifiers and reference accession
  normalisation, including explicit post-normalisation conflict/duplicate
  reporting.
- Stricter scientific matrix guard behaviour (forbid-policy enforcement,
  bool-frame rejection, and fail-fast invalid preprocessing metadata handling).
- Fixed centred site-sequence validation and phosphosite identity collision
  handling, and removed duplicate analysis-ready validation paths.
- Enforced established log2 intensity scale before differential `logFC`
  emission and prevented unaudited intensity-scale establishment.

## Scientific Scope

Bundled runtime references in this release are rat-only. Human and mouse remain
valid enum values, but they require a caller-supplied `ReferenceBundle` for
workflow execution.

KSEA-style activity output is supported as a PhosPy activity method and is
reported as an explicit PhosPy method variant, not as a claim of PhosR
equivalence.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
