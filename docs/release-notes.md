# PhosPy Release Notes

## Version 1.5.2 (2026-05-06)

## Release Overview

PhosPy keeps the supported public shape clear: build an
`AnalysisReadyPhosphoDataset`, run `KinaseWorkflow`, and optionally run
`SignalomeWorkflow` when protein identifiers are available.

This release focuses on stricter preprocessing/scientific boundaries, expanded
signalome and kinase workflow components, and stronger provenance and typing
guarantees.

## Added

- Explicit KSEA z-score activity scoring, shared threshold-membership policy
  handling, and condition-specific `activity_substrate_counts` reporting.
- FASTA-backed site-sequence resolution components with configurable conflict
  policies and durable preprocessing/kinase provenance reporting.
- Opt-in missing-data preprocessing policies for MinProb and KNN imputation,
  plus explicit forbid-path diagnostics and preprocessing readiness reporting.
- Structured identifier normalisation provenance and conflict diagnostics across
  dataset ingestion and reference-table boundaries.
- Expanded signalome clustering components (candidate scoring/selection, module
  selection, tree building, scale guards, and backend diagnostics schemas) with
  explicit policy records.
- Docs and governance additions including ADR-0016/ADR-0017, contributor
  guidance, workflow contracts, and consolidated release-notes pages.

## Changed

- Split preprocessing configuration and processing-state responsibilities into
  focused modules/packages, with an authoritative stage registry and stricter
  diagnostics parsing.
- Refactored kinase and signalome workflow orchestration into dedicated
  runner/result/provenance collaborators for clearer ownership boundaries.
- Promoted high-impact scientific/workflow behaviour toggles to explicit
  enum-backed policy models with stricter public validation boundaries.
- Expanded strict typing and CI quality gates (Pyright coverage,
  benchmark/performance contracts, and broader boundary contract tests).
- Refreshed docs, CLI docs, examples, and MkDocs structure/styling to match
  the current public API surfaces.

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

## Scientific Scope

Bundled runtime references in this release are rat-only. Human and mouse remain
valid enum values, but they require a caller-supplied `ReferenceBundle` for
workflow execution.

KSEA-style activity output is supported as a PhosPy activity method and is
reported as an explicit PhosPy method variant, not as a claim of PhosR
equivalence.

Next: [Quickstart](quickstart.md) or [API Guide](api.md).
