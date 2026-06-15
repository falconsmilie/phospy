# PhosPy Release Notes

## Unreleased

### Deprecations

- Bare motif sequence strings passed to
  `build_motif_library_from_sequences(...)` now emit `DeprecationWarning`.
  Use `ExplicitMotifSequence` values or structured mapping entries so motif
  library validation records stable reference identity and, when available,
  site-residue identity.
- `PeptideToSiteAggregationConfig(strategy="compat_best_p_value")` now emits
  `DeprecationWarning`. It remains available for historical minimum-p-value
  reproduction during the warning window, but new analyses should select
  `fixed_effect_meta`, `random_effect_meta`, or `stouffer_z` according to the
  intended site-level uncertainty model.

## Version 1.5.3 (2026-06-08)

## Release Overview

PhosPy keeps the supported public shape clear: build an
`AnalysisReadyPhosphoDataset`, run `KinaseWorkflow`, and optionally run
`SignalomeWorkflow` when protein identifiers are available.

This release formalises protein-scoped phosphosite row identity. Analysis-ready
datasets and workflows now use encoded `site_key` values as row identity, while
preserving `display_id` as the human-readable label used for reports and
reference-facing interpretation.

The release also tightens workflow/result identity validation, hardens
differential result invariants, and splits several large dataset and
site-sequence internals into focused collaborators.

## Added

- Protein-scoped phosphosite `site_key` row identity, derived `display_id`
  metadata, and builder support for deriving auditable site identity from
  explicit protein context.
- Dedicated dataset identity validators for encoded site keys, display-label
  handling, and protein-scoped site metadata coherence.
- ADR-0024 and dataset-builder documentation covering `site_key` row identity,
  repeated `display_id` labels, workflow output contracts, and kinase reference
  display-ID mapping.
- Regression coverage for `site_key` identity across dataset construction,
  kinase, differential, signalome, bundle reconstruction, and public output
  tables.

## Changed

- Analysis-ready datasets and workflows now operate on encoded `site_key`
  indexes while preserving `display_id` in site-level outputs for
  human-readable reporting.
- Dataset convention normalisation and site-sequence resolution now use focused
  collaborators for metadata normalisation, total-matrix handling, conflict
  resolution, diagnostics, and request building.
- Kinase, differential, and signalome workflow contracts now enforce explicit
  `site_key` identity and metadata coherence.
- Technical-replicate aggregation moved out of the differential validator and
  into dedicated workflow aggregation logic.
- `protein_id` is optional source/workflow metadata at the base dataset
  boundary. Completeness is enforced only by workflows that require protein
  grouping, such as signalome.
- Strict typing expanded across dataset, validation, prediction, differential,
  kinase, and signalome boundaries.

## Fixed

- Unsafe display-ID fallbacks, mismatched explicit builder `site_key` values,
  and semantic `site_key`/metadata incoherence now fail at dataset and public
  result boundaries.
- Duplicate phosphosite handling now resolves rows by protein-scoped `site_key`
  identity instead of display labels, allowing repeated human-readable
  `display_id` values under unique rows.
- Workflow outputs preserve `site_key` and `display_id`, and signalome
  assignment tables validate encoded site keys.
- NumPy integer phosphosite positions are accepted by the site-key builder, and
  strict `position`/`site_position` metadata validation is tighter.
- Differential result invariants and numerical edge cases are hardened,
  including public result identity coherence.

## Scientific Scope

Bundled runtime references in this release are rat-only. Human and mouse remain
valid enum values, but they require a caller-supplied `ReferenceBundle` for
workflow execution.

Analysis-ready row identity is now protein-scoped `site_key`. `display_id`
remains available for interpretation and reporting, but it is not the
analysis-ready row key and may repeat across distinct protein contexts.

KSEA-style activity output is supported as a PhosPy activity method and is
reported as an explicit PhosPy method variant, not as a claim of PhosR
equivalence.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
