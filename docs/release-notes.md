# PhosPy Release Notes

## Version 1.6.0

These notes describe the next planned release. The package metadata and
citation file still identify the current packaged release as `1.5.2` until the
release version is bumped.

## Release Overview

The draft release is centred on stricter phosphosite identity, broader workflow
coverage, and clearer scientific boundaries. Analysis-ready datasets now use
protein-scoped encoded `site_key` values as row identity while preserving
`display_id` for human-readable reporting. The public workflow surface also
expands with importer support, richer differential designs, batch
residualisation, protein-aware preparation, Kinase Library motif scoring,
ssGSEA-style kinase activity, and offline enrichment ORA.

## Major Additions

- Protein-scoped phosphosite `site_key` row identity, derived `display_id`
  metadata, builder support, identity validators, and regression coverage
  across dataset construction, kinase, differential, signalome, bundle
  reconstruction, and public output tables.
- Local reference-bundle builder support, stricter reference manifest and
  provenance validation, approved human/mouse reference-bundle handling, and
  tightened rat bundled-reference licence/provenance metadata.
- Kinase Library-style reference schema loading, motif scoring, reference
  display-ambiguity policy handling, and KinaseWorkflow integration for
  caller-supplied local `KinaseLibraryResource` values. This is not an official
  Kinase Library implementation and does not bundle official Kinase Library
  data.
- ssGSEA-style kinase substrate-enrichment activity scoring, alongside
  method-specific activity diagnostics and preferred method-neutral
  `activity_matrix` and `substrate_count_matrix` outputs.
- MaxQuant phosphosite and FragPipe/PTMProphet importers built on a shared
  phosphosite importer foundation with validation, realistic fixtures,
  dataset handoff tests, docs, and examples.
- Fixed-effect differential design support covering paired policy defaults,
  fixed block terms, design-matrix construction, validation, execution, and
  provenance.
- Linear batch-residualisation preprocessing with configuration, metadata
  resolution, validation, reports, and dataset preprocessing integration.
- Protein-aware preparation preprocessing with site mapping, sample-alignment
  diagnostics, result models, executable workflow docs, and boundary tests.
- Native offline enrichment analysis with gene/PTM set collections, GMT/CSV/TSV
  readers, ORA execution, p-value correction, request/result contracts,
  validation, and public workflow docs.
- Imputation observation metadata, deterministic KNN imputation regression
  coverage, and imputation-aware differential safety/support contracts.
- ADR-0023 through ADR-0028 plus expanded importer, reference-bundle,
  enrichment, workflow, scientific-scope, and testing documentation.

## Changed

- Analysis-ready datasets and workflows operate on encoded `site_key` indexes
  while preserving `display_id` in site-level outputs for reporting.
- Dataset convention normalisation and site-sequence resolution were split into
  focused collaborators for metadata normalisation, total-matrix handling,
  conflict resolution, diagnostics, request building, and reporting.
- Technical-replicate aggregation moved out of the differential validator and
  into dedicated workflow aggregation logic.
- Differential computation, interpretation, public result construction, and
  provenance ownership are now separated by clearer workflow contracts.
- Public request/result validation is stricter, including duplicate
  sample-metadata rejection, duplicate site-key resolution policies,
  DataFrame/payload snapshot helper constraints, and output-table identity
  checks.
- `protein_id` is optional source/workflow metadata at the base dataset
  boundary. Completeness is enforced only by workflows that require protein
  grouping, such as signalome.
- Workflow API docs were split into dedicated kinase, differential, signalome,
  and enrichment pages. README, quickstart, examples, MkDocs navigation,
  maintenance docs, scientific-scope docs, and testing audit assets were
  refreshed.
- Strict typing, Pyright formatting, CI performance contracts, and
  unsupported-support-claim tests were expanded.

## Deprecations

- `load_enrichment_sets_gmt`, `load_enrichment_sets_table`,
  `load_enrichment_sets_csv`, and `load_enrichment_sets_tsv` now emit
  `DeprecationWarning`; use the matching `read_enrichment_sets_*` functions.
- `KinaseActivityResult.activity_scores` and
  `KinaseActivityResult.weighted_activity` now emit `DeprecationWarning`; use
  `KinaseActivityResult.activity_matrix`.
- Bare motif sequence strings passed to
  `build_motif_library_from_sequences(...)` now emit `DeprecationWarning`.
  Use `ExplicitMotifSequence` values or structured mapping entries with stable
  `reference_id`, optional `site_id`, `kinase`, and `sequence` fields.
- `PeptideToSiteAggregationConfig(strategy="compat_best_p_value")` now emits
  `DeprecationWarning`. Keep it only for historical minimum-p-value
  reproduction and migrate new analyses to `fixed_effect_meta`,
  `random_effect_meta`, or `stouffer_z` according to the intended uncertainty
  model.
- `TechnicalReplicateResolver` now emits `DeprecationWarning` when constructed
  or run. Normal callers should use `DifferentialAnalysisWorkflow`; low-level
  callers should use `TechnicalReplicateAggregationPlanner` plus
  `TechnicalReplicateAggregator`.

## Removed

- BREAKING: legacy provenance and saved-bundle compatibility paths were
  removed. PhosPy now supports only the current provenance and bundle schemas.
  Bundles or provenance payloads generated by older development versions may
  need to be regenerated.
- Internal constructed-site-ID fallbacks, site-matrix compatibility wrappers,
  and unused signalome compatibility primitive helpers were removed.

## Fixes and Hardening

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
- Differential result invariants and numerical edge cases were hardened,
  including public result identity coherence, finite-p-value BH adjustment
  semantics, and explicit infinite-value coverage.
- Public API `__all__` compatibility exports and request validation boundaries
  were corrected.
- MaxQuant and FragPipe/PTMProphet importers were hardened with realistic
  edge-case fixtures and dataset handoff regression coverage.
- Differential imputation summaries moved behind a dataset-domain API, and
  observation-mask alignment now has direct dataset invariant coverage.
- Testing audit documentation, performance-contract expectations, bundled rat
  reference metadata, Kinase Library workflow requirements, and motif-scoring
  support-boundary docs were refreshed.

## Scientific Scope

Bundled runtime references remain rat-only for `ReferencePreset.AUTO`.
Human and mouse remain valid organisms, but workflows require an explicit
caller-supplied `ReferenceBundle` unless a future release adds approved
redistributable packaged data.

`site_key` is the analysis-ready row identity. `display_id` remains available
for interpretation and reporting, but it is not the row key and may repeat
across distinct protein contexts.

Kinase prediction, Kinase Library motif scoring, KSEA-style activity, and
ssGSEA-style substrate enrichment are explicit PhosPy workflow methods. They
do not claim calibrated causal kinase inference, full PhosR kinase-activity
equivalence, validated Kinase Library parity, or PTM-SEA support.

Enrichment support is offline ORA over caller-supplied collections and explicit
backgrounds. It does not bundle GO, KEGG, Reactome, PTM-SEA, PTMsigDB, Enrichr,
gseapy, clusterProfiler, GSEA, or online-service behaviour.

Fixed-effect batch/covariate/block differential designs and linear batch
residualisation are executable PhosPy features. They are not correlated
repeated-measure modelling, mixed-effect modelling, limma
`duplicateCorrelation`, ComBat, RUV, SPS/RUV-III, or limma
`removeBatchEffect` parity claims.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
