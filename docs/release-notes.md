# PhosPy Release Notes

## Version 1.6.0

These notes describe Version 1.6.0. The package metadata in `pyproject.toml`
and the citation metadata in `CITATION.cff` both declare `1.6.0`.

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
- Local reference-bundle builder support with separate upstream
  `source_version` and local PhosPy `reference_version` fields. When
  `reference_version` is omitted, the builder emits a deterministic
  `local-snapshot-sha256-...` value derived from the two local source-file
  SHA-256 fingerprints.
- Stricter reference manifest and provenance validation, caller-supplied
  human/mouse reference-bundle handling, and tightened rat bundled-reference
  licence/provenance metadata.
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
- ADR-0023 through ADR-0029 plus expanded importer, reference-bundle,
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
- `protein_group_id` is Signalome-owned grouping metadata. It is optional at
  the base dataset boundary, and Signalome enforces completeness when it needs
  protein grouping. Legacy `protein_id` remains a migration alias for that
  Signalome field.
- Workflow API docs were split into dedicated kinase, differential, signalome,
  and enrichment pages. README, quickstart, examples, MkDocs navigation,
  maintenance docs, scientific-scope docs, and testing audit assets were
  refreshed.
- Strict typing, Pyright formatting, CI performance contracts, and
  unsupported-support-claim tests were expanded.
- Peptide-to-site differential evidence keeps the preferred supported lane:
  peptide evidence resolution at sample-intensity level before
  `DifferentialAnalysisWorkflow`. The post-hoc peptide differential
  estimate-combination lane is withdrawn from public support under
  `unsupported_withdrawn_posthoc_estimate_combination_v1`; its retained
  compatibility shell fails closed because coherent combined effect/inference
  and executable mapping semantics are not implemented.

## Deprecations

- `load_enrichment_sets_gmt`, `load_enrichment_sets_table`,
  `load_enrichment_sets_csv`, and `load_enrichment_sets_tsv` now emit
  `DeprecationWarning`; use the matching `read_enrichment_sets_*` functions.
- `KinaseActivityResult.activity_scores` and
  `KinaseActivityResult.weighted_activity` now emit `DeprecationWarning`; use
  `KinaseActivityResult.activity_matrix`.
- `KinaseActivityResult.legacy_condition_statistics_table_dataframe()` emits
  `DeprecationWarning`; use statistics tables with `profile_id`.
- Bare motif sequence strings passed to
  `build_motif_library_from_sequences(...)` now emit `DeprecationWarning`.
  Use `ExplicitMotifSequence` values or structured mapping entries with stable
  `reference_id`, optional `site_id`, `kinase`, and `sequence` fields.
- Post-hoc peptide-to-site differential estimate combination is withdrawn
  rather than retained as a supported typed route. Future public support
  requires executable mapping semantics, a coherent combined estimand, an
  inferential result, dependence handling, multiple-testing semantics,
  provenance semantics, docs, and tests.
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
- The old public peptide-to-site differential aggregation shell is no longer a
  supported production route. Compatibility access is internal/experimental and
  fails closed before calculation.

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
- Post-hoc peptide-to-site differential estimate combination no longer executes
  through public compatibility access. Mapping policies such as equal splitting
  or statistical-model exclusion must not silently execute as ordinary evidence.
- Public API `__all__` compatibility exports and request validation boundaries
  were corrected.
- MaxQuant and FragPipe/PTMProphet importers were hardened with realistic
  edge-case fixtures and dataset handoff regression coverage.
- Differential imputation summaries moved behind a dataset-domain API, and
  observation-mask alignment now has direct dataset invariant coverage.
- Testing audit documentation, performance-contract expectations, bundled rat
  reference metadata, Kinase Library workflow requirements, and motif-scoring
  support-boundary docs were refreshed.
- Reference-bundle release hardening now validates exact bytes in the source
  tree and built wheels, rejects unknown manifest/evidence fields, rejects
  non-Boolean or JSON null raw `redistribution_allowed` values, requires an
  explicit `verified_at` date for approved bundled evidence, and treats file
  hashes as integrity checks rather than redistribution approval.
- Release artifact verification now installs both the built wheel and sdist in
  isolated temporary environments outside the checkout, verifies installed
  import origins, checks bundled rat reference resources against manifest
  SHA-256 values, and executes representative public dataset, differential,
  and kinase workflow contracts.

## Scientific Scope

Bundled runtime references remain rat-only for `ReferencePreset.AUTO`. The only
approved packaged runtime reference is the exact rat `l6_native` snapshot
derived from upstream PhosR 1.20.0 package data. Its approval is scoped only to
the exact packaged files in that committed PhosPy snapshot. It does not approve
future bundles, other rat bundles, other organisms, or future PhosR/PhosPy
snapshots, and it does not claim independent direct permission from
PhosphoSitePlus, PRIDE, Kinase Library, or another upstream scientific
database.

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
Typed selected/background identifier-set provenance is optional and does not
change enrichment statistics. PhosPy-derived quantitative provenance preserves
typed derivation metadata plus declared-versus-observed input intensity-scale
evidence: declared-only evidence emits the documented caveat, while observed
transformation evidence is recorded without that caveat.

Kinase and signalome provenance distinguish causal site-row attrition from
compatibility metrics. `row_attrition` records stage-local site-row removals
only when rows are actually removed by that stage. `row_attrition_metrics`
remains available for legacy diagnostics, including kinase site/kinase-pair
counts, but pair loss is not encoded as site-row loss.

Fixed-effect batch/covariate/block differential designs, linear batch
residualisation, and native SPS/RUV-style preprocessing correction through
`SpsRuvBatchCorrectionConfig` are executable PhosPy features. They are not
correlated repeated-measure modelling, mixed-effect modelling, limma
`duplicateCorrelation`, ComBat, PhosR-equivalent RUV/SPS/RUV-III, or limma
`removeBatchEffect` parity claims.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
