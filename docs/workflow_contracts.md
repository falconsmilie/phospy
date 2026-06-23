# Workflow Contracts

This page is the short contract map for PhosPy workflows. It keeps the public
boundaries visible and links to the detailed API pages that own request,
configuration, result, example, provenance, and limitation details.

For current support status, use
[Scientific Coverage](scientific-coverage.md). For fixture-backed comparison
evidence, use [Parity](parity.md). Roadmap ADRs are not current support claims.

## Public Workflow Map

| Area | Detailed API page | Public entrypoint | Request | Result |
| --- | --- | --- | --- | --- |
| Dataset building | [Dataset Build API](api/dataset-build-workflow.md) | `AnalysisReadyDatasetBuilder.run(...)` | `DatasetBuildRequest` | `AnalysisReadyPhosphoDataset` |
| Differential analysis | [Differential Analysis Workflow](api/differential-analysis.md) | `DifferentialAnalysisWorkflow.run(...)` | `DifferentialAnalysisRequest` | `DifferentialAnalysisResult` |
| Enrichment | [Enrichment Workflow](api/enrichment.md) | `EnrichmentWorkflow.run(...)` | `EnrichmentWorkflowRequest` | `EnrichmentWorkflowResult` |
| Kinase | [Kinase Workflow](api/kinase.md) | `KinaseWorkflow.run(...)` | `KinaseWorkflowRequest` | `KinaseWorkflowResult` |
| Signalome | [Signalome Workflow](api/signalome.md) | `SignalomeWorkflow.run(...)` | `SignalomeWorkflowRequest` | `SignalomeWorkflowResult` |

Public workflow classes expose `run(...)` only. User-facing docs should focus on
workflow classes, request objects, config objects, and result objects.

## Shared Boundary Rules

`AnalysisReadyPhosphoDataset` is the strict analysis-ready dataset boundary for
downstream dataset, differential, kinase, and signalome work.

Key identity rules:

- `dataset.phospho.index` and `dataset.site_metadata.index` are `site_key`.
- `display_id` is a human-readable label, not row identity.
- Duplicate `display_id` values remain valid when the corresponding `site_key`
  values differ.
- Duplicate rows that resolve to the same `site_key` are a scientific ambiguity
  and fail by default unless an explicit non-error duplicate-site preprocessing
  policy is chosen.
- Site-level workflow outputs that materialize identity include both `site_key`
  and `display_id` where applicable.

Requests are command payloads. Constructing a request records intent; the
builder or workflow run validates scientific compatibility before execution.

Result models are typed containers. Public helpers such as `to_dataframe()`,
`*_dataframe()`, `table`, `result_table`, and `to_payload()` return defensive
in-memory snapshots for inspection or handoff only. They are not exporters,
formatters, plotting helpers, report generators, or places to run additional
scientific post-processing.

## Dataset Builder Contract

The builder validates input tables, applies explicit preprocessing policy, and
returns a complete `AnalysisReadyPhosphoDataset`.

Important user-facing assumptions:

- `phospho` is a numeric site-by-sample matrix.
- `site_metadata` aligns to `phospho.index` at ingestion.
- `sample_metadata`, when supplied, is passive metadata. It does not infer
  differential conditions, replicates, batches, or blocks.
- Missing values are forbidden by default or handled by an explicit
  preprocessing policy.
- missing-data handling runs before normalisation in preprocessing stage order.
- row-median imputation is deterministic.
- row-median imputation is not left-censored imputation.
- Localisation should be configured before site-level scientific workflows when
  localisation confidence matters.
- `protein_id` is optional at the base dataset boundary, but signalome requires
  complete non-empty values for interpreted sites.

`DatasetPreprocessingConfig` owns transforms, normalisation, missing-data
handling, site construction, site-sequence resolution, total-protein correction,
protein-aware preparation, optional `linear_residualize_batch`, comparison
building, localisation policy, and `ruv_readiness` reporting.

`linear_residualize_batch` fixed-effect residualisation preserves condition
effects by design and requires explicit batch and condition metadata.
Confounded batch/condition designs are rejected. It is not ComBat, not RUV, not
limma `removeBatchEffect` parity, and not mixed-effects modelling. It does not
solve all batch-effect problems.

Do not interpret `ruv_readiness` as RUV support. It is metadata readiness
reporting only; it does not select SPS controls, run RUV/SPS/RUV-III
correction, or produce PhosR-equivalent corrected output.

`DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")` is
preparation-only. It does not modify phosphosite values, does not subtract total
protein, does not run joint PTM/protein differential modelling, does not adjust
differential models, and does not claim MSstatsPTM-style inference. Current
`DifferentialAnalysisWorkflow` does not consume `ProteinAwarePreparationResult`.

## Differential Analysis Contract

`DifferentialAnalysisWorkflow` consumes an `AnalysisReadyPhosphoDataset`, an
explicit `ExperimentalDesign`, explicit `Contrast` objects, and
`DifferentialAnalysisConfig`.

Important user-facing assumptions:

- Conditions, batches, blocks, and covariates are not inferred from sample
  names.
- Batch can be modelled as an ordinary fixed-effect covariate when the design is
  valid, full rank, and contrasts are estimable. This is not batch correction.
- Fixed-block paired designs are supported only when
  `paired_design_policy="fixed_block"` and every block has complete coverage for
  each requested contrast.
- Block terms are fixed effects, not random effects. This is not limma
  `duplicateCorrelation`, not mixed-effects modelling, and no mixed effects are
  fitted.
- Incomplete or partially covered blocks are rejected; PhosPy does not drop
  those blocks or samples.
- Upstream-imputed datasets are rejected by default. The explicit
  `withhold_imputed_features` policy uses dataset-owned imputation observation
  metadata and excludes withheld rows from the Benjamini-Hochberg denominator.

Each contrast result table is indexed by the input `site_key` and includes the
minimum public identity columns documented in the
[Differential Analysis Workflow](api/differential-analysis.md).

## Enrichment Contract

`EnrichmentWorkflow` runs offline over-representation analysis against
caller-supplied `GeneSetCollection`, `PtmSetCollection`, or homogeneous
`EnrichmentSetCollection` inputs.

Important user-facing assumptions:

- Provide exactly one selected-identifier source: `selected_identifiers` or
  `input_table` plus `identifier_column`.
- `identifier_kind`, selected identifiers, set members, and
  `background_universe` must use the same namespace.
- The background universe is explicit and required.
- Gene-level and site-level enrichment require explicit identifier semantics and
  are not interchangeable.
- `EnrichmentWorkflow` consumes selected foreground identifiers for ORA, not a
  complete ranked list. Future ranked-list enrichment is deferred by
  [ADR-0030](adr/adr_0030_ranked_list_enrichment_prerequisites.md).
- Enrichment ratio is a descriptive overlap summary, not a pathway activity
  score.
- Adjusted p-values describe statistical evidence under the ORA model and the
  selected correction method.
- ORA does not prove pathway activation, regulation, biological causality, or a
  mechanism.
- GO, KEGG, Reactome, PTM-SEA, and PTMsigDB resources are not bundled by this
  feature unless the caller supplies them as ordinary local collections.
- Enrichr, gseapy, clusterProfiler, and similar online calls are not native core
  workflow behavior.
- ORA is not GSEA, ssGSEA, or PTM-SEA support.

## Kinase Contract

`KinaseWorkflow` consumes an analysis-ready dataset and either a
`ReferencePreset` or explicit `ReferenceBundle`.

Important user-facing assumptions:

- Quantified rows are keyed by `site_key`; reference display IDs are projected
  through dataset `display_id` metadata.
- `reference_display_ambiguity_policy="error"` rejects one-display-to-many
  `site_key` projection by default.
- Site sequences are required for scoring rows.
- Scores are relative support values within a run, not calibrated
  probabilities or proof of causal regulation.
- Kinase Library-style workflow modes require a compatible caller-supplied
  local `KinaseLibraryResource`. PhosPy does not bundle official Kinase Library
  data and does not claim validated Kinase Library parity.
- Activity score output is optional and method-specific.
  `activity_result.activity_matrix` is the preferred method-neutral kinase
  activity score matrix.
- Activity scores depend on substrate coverage and reference evidence; sparse
  support weakens interpretation.
- KSEA-style activity scores are not equivalent to PhosR kinase activity
  inference, and ssGSEA-style activity-like scores are not PTM-SEA support.
- Causal kinase activity claims require external validation and study design
  support.

## Signalome Contract

`SignalomeWorkflow` consumes a `KinaseWorkflowResult`.

Important user-facing assumptions:

- The upstream kinase result must provide aligned prediction and downstream
  score matrices.
- Signalome aligns by `site_key`; it does not reinterpret display labels as row
  identity.
- Complete non-empty `dataset.site_metadata.protein_id` values are required for
  interpreted sites as signalome-specific protein grouping metadata.
- Module and network outputs are derived summaries, not causal proof.
- `candidate_scoring_policy="sampled"` approximates candidate module-count
  scoring only. It is not a general bypass for all signalome scale limits.
- Mixed corrected/uncorrected total-protein quantitative meaning is rejected by
  default unless explicitly allowed.

## Provenance and Reproducibility

Workflow provenance records resolved request/config choices, table
fingerprints, relevant scientific policy records, diagnostics, and environment
metadata where available. Environment provenance supports reproducibility
audits, but it does not guarantee bitwise identical numeric outputs across
different machines or dependency builds.

Saved workflow bundles and provenance payloads are supported only for the
current PhosPy schema. Legacy saved-result compatibility has been intentionally
removed. Regenerate older development-version outputs instead of relying on
schema repair during loading. Provenance remains supported for current outputs,
and current table-hash semantics are unchanged.

## Where Details Live

- [API Guide](api/guide.md) links to supported imports and workflow pages.
- [Dataset Build API](api/dataset-build-workflow.md) owns dataset request,
  preprocessing configuration details, workflow-specific request,
  config, result, example, provenance, and limitation details.
- [Scientific Coverage](scientific-coverage.md) is the maintained scope and
  coverage matrix.
- [Parity](parity.md) tracks fixture comparison evidence.
- [Performance Contracts](performance.md) covers scale limits and guardrails.
- [ADR Index](adr/index.md) stores maintainer decision records.
