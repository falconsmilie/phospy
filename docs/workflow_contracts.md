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
- `site_sequence` is required at the analysis-ready dataset boundary as
  plausible sequence evidence, but sequence-aware workflows may require stricter
  centered window contracts before execution.
- `protein_id` is optional at the base dataset boundary, but signalome requires
  complete non-empty values for interpreted sites.

`DatasetPreprocessingConfig` owns transforms, normalisation, missing-data
handling, site construction, site-sequence resolution, total-protein
correction, protein-aware preparation, comparison building, and localisation
policy.

Batch-correction preprocessing surfaces are explicit: optional
`linear_residualize_batch` fixed-effect residualisation and native SPS/RUV-style
correction through `SpsRuvBatchCorrectionConfig`.

`ruv_readiness` diagnostics are report-only metadata readiness reporting.

`linear_residualize_batch` fixed-effect residualisation preserves condition
effects by design and requires explicit batch and condition metadata.
Confounded batch/condition designs are rejected. It is not ComBat, not RUV, not
limma `removeBatchEffect` parity, and not mixed-effects modelling. It does not
solve all batch-effect problems.

Do not interpret `ruv_readiness` as RUV support. It is report-only metadata
readiness reporting; it does not select SPS controls, run correction, or produce
PhosR-equivalent corrected output. Native SPS/RUV-style correction is
executable only through the separate explicit `SpsRuvBatchCorrectionConfig`
preprocessing config.

Native SPS/RUV-style correction also stays in dataset preprocessing. The native
PhosPy SPS/RUV-style preprocessing correction estimates unwanted factors from
eligible caller-supplied control `site_key` residuals after protected-design
handling, applies the correction to the phosphosite matrix before downstream
workflows consume it, and records diagnostics plus `BatchCorrectionProvenance`.
Batch terms are resolved for validation and diagnostics, including
batch-associated-variance summaries; they are not directly residualized as
fixed effects by the native correction. Required inputs are aligned
`sample_metadata` with batch and protected condition columns, replicate
metadata only when the caller wants it recorded for provenance and diagnostics,
an explicit `ControlSiteSet`, `CorrectionMissingnessPolicy`,
`n_unwanted_factors`, diagnostics, and provenance. Multiple protected condition
columns mean the native correction protects joint condition strata such as
`condition=treated|timepoint=early`; it does not fit additive
protected-condition terms. Providing
`replicate_column` in the native lane does not enable replicate-aware RUV-III
correction semantics; replicate labels are not used for numerical
unwanted-factor estimation and do not enable RUV-III semantics. Temporary
imputation is correction mechanics only:
observation masks preserve which cells were originally observed, and imputed
temporary values must not be treated as observed evidence. Recognized
native-correction temporary-imputation policy/mechanics labels are `none` and
`row_median_temporary`; `row_median_temporary` is not public-workflow
permission to correct incomplete matrices. The public native workflow requires
a complete correction-stage matrix and rejects actual missing values (NaNs)
before executor invocation. `minprob_temporary` and `knn_temporary` are
rejected until supported semantics are implemented. Upstream-imputed cells
remain tracked through observation masks and are not treated as observed evidence.
Differential batch covariates
remain ordinary downstream model terms; they do not replace preprocessing
correction and are not removed from the differential design when the user
chooses to model them. The `ruv_iii_style` method label is not executable
unless a future feature implements replicate-aware RUV-III semantics.

Externally supplied `CorrectedPreprocessingOutput` must enter only at a safe
dataset preprocessing boundary. It cannot be combined with configured
downstream matrix-consuming preprocessing stages such as total-protein
correction, site-matrix construction, normalisation, or comparison building.
Provide external corrected output as the only matrix-changing preprocessing
input, or use native SPS/RUV-style `SpsRuvBatchCorrectionConfig` inside the
preprocessing pipeline when downstream preprocessing stages also need to run.

`ControlSiteSet` metadata is validated before native SPS/RUV-style correction.
Caller-supplied controls must provide auditable control-source metadata or an
explicit `metadata_missing_reason` entry for each unavailable caller-local
field. Audited fields include organism, identifier namespace, source identity,
source version, license, and redistribution. Formal or external source names
require source version. Packaged control references, if introduced later, must
provide complete organism, namespace, source name, source version, license, and
redistribution metadata. Validation rejects incomplete packaged metadata,
incompatible organism or namespace metadata when present, ambiguous accepted
control metadata, and caller controls missing audit metadata without rationale.

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
- Suspicious declared-only log2 intensity scale is rejected by default when
  declaration provenance records diagnostics. The explicit
  `allow_suspicious_declared_input_scale` differential override is recorded in
  policy provenance when used.

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
- Kinase Library-style motif scoring requires a fixed centered window matching
  the supplied `KinaseLibraryResource.sequence_window`. The selected sequence
  must have the expected length, center index, `S/T/Y` center residue, supported
  alphabet, accepted padding/lowercase/modified-symbol policy, and known source.
  Display IDs are not accepted as sequence context, and dataset/reference
  sequence conflicts must be resolved by the request conflict policy before
  scoring.
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
- Signalome enforces centered phosphosite sequence context for sequence-aware
  upstream identity, but it does not apply the fixed Kinase Library motif-window
  contract.
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
