# Scientific Coverage

PhosPy implements selected PhosR-style phosphoproteomics workflows.
"PhosR-style" and "PhosR-inspired" in this repository always mean
feature-scoped, evidence-scoped comparison lanes. They do not mean full PhosR
package equivalence.

## Current Supported Public Lanes

PhosPy currently supports these public analysis lanes:

1. build an `AnalysisReadyPhosphoDataset`
2. run differential phosphorylation analysis
3. run kinase scoring and prediction
4. optionally run signalome analysis from the kinase result
5. run offline over-representation enrichment against caller-supplied gene-set
   or PTM-set collections

Public input-preparation support additionally includes generic column-mapped
phosphosite import, MaxQuant phosphosite import, and
FragPipe/Philosopher/PTMProphet phosphosite import. Importers emit
`PhosphositeImportResult` candidates and dataset-builder requests; they do not
construct analysis-ready datasets, infer differential design, or bypass builder
validation. MaxQuant and FragPipe/PTMProphet fixture coverage is targeted to
selected localisation, mapping, flagging, grouped-row, intensity-column, and
ambiguity edge cases; it is not full vendor-tool parity.
Dedicated Spectronaut and DIA-NN semantic importers are not currently
supported. Manually mapping a compatible table through
`MappedPhosphositeTableImporter` remains generic mapped-table support, not
Spectronaut/DIA-NN support or upstream statistical result import.

The prediction science layer includes a pure Kinase Library-style motif scorer.
`KinaseWorkflow` exposes it through explicit scoring modes only; the default
kinase lane remains the PhosR-inspired rank-weighted scoring mode implemented
by PhosPy. It uses available substrate/reference evidence, profile support,
motif support when available, and rank-weighted combination under
minimum-substrate/support rules. It is not an exact PhosR implementation and is
not intended to provide numerical parity with PhosR. Kinase Library workflow
scoring is still a kinase workflow lane, not a fully independent official
Kinase Library implementation. It requires the normal workflow
reference context plus a caller-supplied `KinaseLibraryResource` with
compatible organism, residue-class lanes, score matrices, sequence-window
definition, and provenance. No official Kinase Library compatibility or parity
claim is made.

Sequence support uses a layered readiness contract. Analysis-ready sequence
evidence means `AnalysisReadyPhosphoDataset` carries a required, plausible
`site_sequence` value for each row. Workflow-specific sequence-context readiness
is validated later by sequence-aware workflows against the selected
dataset/reference sequence and active scoring contract. Every current kinase
scoring mode requires centered phosphosite context; Kinase Library
resource-backed motif modes additionally require the exact fixed centered
window, source, residue lane, and alphabet/padding policy supplied by the
caller-provided resource. Not every `site_sequence` value that passes dataset
construction is motif-ready.

Reference-context compatibility is conservative by default: kinase and signalome
workflows require known matching dataset/reference contexts. The explicit
`ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT` config override
permits unknown context only with result caveats; known mismatches never pass.

Native enrichment support is offline over-representation analysis (ORA) over
caller-supplied `GeneSetCollection`, `PtmSetCollection`, or homogeneous
`EnrichmentSetCollection` inputs. The background universe is explicit and
required; the caller supplies both the selected identifiers and background, and
PhosPy does not infer either from a dataset, reference bundle, or set
collection. Enrichment ratios are descriptive, and adjusted p-values describe
statistical evidence under the ORA model rather than pathway activation,
regulation, mechanism, or biological causality. PhosPy does not bundle GO,
KEGG, Reactome, PTM-SEA, or PTMsigDB resources for this feature, and the core
workflow does not call Enrichr, gseapy, clusterProfiler, or other online
services. Those online calls are not native core workflow behavior. ORA does
not imply GSEA, ssGSEA, or PTM-SEA support.

The shared multiple-testing helper adjusts only finite p-values for supported
correction methods. Its denominator is the number of finite p-values passed to
the helper, and non-finite input positions remain missing in the adjusted
output. This is the generic multiple-testing helper contract; workflow-specific
eligibility can be stricter.

Differential analysis requires analysis-ready numeric inputs plus valid
`ExperimentalDesign` and `Contrast` metadata. It does not infer design from
sample names and does not replace upstream preprocessing requirements.
`AnalysisReadyPhosphoDataset` requires a complete matrix, but complete does not
mean fully observed when upstream imputation filled cells. By default,
`DifferentialAnalysisWorkflow` rejects datasets where
`dataset.processing_state.missing_data.imputed` is true because current
differential model fitting does not treat imputed cells as fully observed
measurements. An explicit `imputed_value_policy="withhold_imputed_features"`
path is supported only when dataset-owned imputation observation metadata is
present. That policy marks features as `tested`,
`withheld_high_imputation`, or `withheld_insufficient_observed_values`, reports
per-feature `imputed_cell_count`, `observed_cell_count`, `imputed_fraction`, and
`imputation_policy`, and fits statistics only for tested rows. It is not
observed-only fitting and does not use feature-specific residual degrees of
freedom.
Differential designs may explicitly declare fixed-effect covariates on
`ExperimentalDesign`: batch, categorical covariates, and continuous covariates.
Modelled fixed-effect covariates are included in the fitted differential design
matrix and recorded in result provenance. These terms are not inferred from
passive sample metadata. A batch fixed effect is a model covariate; it is not
batch correction and is not ComBat, RUV, `removeBatchEffect`, limma
`duplicateCorrelation`, or mixed-effects modelling.

Dataset preprocessing KNN imputation is deterministic. PhosPy computes
`nan_euclidean` neighbour distances for `impute_knn` and breaks exact-distance
ties by stable row identity rather than input row order. Reproducibility-marked
unit tests cover repeated output, tie-heavy input, diagnostics, and provenance
for the configured `k`, distance, and row-filtering policy.

Separately, dataset preprocessing supports one opt-in residualisation method
under the `batch_correction` configuration group: `linear_residualize_batch`.
It is fixed-effect residualisation of batch terms. It preserves condition
effects by design: condition terms are included in the residualisation design.
Batch and condition metadata must be supplied explicitly through
`sample_metadata` and `DatasetBatchCorrectionConfig`. Confounded
batch/condition designs, rank-deficient designs, singleton batches, and designs
without residual degrees of freedom are rejected before correction. This support
is not ComBat, not RUV, not limma `removeBatchEffect` parity, and not
mixed-effects modelling.

### RUV/SPS/RUV-III Batch-Correction Posture

Supported today:

- `linear_residualize_batch`, a limited fixed-effect residualisation method
  under dataset preprocessing.
- `SpsRuvBatchCorrectionConfig`, a native SPS/RUV-style preprocessing
  correction lane that requires caller-supplied controls, protected design
  metadata, explicit control-source audit metadata or field-level
  missing-metadata rationale, explicit missingness policy, unwanted-factor
  count, diagnostics, and provenance. Recognized temporary-imputation
  policy/mechanics labels are `none` and `row_median_temporary`;
  `row_median_temporary` is not public-workflow permission to correct
  incomplete matrices. The public native workflow requires a complete
  correction-stage matrix and rejects actual missing values (NaNs) before
  executor invocation. `row_median_temporary` is a recognized policy/mechanics
  label, not permission for actual NaNs to pass through the public native
  workflow. Upstream-imputed cells remain tracked through
  observation-mask provenance and are not treated as observed evidence. It
  estimates unwanted factors from eligible control-site residuals after
  protected-design handling.
  Batch terms are resolved for validation and diagnostics, not directly
  residualized as fixed effects by the native correction. Supplied replicate
  metadata is validated and recorded, but it is not used for numerical
  unwanted-factor estimation and does not enable replicate-aware RUV-III
  semantics.
- `ruv_readiness`, when enabled, as diagnostic/report-only metadata readiness
  reporting.

Not supported today:

- no PhosR-equivalent SPS/RUV-III batch correction.
- no executable `ruv_iii_style` correction unless a future feature implements
  replicate-aware RUV-III semantics.
- no treatment of linear residualisation as equivalent to SPS/RUV/RUV-III
  correction.
- no treatment of linear residualisation as equivalent to PhosR-style RUV/SPS
  correction.
- no treatment of report-only `ruv_readiness` or similar diagnostics as
  correction support.

The native SPS/RUV-style lane is a PhosPy implementation. It is not a claim of
PhosR-equivalent SPS/RUV-III correction.

At dataset-construction boundary, PhosPy uses a protein-scoped analysis-ready
row key (`site_key`) and keeps `display_id` (for example `GENE;SITE;`) as a
human-readable label. `site_key` is required to be unique, while `display_id`
may repeat once `site_key` is the row identity. Direct analysis-ready datasets
must use `site_key` indexes and include auditable protein context metadata
(`organism`, `protein_namespace`, `protein_identifier`, and `site`); they must
not silently fall back to display-site identity. Builder input may accept legacy
display-indexed shape only when enough protein context exists to derive
`site_key`. Workflows operate on `site_key`, and site-level outputs that
materialize row identity include both `site_key` and `display_id`. `display_id`
is a human-readable label and may repeat. Rows that resolve to the same
`site_key` are a scientific ambiguity; the default duplicate-site policy fails,
and non-error policies are deliberate preprocessing choices that change which
evidence enters the analysis-ready dataset. Differential result tables and
direct public `DifferentialAnalysisResult` construction require encoded
protein-scoped `site_key` indexes plus explicit `site_key`, `display_id`,
`organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`
columns. Workflow-created differential results preserve that required protein
context and optional protein metadata such as `protein_id` when present.
Validation fails rather than repairing display-indexed, display-keyed,
arbitrary-keyed, or stat-only public result tables. Kinase references may use
display IDs only through the explicit
reference-to-dataset mapping layer; references remain reference/display
identifiers. Analysis-ready datasets may carry `site_metadata.protein_id` as
optional metadata, and it may be absent or incomplete at the base dataset
boundary. Signalome requires complete `site_metadata.protein_id` values as
algorithm-specific protein grouping metadata for protein-level module and
protein-site context summaries; that field is not core protein identity,
is not the dataset row identity, and is not encoded in `site_key`. See
[ADR-0024: Protein-Scoped Phosphosite Row Identity](adr/adr_0024_protein_scoped_phosphosite_row_identity.md).

### Differential Parity Envelope (Current Release)

Current differential analysis is not full PhosR or limma parity. Supported
designs are limited to tested design and contrast envelopes; broader limma or
PhosR behavior must not be inferred unless explicit parity tests and public docs
name that behavior.

`DifferentialAnalysisWorkflow` parity claims are currently scoped to:

- two-condition unpaired designs with biological-replicate rows
- simple condition-vs-condition contrasts with one `+1` and one `-1` term
- empirical-Bayes modes: `method="standard"` and `method="robust"` with
  optional `trend=True`
- Benjamini-Hochberg multiple-testing adjustment (`adj.P.Val`)

Validated PhosPy fixed-effect support additionally includes execution for:

- explicit batch-as-fixed-effect covariates
- declared categorical fixed-effect covariates
- declared continuous fixed-effect covariates
- explicit `paired_design_policy="fixed_block"` designs where every block has
  complete coverage for each requested condition contrast
- rank and contrast-estimability validation before execution

These fixed-effect covariate models are support for ordinary fixed terms in the
linear design. They are not batch-correction methods, correlated-replicate
modelling, or mixed-effects modelling.
Fixed-effect covariates are not full batch correction, and fixed-effect batch
terms are not substitutes for full batch correction or mixed-effect modelling.

Fixed-block paired designs are supported only when the caller explicitly sets
`paired_design_policy="fixed_block"` and supplies `SampleDesignRecord.block_id`
for every analysed sample. The block terms are ordinary fixed effects in the
design matrix; block terms are fixed effects, not random effects. This is not
limma `duplicateCorrelation`, not mixed-effects modelling, and not random
subject modelling; no mixed effects are fitted. For
every requested condition contrast, every block must contain both numerator and
denominator conditions.
Incomplete or partially covered blocks are rejected before execution; PhosPy
does not drop those blocks or samples to make the design fit. Simple unpaired
workflows keep the default `paired_design_policy="reject"` behavior and are not
changed by fixed-block support.

Explicitly unsupported in this release:

- limma `duplicateCorrelation`-style correlated-replicate modelling
- mixed-effect differential modelling

Contract difference vs limma/PhosR surface:

- analysis-ready inputs must be complete at boundary; missing values are
  rejected before differential execution instead of being handled inside
  differential model fitting.
- differential execution validates generated `P.Value` values are finite and
  within `[0, 1]` before Benjamini-Hochberg adjustment.
- upstream-imputed analysis-ready inputs are rejected by default. The explicit
  `withhold_imputed_features` policy can withhold rows above
  `imputed_value_max_fraction` or with insufficient originally observed samples
  in contrast conditions. Withheld rows are excluded from model fitting and
  Benjamini-Hochberg adjustment, and receive missing test statistics.

Bundled runtime references in the current release are rat-only. Human and mouse
analysis can be run by passing an explicit `ReferenceBundle` in Python. No
packaged human or mouse reference lane is committed in this release because no
approved redistributable source bundle with complete license and provenance
metadata is included. The packaged rat `l6_native` lane now has a
hash-verifiable reference manifest covering the logical bundle and all packaged
files. That manifest records a PhosR 1.20.0-derived lineage and typed
exact-file upstream-package license evidence for this exact PhosPy-packaged
snapshot. The
approval does not generalize to other upstream databases, future PhosR versions,
future PhosPy bundles, or arbitrary external datasets. Bundled release
eligibility is governed by `redistribution_status`, not by broad source-name
assumptions. The release checks reject bundled manifests that are missing
required metadata, omit packaged files, fail file hash verification, omit
typed exact-file redistribution evidence, or declare a
non-release-eligible status. Unresolved bundled references block release.
External-only references must remain caller-supplied and must not be shipped as
bundled data. Codex agents and human developers must not mark references
`approved` without verified evidence in the manifest for the exact files being
packaged.

## Roadmap Visibility and Guardrails

[ADR-0025: Competitive Phosphoproteomics Workflow Coverage Roadmap](adr/adr_0025_competitive_phosphoproteomics_workflow_coverage.md)
records intended future direction for references, kinase inference, importers,
richer differential designs, enrichment, visualisation, and possible CLI
support.

Roadmap entries are not current feature claims. A roadmap item becomes supported
only when executable implementation, public contracts, documentation, and tests
exist, and this page is updated to the correct scope category.

| Roadmap area | Current status | Direction, not current support |
| --- | --- | --- |
| References | Bundled runtime references are rat-only. The only approved packaged runtime reference is the exact rat `l6_native` snapshot with a hash-verifiable manifest and typed exact-file upstream PhosR 1.20.0 license evidence. Human and mouse workflows require an explicit caller-supplied `ReferenceBundle`. | Broader reference handling should use explicit provenance, compatibility checks, external bundle validation, and file-level hashes. New bundled data requires verified typed redistribution evidence, provenance, docs, and tests before `_BUNDLED_DEFAULTS` is updated. Packaged manifests with `unresolved` or `external_only` redistribution status block release. File hashes alone do not establish redistribution approval. |
| Kinase inference | Kinase scoring/prediction and three explicit activity methods are executable. Scores are relative support or substrate-set summaries, not calibrated causal inference. | Additional kinase inference or activity methods should be added one method at a time with stable scientific policy records and method-specific validation. |
| Importers | PhosPy supports analysis-ready tables, generic table I/O contracts, generic column-mapped phosphosite import, MaxQuant phosphosite import, and FragPipe/Philosopher/PTMProphet phosphosite import. These are input-preparation adapters that feed dataset-builder validation. They are not broad support for all vendor, search-engine, upstream statistical outputs, Spectronaut, or DIA-NN. | Additional semantic importers should produce typed tables or requests that still pass builder and workflow validation; they must not bypass site identity or provenance contracts. Spectronaut and DIA-NN remain future/demand-driven candidates, not current support. |
| Richer differential designs | Current parity-protected differential lane is two-condition unpaired simple contrasts. Fixed-effect batch, categorical covariate, continuous covariate, and explicit complete fixed-block terms are executable as ordinary fixed covariates with completeness, rank, and estimability validation. Correlated repeated-measure, `duplicateCorrelation`-style, and mixed-effect modelling are not executable in this release. | Additional richer designs require explicit design/result contracts, provenance, validation, and parity or method-specific evidence before any support claim. |
| Batch correction | Current executable preprocessing support includes `linear_residualize_batch` fixed-effect residualisation and native SPS/RUV-style correction through `SpsRuvBatchCorrectionConfig`; `ruv_readiness` remains diagnostic/report-only metadata readiness reporting. PhosR-equivalent SPS/RUV-III batch correction is not supported today. | Additional batch-correction methods require explicit control-site, design, imputation, diagnostics, provenance, downstream eligibility, validation, and test contracts. |
| Enrichment | Native support is offline ORA over caller-supplied gene-set or PTM-set collections with explicit identifier semantics and background. No online resource access is bundled or executed. | Additional enrichment methods or curated resource integrations require explicit data provenance, redistribution approval where relevant, identifier semantics, validation, docs, and tests before support is claimed. |
| Visualisation | Core PhosPy has no first-class visualisation workflow/API. | Visualisation should consume validated result objects and must not become a hidden analysis engine or source of scientific truth. |
| CLI workflow support | Scientific workflow execution through a CLI is not currently supported; the Python API is the supported interface. | Any future CLI must be a thin wrapper over Python API requests/workflows and satisfy ADR-0022 reintroduction criteria before support is claimed. |

## Scope Categories

Every scientific scope claim in public docs must map to one category:

| Category | Meaning |
| --- | --- |
| `parity-gated` | Executable lane with active fixture-backed parity checks in release checks |
| `validated PhosPy implementation` | Executable lane validated by PhosPy contract/unit/integration tests; not a PhosR-equivalence claim by itself |
| `experimental` | Executable but intentionally provisional/approximate behavior with explicit caveats |
| `open gap` | Not currently executable in the supported public workflow lane |
| `deliberate scope difference` | Intentionally different from PhosR surface or intentionally narrowed contract |
| `not planned` | Intentionally outside supported scope |

## Scientific Scope Matrix (Single Source of Truth)

<!-- Documentation smoke-test marker: ## Scientific Scope Matrix (Single Source Of Truth) -->

This matrix is the maintained user-facing scope source for PhosPy. Parity is
feature-specific and evidence-scoped. Full PhosR package equivalence is not
claimed.

| Area | Scope category | Current executable support | Evidence and release checks | Limits and non-claims |
| --- | --- | --- | --- | --- |
| Differential analysis | `parity-gated` | `DifferentialAnalysisWorkflow` for two-condition unpaired simple contrasts with empirical-Bayes `standard`/`robust` and optional `trend`; fixed-effect batch, categorical covariate, continuous covariate, and complete fixed-block terms are executable as ordinary design covariates. Upstream-imputed datasets remain rejected by default; `imputed_value_policy="withhold_imputed_features"` is an explicit validated PhosPy policy when imputation observation metadata is present. | `tests/parity/test_differential_analysis_parity.py`, `tests/parity/test_differential_limma_parity.py`, plus unit/integration design, fixed-effect provenance, result-contract tests, and `tests/unit/test_differential_imputation_policy.py` | Fixed-effect batch terms are not batch correction. Fixed-block terms require complete within-block contrast coverage and full-rank/estimable designs. Correlated repeated-measure, limma `duplicateCorrelation`-style, mixed-effect, and random subject-effect designs are rejected in this release. Missing values are rejected at analysis-ready boundary before model fitting. Imputed cells are not treated as fully observed by default. `withhold_imputed_features` withholds high-imputation or insufficient-observation rows, reports status, and excludes withheld rows from model fitting and the Benjamini-Hochberg denominator. It is not observed-only fitting. |
| Kinase scoring | `parity-gated` | `KinaseWorkflow` default `scoring_mode="phosr_rank_weighted"` PhosR-inspired profile/motif scoring and rank-weighted fusion implemented by PhosPy | `tests/parity/test_kinase_workflow_parity.py`, `tests/parity/test_prediction_science_parity.py`, `tests/parity/test_l6_prediction_parity.py` | Relative support scoring only; not calibrated causal inference. The mode is not an exact PhosR implementation and is not intended to provide numerical parity with PhosR. Kinase Library scoring is not the default parity lane. |
| Kinase Library motif scoring | `validated PhosPy implementation` | Pure science-layer `KinaseLibraryMotifScorer` / `score_kinase_library_motifs`, plus opt-in `KinaseWorkflow` modes `kinase_library_motif` and `combined_profile_motif` for supplied Kinase Library-style resources | `tests/unit/test_kinase_library_motif_scoring.py`, `tests/unit/test_kinase_library_workflow_requirements.py`, `tests/science/test_kinase_library_motif_scoring_science.py`, `tests/integration/test_kinase_library_workflow_scoring.py` | Workflow mode still requires resolved `ReferenceBundle` context with `kinase_substrate_map` overlap, eligible kinases at `min_substrates`, and resolved site sequences. It also requires an explicit compatible local `KinaseLibraryResource`. Missing matching residue-class lanes fail validation; they do not fall back to PhosR-inspired motif scoring. Workflow motif scores are normalized to unit interval per kinase matrix for within-run ranking support; raw science-layer motif scores preserve provider scale. Scores are not probabilities. Ser/Thr and Tyr matrix lanes are not interchangeable. No official Kinase Library parity claim is made. |
| Kinase prediction | `parity-gated` | Deterministic and adaptive kinase prediction in `KinaseWorkflow` | `tests/parity/test_public_predmat_parity.py`, `tests/parity/test_l6_prediction_parity.py`, `tests/parity/test_adaptive_prediction_parity.py`, `tests/parity/test_adaptive_replay_parity.py` | Prediction scores are ranking support, not probabilities. |
| Kinase activity scoring | `validated PhosPy implementation` | Supported exploratory activity score methods: `simplified_weighted_substrate_activity_v1`, `ksea_zscore_activity_v1`, and `ssgsea_substrate_enrichment_activity_v1` | Unit activity tests (`tests/unit/test_activity_science.py`), workflow activity tests (`tests/workflows/test_kinase_activity_ssgsea.py`), and parity activity gate (`tests/parity/test_activity_stage_parity.py`) | Scores depend on substrate coverage and reference evidence; sparse support weakens interpretation. KSEA-style activity scores are not a claim of full PhosR kinase activity equivalence. ssGSEA-style activity-like scores are a PhosPy rank-walk implementation and are not a PTM-SEA parity claim. Causal kinase activity claims require external validation. |
| Signalome analysis | `parity-gated` | `SignalomeWorkflow` module assignment, network outputs, and protein-site context | `tests/parity/test_signalome_workflow_parity.py`, `tests/parity/test_signalome_clustering_backend_parity.py` | Derived summaries, not causal proof. Requires explicit signalome protein grouping metadata in `site_metadata.protein_id`. |
| Signalome sampled candidate scoring policy | `experimental` | `SignalomeConfig.sampled_candidate_scoring()` approximates candidate module-count scoring | Parity/contract coverage through signalome parity tests and workflow contract checks | Approximation applies to candidate scoring only; tree generation remains exact-policy governed. |
| Sequence context | `parity-gated` | `site_sequence` is required as analysis-ready sequence evidence, and kinase/signalome validators separately enforce workflow-specific centered sequence-context readiness before sequence-aware scoring/prediction. | `tests/parity/test_l6_prediction_parity.py`, `tests/parity/test_prediction_science_parity.py` | Base dataset validation is plausibility-level and does not make every `site_sequence` motif-ready. Kinase Library-style motif scoring requires the exact centered window from the selected `KinaseLibraryResource`, with known source and compatible residue lane. |
| Localisation handling | `validated PhosPy implementation` | Localisation confidence validation and fail-fast threshold policies are supported | `tests/unit/test_localisation_policy_preprocessing.py`, `tests/unit/test_validator_boundaries.py` | No full localisation-filter workflow parity claim in this release. |
| Phosphosite importers | `validated PhosPy implementation` | Generic `MappedPhosphositeTableImporter`, `MaxQuantPhosphositeImporter`, and `FragPipePTMProphetImporter` translate upstream phosphosite tables into `PhosphositeImportResult` candidates and dataset-builder requests | `tests/unit/test_maxquant_phosphosite_importer.py`, `tests/unit/test_fragpipe_ptmprophet_importer.py`, `tests/integration/test_maxquant_importer_dataset_integration.py`, `tests/integration/test_fragpipe_importer_dataset_integration.py` | Importers do not construct analysis-ready datasets, infer sample groups, infer contrasts, infer batches or blocks, infer differential design, or bypass builder validation. Targeted MaxQuant and FragPipe/PTMProphet adapters are not broad vendor/search-engine parity, Spectronaut/DIA-NN support, or upstream statistical result import. Generic mapped import of a caller-mapped compatible table is not a dedicated Spectronaut or DIA-NN importer. |
| Missing values | `parity-gated` | Missing-data policy execution in preprocessing and downstream score preconditioning | `tests/parity/test_preprocessing_science_parity.py`, unit missing-data tests | Policy choice changes retained rows and downstream behavior. |
| Imputation | `validated PhosPy implementation` | Supported preprocessing policies include `row_median`, `minprob`, `knn`; imputed datasets expose typed per-feature observation metadata (`imputed_cell_count`, `observed_cell_count`, `imputed_fraction`) and an aligned defensive observed-cell mask export. Differential analysis can consume that metadata only through explicit `imputed_value_policy="withhold_imputed_features"`. | Unit preprocessing/scientific invariant tests plus `tests/unit/test_differential_imputation_policy.py` | Policy-dependent behavior; not blanket PhosR-equivalent imputation. Differential analysis still rejects upstream-imputed datasets by default. The withhold policy reports deterministic row statuses and excludes withheld rows from testing; it does not implement observed-only fitting. |
| Normalisation | `parity-gated` | Supported methods: `none`, `median_center`, `quantile` with stage-order provenance | `tests/parity/test_preprocessing_science_parity.py`, unit preprocessing tests | Method-specific claims only; no blanket normalisation equivalence claim. |
| Total-protein subtraction: `subtract_log_total` | `validated PhosPy implementation` | Dataset preprocessing can subtract matched log-scale total-protein abundance from log-scale phosphosite abundance (`log2_phospho - log2_total`) | Unit total-protein correction tests and dataset integration tests | Direct transformation only. Requires total-protein input and compatible log2 preprocessing. Not protein-aware modelling, not normalisation, not joint PTM/protein inference, and not MSstatsPTM equivalence. |
| Protein-aware preparation | `validated PhosPy implementation` | `DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")` builds `ProteinAwarePreparationResult` and `ProteinAwarePreparationReport` with matched phosphosite/protein pairs, sample-aligned protein covariates, eligibility rows, mapping diagnostics, sample-alignment diagnostics, transformation-state diagnostics, and explicit limitations | Unit protein-aware preparation/mapping/sample-alignment tests and dataset integration diagnostics tests | Preparation-only model-input preparation. Does not modify phosphosite values, subtract total protein, normalise intensities, run joint PTM/protein modelling, adjust differential models, or claim MSstatsPTM-style inference or equivalence. Current `DifferentialAnalysisWorkflow` does not consume the prepared covariate matrix. |
| Batch correction: `linear_residualize_batch` | `validated PhosPy implementation` | Dataset preprocessing supports opt-in `linear_residualize_batch` fixed-effect residualisation with explicit batch/condition metadata, condition-preserving design, design adequacy validation, and typed reports | Unit batch-correction engine, metadata, validation, and report tests; dataset integration tests for no-op, applied correction, invalid metadata/designs, and downstream differential consumption | This is the limited fixed-effect residualisation method under the `batch_correction` config group. Confounded batch/condition designs are rejected because preserving condition effects would otherwise be impossible. This is not ComBat, not RUV, not native SPS/RUV-style correction, not limma `removeBatchEffect` parity, not limma `duplicateCorrelation`, and not mixed-effects modelling. It does not solve all batch-effect problems. Differential batch fixed effects are model covariates, not preprocessing correction. |
| Batch correction: `SpsRuvBatchCorrectionConfig` | `validated PhosPy implementation` | Dataset preprocessing supports native SPS/RUV-style correction only through an explicit structured config with caller-supplied controls, batch column, protected condition terms, optional replicate metadata for provenance and diagnostics only, missingness policy, unwanted-factor count, diagnostics, provenance, and the supported stage placement after missing-data handling and before downstream preprocessing consumers. Supplied replicate metadata is validated and rejected when labels are all the same, all unique, perfectly confounded with batch, or perfectly confounded with protected condition metadata. Replicate metadata is not used for numerical unwanted-factor estimation and does not enable RUV-III or replicate-aware RUV-III semantics. The public native workflow requires a complete correction-stage matrix and rejects actual missing values (NaNs) before executor invocation. Recognized temporary-imputation policy/mechanics labels are `none` and `row_median_temporary`; `row_median_temporary` is not public-workflow permission to correct incomplete matrices and is a recognized policy/mechanics label, not permission for actual NaNs to pass through the public native workflow. Upstream-imputed cells remain tracked through observation-mask provenance and are not treated as observed evidence. | Unit config/validator tests, SPS/RUV-style executor tests, workflow orchestration/provenance tests, and dataset integration tests for public preprocessing execution | This is not PhosR-equivalent SPS/RUV-III parity, not executable RUV-III support, not replicate-aware RUV-III semantics, not ComBat, not limma `removeBatchEffect`, and not a hidden control-selection or online-control lookup feature. Controls must be explicit `site_key` annotations with control-source metadata or field-level missing-metadata rationale: caller-supplied controls audit organism, identifier namespace, source identity, source version, license, and redistribution; packaged controls, if added, require complete organism, namespace, source, version, license, and redistribution metadata. Correction remains in dataset preprocessing, not downstream workflows. Unsupported `stage_order` values are rejected rather than recorded as if they ran. `minprob_temporary` and `knn_temporary` are rejected for native correction until supported semantics are implemented. |
| Joint PTM/protein modelling and MSstatsPTM-style inference | `open gap` | No executable joint phosphosite/total-protein differential modelling lane is supported | N/A | Do not interpret total-protein subtraction or protein-aware preparation as MSstatsPTM-style inference, protein-adjusted differential modelling, or equivalence to MSstatsPTM. |
| PhosR-equivalent SPS/RUV-III, ComBat, and `removeBatchEffect` parity | `open gap` | No PhosR-equivalent SPS/RUV-III, ComBat, or limma `removeBatchEffect` parity lane is supported. `ruv_readiness` is report-only RUV-readiness metadata. | Readiness/report tests only; no correction parity claim | Do not interpret `ruv_readiness` as RUV support. Do not interpret `linear_residualize_batch` fixed-effect residualisation as ComBat, RUV, SPS/RUV-III, limma `removeBatchEffect`, or mixed-effects support. Do not interpret native SPS/RUV-style correction as PhosR parity. |
| Enrichment | `validated PhosPy implementation` | `EnrichmentWorkflow` runs offline over-representation analysis over caller-supplied `GeneSetCollection`, `PtmSetCollection`, or homogeneous `EnrichmentSetCollection` inputs with explicit identifier kind, selected identifiers, background universe, and multiple-testing correction | `tests/unit/test_enrichment_ora.py`, `tests/unit/test_enrichment_workflow_validation.py`, `tests/unit/test_public_contract_enrichment.py`, `tests/workflows/test_enrichment_workflow.py` | Requires user-supplied collections and explicit background. GO, KEGG, Reactome, PTM-SEA, and PTMsigDB resources are not bundled unless supplied as ordinary local collections. Online Enrichr, gseapy, clusterProfiler, and similar calls are not native core workflow behavior. ORA is not GSEA, ssGSEA, or PTM-SEA support. Gene-level and site-level enrichment require explicit identifier semantics. |
| Visualisation | `deliberate scope difference` | No first-class visualization workflow/API in core PhosPy | N/A | Visualization is intentionally out of current scientific parity scope. |
| Supported bundled organisms and references | `deliberate scope difference` | Bundled runtime references are rat-only for `ReferencePreset.AUTO` in this release. The exact rat `l6_native` snapshot is derived from upstream PhosR 1.20.0 package data. | Runtime behavior, reference compatibility tests, manifest `redistribution_status` release checks, typed exact-file rat manifest evidence, and workflow docs | Human/mouse are valid organisms but require explicit caller-supplied `ReferenceBundle` unless a future release commits approved redistributable packaged data with verified typed evidence. External-only references must not be shipped as bundled data, unresolved bundled references block release, and the rat `l6_native` bundle should not be treated as redistribution approval for other reference data, future bundles, other rat bundles, other organisms, or an independent direct permission claim from PhosphoSitePlus, PRIDE, Kinase Library, or another upstream scientific database. |
| Full PhosR package equivalence claim | `not planned` | Not claimed | Guardrail documentation in this matrix and `docs/parity.md` | Any implication of global PhosR parity is out of scope. |

## Release-Check Scientific Checks

Release-bearing scientific checks are documented and executable through the
maintainer command, `make release-check`. Default `pytest` is not sufficient for
publishing because it does not run threshold-bearing parity lanes, performance
contracts, checked-in reference validation, distribution metadata checks, or
packaged-reference validation. The maintained commands/workflows are:

- Local release command: `make release-check`
- `make release-check` executes lint, type checking, `pytest -m "not parity"`,
  `pytest tests/parity -m parity -s`,
  `pytest tests/performance -m "performance or release_gate"`,
  `python scripts/validate_reference_bundle_index.py --repo-root .`, and
  `make build`.
- `make build` starts from an empty `dist/`, builds one wheel and one sdist,
  runs metadata checks, and validates packaged reference manifests and declared
  file hashes in both archives. It does not require Git metadata.
- The publish workflow runs `make release-check` once on the checked-out tag and
  publishes the freshly built `dist/` artifacts through trusted publishing.
- CI still runs clean constrained `[dev,test]` installs and the full default
  source suite on Python 3.10, 3.11, and 3.12.
- This process provides normal CI/build confidence, not formal
  exact-source/exact-artifact attestation.
- CI parity workflows:
  - `.github/workflows/ci.yml` job `activity-parity-gate` runs `pytest tests/parity/test_activity_stage_parity.py -m "parity and activity_parity" -s`
  - `.github/workflows/ci.yml` job `parity-tests` runs `pytest tests/parity -m "parity and not parity_diagnostic" -s`
  - `.github/workflows/ci.yml` job `parity-diagnostics` runs `pytest tests/parity -m "parity_diagnostic" -s` with `continue-on-error: true`

## Interpretation Limits

- Weighted activity-like score output
  (`simplified_weighted_substrate_activity_v1`) is a heuristic summary over
  predicted substrates above threshold/top-N support.
- KSEA-style kinase activity score output (`ksea_zscore_activity_v1`) applies
  unweighted substrate-set enrichment z-scores after evidence thresholding and
  reports p-values (and q-values when enabled).
- ssGSEA-style activity-like score output
  (`ssgsea_substrate_enrichment_activity_v1`) applies a deterministic rank-walk
  enrichment score over phosphosite effect values using explicit
  kinase-substrate membership and reports seeded empirical permutation p-values
  when requested.
- Activity scores depend on substrate coverage and reference evidence; missing
  or sparse substrate support weakens interpretation.
- Causal kinase activity claims require external validation and study design
  support.
- KSEA-style activity scores are not equivalent to full PhosR kinase activity
  inference.
- ssGSEA-style activity-like scores are not a PTM-SEA parity claim.
- Rank-weighted fusion scores in `scoring_mode="phosr_rank_weighted"` are
  PhosR-inspired PhosPy scores. They combine profile-correlation and
  motif-frequency evidence using rank-derived weights, and they are not an
  exact PhosR implementation or numerical compatibility mode.
- Kinase Library-style science-layer motif scores are raw position-specific
  matrix sums on the caller-supplied score scale. Optional percentiles/ranks are
  empirical summaries against caller-supplied reference distributions only.
- Kinase Library workflow motif scores are normalized per kinase matrix to a
  unit interval for within-run ranking support. They are not calibrated
  probabilities and do not imply activity without an explicit activity method.
- An `AnalysisReadyPhosphoDataset` with `site_sequence` has required sequence
  evidence, not guaranteed motif-ready sequence context. Sequence-aware
  workflows may still reject rows or requests that fail centered-context,
  sequence-source, residue, alphabet, padding, or conflict-policy requirements.
- In `KinaseWorkflow`, `scoring_mode="kinase_library_motif"` still runs inside
  normal kinase workflow orchestration. Reference discovery and display-to-row
  projection happen before scoring; profile context from the resolved
  kinase-substrate map is still required to establish eligible kinases. The
  Kinase Library-style resource supplies motif matrices for the authoritative
  workflow score matrix, not substrate-map membership or activity inference.
- Signalome module/network scores are derived summaries, not probabilities,
  calibrated confidence values, or causal proof.
- Missing kinase correlations stay missing. `0.0` means a finite near-zero
  correlation was estimated.
- Differential phosphorylation results depend on valid design matrices, contrast
  definitions, replicate structure, and upstream preprocessing quality.
- Differential analysis does not resolve peptide/site ambiguity, localisation
  confidence, imputation, normalisation, or broad batch-effect correction.
  Batch-related preprocessing is upstream dataset-builder work: either opt-in
  `linear_residualize_batch` fixed-effect residualisation or native
  SPS/RUV-style correction through `SpsRuvBatchCorrectionConfig`.
- Differential analysis rejects upstream-imputed datasets by default. The
  explicit `withhold_imputed_features` policy consumes dataset-owned imputation
  observation metadata to mark `tested`, `withheld_high_imputation`, and
  `withheld_insufficient_observed_values` rows. Withheld rows receive missing
  `logFC`, `t`, `P.Value`, and `adj.P.Val`, and are excluded from the
  Benjamini-Hochberg denominator. This policy does not implement observed-only
  fitting or feature-specific residual degrees of freedom.
- Shared multiple-testing adjustment ranks or counts finite p-values only; the
  denominator is the finite p-value count passed to the helper. Non-finite
  positions remain missing in adjusted output. Differential workflow output is
  stricter for tested rows: generated `P.Value` values must be finite and in
  `[0, 1]` before BH is called, while imputation-withheld rows are excluded
  before adjustment.
- Dataset preprocessing `linear_residualize_batch` is opt-in fixed-effect
  residualisation. It preserves condition effects by including condition terms
  in the residualisation design, rejects confounded batch/condition designs, and
  records a typed `BatchCorrectionReport`. It is not ComBat, RUV,
  `removeBatchEffect` parity, `duplicateCorrelation`, or mixed-effects
  modelling, and it does not solve all batch-effect problems.
- Native SPS/RUV-style batch correction is exposed only through
  `SpsRuvBatchCorrectionConfig` with explicit controls, protected design terms,
  missingness policy, factor count, diagnostics, and provenance. It is not
  PhosR-equivalent SPS/RUV-III parity. `ruv_readiness` diagnostics are
  report-only metadata readiness signals and must not be interpreted as
  correction support.
- Native SPS/RUV-style support is preprocessing-owned. Differential, kinase,
  enrichment, and signalome workflows must not own correction logic.
- Total-protein correction `subtract_log_total` is a direct
  `log2_phospho - log2_total` transformation for matched rows. It is not
  protein-aware differential modelling, not normalisation, and not
  MSstatsPTM-style inference.
- Protein-aware preparation is preparation-only model-input preparation. It
  produces aligned phosphosite/protein input contracts and diagnostics only. It
  records mapping, missing-total, sample-alignment, transformation-state, and
  limitation fields. It does not modify phosphosite values, does not subtract
  total protein, does not run joint PTM/protein differential modelling, does
  not adjust differential models, and does not claim MSstatsPTM-style inference
  or equivalence.
- Enrichment workflow support is offline ORA against supplied gene-set or
  PTM-set collections. It uses the caller's selected identifiers, explicit
  background universe, and identifier semantics.
- Enrichment ratios are descriptive overlap summaries, and adjusted p-values
  describe statistical evidence under the ORA model and selected correction
  method. They do not prove pathway activation, regulation, mechanism, or
  biological causality.
- Enrichment workflow support does not bundle or fetch GO, KEGG, Reactome,
  PTM-SEA, or PTMsigDB resources. Online Enrichr, gseapy, clusterProfiler, and
  similar calls are not native core workflow behavior.
- Enrichment ORA results must not be interpreted as GSEA, ssGSEA, or PTM-SEA
  support.
- Joint PTM/protein differential modelling is not implemented in this release.
  `DifferentialAnalysisWorkflow` does not consume
  `ProteinAwarePreparationResult`.
- Fixed-effect covariates in differential analysis are ordinary fixed terms in
  the design matrix. Batch can be modelled this way, but this does not remove
  batch effects from data and does not implement ComBat, RUV,
  `removeBatchEffect`, `duplicateCorrelation`, correlated repeated-measure
  models, or mixed-effects models. Explicit fixed-block terms are supported
  only when every block is complete for the requested condition contrasts and
  the resolved design is full rank with estimable contrasts.
- Adjusted p-values control false discovery rate according to the implemented
  correction method; they do not validate biological causality.

## Scientific Policy Records

Workflow provenance includes machine-readable `scientific_policies` records.
Each record carries:

- stable policy ID
- name and version
- plain-language description
- active parameters
- scientific assumptions
- output scale/meaning

Ownership of scientific policy modules is domain-scoped:

- shared models: `phospy.provenance.scientific_policy_models`
- prediction: `phospy.science.prediction.scientific_policies`
- activities: `phospy.science.activities.scientific_policies`
- preprocessing: `phospy.science.datasets.preprocessing.scientific_policies`
- signalome workflow: `phospy.workflows.signalome.scientific_policies`
- signalome clustering: `phospy.science.signalomes.clustering.scientific_policies`
- differential aggregation: `phospy.science.differential.aggregation.scientific_policies`

Differential outputs now expose structured policy provenance through
`DifferentialAnalysisResult.policy_provenance`, including:

- design formula/description, condition columns, covariate columns, and
  covariate kinds
- paired-design policy, `block_id` field name, block count, included block
  levels, block column names, condition-coverage rule, and fixed-block
  limitations when applicable
- rank and contrast-estimability validation status
- explicit contrast definitions and contrast vectors
- replicate/group requirements and technical-replicate lineage
- empirical-Bayes moderation settings
- p-value and adjusted p-value methods
- missing-value handling policy, differential imputed-value policy, imputation
  fraction threshold, and adjusted-p-value scope
- unsupported-design rejection policy and intentionally rejected unsupported
  design features (`duplicateCorrelation`-style correlated-replicate and
  mixed-effect or random subject-effect modelling)

Protein-aware preparation records preparation-only model-input preparation
provenance without implying a downstream model. `ProteinAwarePreparationReport`
records the preparation mode, protein-mapping policy, eligibility counts,
sample-alignment and
transformation-state diagnostics, missing-total and ambiguous-mapping
diagnostics, and explicit limitations. Dataset run provenance also records the
active preparation summary in
`workflow_parameters["protein_aware_preparation"]` when preparation runs. These
records state that the stage does not modify phosphosite values, does not
subtract total protein, does not normalise intensities, does not run joint
PTM/protein differential modelling, does not adjust differential models, and
does not claim MSstatsPTM equivalence.

Enrichment workflow provenance records the ORA method, identifier column and
kind, collection kind, analysis level, explicit background universe size,
selected identifier count, selected identifier source, set collection
source/name/version metadata when available, multiple-testing correction,
offline/no-online-resource policy, and explicit limitations. These records
state what was executed; they do not imply bundled GO, KEGG, Reactome, PTM-SEA,
online enrichment service access, GSEA, ssGSEA, or PTM-SEA support.

When callers supply typed selected/background identifier-set provenance,
enrichment result provenance also records compact per-role source type, source
label, normalized identifier count, upstream workflow/result IDs, and typed
input-intensity-scale evidence. This provenance is optional for legacy/manual
identifier lists. It is required when a caller labels an identifier set as
PhosPy-derived quantitative, because the enrichment result must preserve whether
the upstream quantitative intensity scale was observed through transformation
or declared by the user. Declared evidence produces a role-specific caveat;
observed transformation evidence is recorded without that declared-only caveat.

## Saved Output and Provenance Schemas

Only current PhosPy-generated bundles and provenance payloads are supported.
Legacy saved-result compatibility has been intentionally removed. Regenerate
older development-version outputs rather than relying on schema repair during
loading.

Provenance remains supported for current outputs. Current exact and tolerance
table-hash semantics are unchanged.

### `profile_correlation_shifted_unit_v1`

- What it does:
  transforms profile correlations from `[-1, 1]` to `[0, 1]` using `(r + 1) / 2`.
- Assumptions:
  positive correlation increases support.
- Parameters:
  transform formula, clipping to `[0, 1]`, preserve undefined values as missing.
- Output meaning:
  relative support score; larger means stronger positive agreement.
- Output does not mean:
  calibrated probability or direct evidence of inhibition/activation.
  Negative correlations are treated as lower support, not explicit inhibitory
  evidence.

### `kinase_profile_scoring_v1`

- What it does:
  records kinase profile-construction and scoring behavior, including
  self-inclusion vs leave-one-out semantics.
- Assumptions:
  profile rows can include the same substrate site later scored unless a
  leave-one-out policy is explicitly enabled.
- Parameters:
  profile missing-value strategy, self-inclusion behavior, leave-one-out flag,
  and scoring substrate floors.
- Output meaning:
  explicit provenance of the profile-scoring policy context used for
  downstream support scores.

### `motif_profile_rank_fusion_v1`

- What it does:
  records PhosR-inspired rank-weighted scoring that fuses motif-frequency and
  profile-correlation evidence using rank-derived logarithmic weights.
- Assumptions:
  motif-library size and quantified-substrate count proxy evidence strength.
- Parameters:
  motif/profile weight formulas and fallback/diagnostic flags.
- Output meaning:
  relative downstream support for kinase-site ranking.
- Output does not mean:
  exact PhosR implementation, numerical parity with PhosR, statistical
  enrichment p-value, or calibrated confidence.

### `kinase_library_motif_scoring_v1`

- What it does:
  scores phosphosite sequence windows against Kinase Library-style
  position-specific motif matrices.
- Assumptions:
  Ser/Thr and Tyr residue-class lanes are distinct and must not be substituted
  for each other.
- Parameters:
  score scale, residue classes, upstream/downstream window size, sequence-window
  semantics, and whether reference distributions were supplied.
- Output meaning:
  raw provider-scale motif score sums in the science layer; workflow integration
  records the provider scale and exposes normalized unit-interval motif support
  scores with site and kinase diagnostics.
- Output does not mean:
  calibrated probability, causal kinase-substrate proof, or activity inference
  by itself.

### `simplified_weighted_substrate_activity_v1`

- What it does:
  computes a prediction-weighted kinase activity score and thresholded
  substrate-mean activity-like summary.
- Assumptions:
  predicted substrate support can summarize relative candidate kinase support
  in-run for exploratory interpretation.
- Parameters:
  threshold, `min_substrates`, `top_n_substrates`, and explicit scoring rules.
- Output meaning:
  relative sample-by-kinase activity-like summaries.
- Output does not mean:
  full KSEA-style enrichment statistics, validated kinase activation, or causal
  pathway activity.

### `ksea_zscore_activity_v1`

- What it does:
  computes a KSEA-style z-score inferred kinase activity score.
- Assumptions:
  kinase substrate membership is unweighted after evidence thresholding.
- Parameters:
  evidence threshold, minimum substrates, z-score formula, p-value method, and
  optional q-value adjustment.
- Output meaning:
  substrate-set enrichment z-scores with accompanying p-values.
- Output does not mean PhosR-equivalent kinase activity inference, validated
  kinase activation, or causal pathway activity.

### `ssgsea_substrate_enrichment_activity_v1`

- What it does:
  computes a PhosPy ssGSEA-style rank-walk kinase substrate-set enrichment
  score over site-level phosphosite effect/statistic values.
- Assumptions:
  explicit kinase-substrate membership defines the tested set, and substrate
  concentration near one end of the ranked effect list summarizes candidate
  kinase support.
- Parameters:
  minimum substrates, ranking direction, optional seeded permutation count,
  optional permutation random seed, and optional q-value adjustment.
- Output meaning:
  rank-based substrate-supported kinase scores with optional empirical
  permutation p-values.
- Output does not mean:
  PTM-SEA parity, calibrated causal activity, broader pathway enrichment, or
  validated kinase activation.

### `candidate_substrate_selection_v1`

- What it does:
  records candidate substrate filtering for kinase prediction.
- Assumptions:
  top-k ranking, score-threshold filtering, and inclusion floor jointly define
  usable candidate support.
- Parameters:
  `top_k`, threshold rule, threshold value, inclusion floor, and site
  restriction behavior.
- Output meaning:
  explicit provenance of the candidate-selection rule that gates kinase ranking
  and prediction outputs.

### `signalome_module_candidate_score_v1`

- What it does:
  ranks candidate module counts using within-cluster correlation summaries.
- Assumptions:
  stronger within-cluster profile coherence indicates better candidate module
  structure.
- Parameters:
  requested/resolved candidate-scoring policies, mode, guards, and skip/evaluated
  diagnostics.
- Output meaning:
  candidate module-count support score used for ranking/selection.
- Output does not mean:
  biological certainty or causal regulation evidence.

### `signalome_missing_value_clustering_v1`

- What it does:
  records missing-value handling for clustering distance/tree inputs.
- Assumptions:
  non-finite values are normalized to missing; missing values are imputed for
  clustering internals.
- Parameters:
  missing-value policy name, applicability scope, and whether imputed values
  appear in output tables.
- Output meaning:
  explicit provenance for clustering-matrix preparation rules.

### `signalome_score_preconditioning_v1`

- What it does:
  records row-retention policy for downstream score preconditioning before
  signalome execution.
- Assumptions:
  all-missing rows are unsupported and can be dropped or treated as boundary
  errors depending on policy.
- Parameters:
  preconditioning policy, row-retention rule, and input/dropped/retained row
  counts.
- Output meaning:
  explicit provenance for score-row retention behavior that can change site
  coverage and assignments.

### `protein_module_from_site_membership_v1`

- What it does:
  derives protein module IDs from site-cluster membership incidence patterns.
- Assumptions:
  shared site-cluster membership reflects shared protein-level module context.
- Parameters:
  membership-vector representation and module ID assignment rule.
- Output meaning:
  integer protein module IDs for grouping.
- Output does not mean:
  direct mechanistic proof of shared regulation.

### `preprocessing_stage_order_v1`

- What it does:
  records explicit preprocessing stage order used to construct analysis-ready
  dataset inputs.
- Assumptions:
  stage order is scientifically meaningful and can change transformed values,
  row retention, and derived comparison outputs.
- Parameters:
  configured stage order, default order, and supported stage order metadata.
- Output meaning:
  explicit provenance for preprocessing execution order.

### `peptide_to_site_aggregation_v1`

- What it does:
  records how peptide-level differential statistics are aggregated to site-level
  summaries.
- Assumptions:
  aggregation strategy and variance rules change site-level uncertainty and
  significance behavior.
- Parameters:
  aggregation strategy, minimum peptides per site, missing-variance policy, and
  weighting mode.
- Output meaning:
  explicit provenance for site-level differential summary construction.

## Where Details Live

- [API Guide](api/guide.md) links to each workflow-specific API page.
  config, and result objects.
- [Scientific Coverage](scientific-coverage.md) is the maintained scope and
  coverage matrix.
- [Parity](parity.md) tracks PhosR comparison evidence, fixture locations, and
  parity test references.
- [ADR-0025](adr/adr_0025_competitive_phosphoproteomics_workflow_coverage.md)
  records future coverage direction and guardrails; it is not a current
  support claim by itself.
- [Performance Contracts](performance.md) covers scale limits.
- [ADR Index](adr/index.md) stores maintainer decision records.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

Xiao, D., Yang, P., & Kim, H. J. (2026). *PhosR* (R package manual).
Bioconductor. https://bioconductor.org/packages/release/bioc/manuals/PhosR/man/PhosR.pdf

Xiao, D., Yang, P., & Kim, H. J. (2026). *An introduction to PhosR package*
(Bioconductor vignette). https://bioconductor.org/packages/release/bioc/vignettes/PhosR/inst/doc/PhosR.html
