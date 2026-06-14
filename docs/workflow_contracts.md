# Workflow Contracts

This page gives scientist-facing contracts for each public workflow lane:

1. `AnalysisReadyDatasetBuilder`
2. `DifferentialAnalysisWorkflow`
3. `KinaseWorkflow`
4. `SignalomeWorkflow`
5. `EnrichmentWorkflow`

For support status labels (`parity-gated`, `validated PhosPy implementation`,
`experimental`, `open gap`, `deliberate scope difference`, `not planned`), use
`docs/scientific-coverage.md` as the single maintained matrix.

For roadmap governance, see
[ADR-0025: Competitive Phosphoproteomics Workflow Coverage Roadmap](adr/adr_0025_competitive_phosphoproteomics_workflow_coverage.md).
Roadmap items are not executable workflow contracts until implementation,
documentation, and tests update this page and the scientific coverage matrix.

For executable usage, see:

- `examples/dataset_builder_demo.py`
- `examples/kinase_workflow_demo.py`
- `examples/signalome_workflow_demo.py`
- `docs/api/guide.md#enrichment-contract-boundary` for a minimal offline
  enrichment example

For explicit performance and scale contracts (input sizes, guardrails,
approximation behavior, and failure modes), see `docs/performance.md`.

Environment provenance supports reproducibility audits (runtime versions,
platform, backend, and execution settings), but it does not guarantee bitwise
identical numeric outputs across different machines or dependency builds.

## Result Snapshot Helper Boundary

Result models are typed containers. Public helpers such as `to_dataframe()`,
`*_dataframe()`, `table`, `result_table`, and `to_payload()` return defensive
in-memory snapshots for inspection or handoff only. They are not exporters,
formatters, plotting helpers, report generators, or places to run additional
scientific post-processing.

File writing and bundle publication belong under IO/exporter modules. Plotting
and presentation formatting belong to plotting or reporting adapters outside
the result models.

## AnalysisReadyDatasetBuilder Contract

### Required Input Tables

- `phospho` (required): numeric phosphosite-by-sample matrix.
- `site_metadata` (required): row metadata aligned to `phospho.index`.
- `sample_metadata` (optional): sample metadata aligned to `phospho.columns`.
- `total` (optional): total-protein matrix aligned to `phospho.columns` when total-protein correction or protein-aware preparation is used.

### Identifiers and Alignment Assumptions

- Builder input `site_metadata.index` must match `phospho.index`.
- Builder input may use display labels in standard phosphosite form (for
  example `TSC2;S939;`) only when enough protein context exists to derive
  `site_key`.
- The output `AnalysisReadyPhosphoDataset` uses `site_key` for
  `phospho.index` and `site_metadata.index`; `display_id` remains a
  human-readable label and may repeat.
- Direct `AnalysisReadyPhosphoDataset` construction must already use encoded
  `site_key` indexes and required identity metadata: `site_key`, `display_id`,
  `organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`,
  and `site_sequence`.
- `protein_id` is not part of the base analysis-ready row identity contract.
  Analysis-ready datasets may carry it as optional workflow/source metadata, and
  it may be absent or incomplete until a workflow such as signalome requires it.
- `sample_metadata.index` and `total` sample index must match `phospho.columns` when provided.
- `sample_metadata` column names must be unique when provided.
- If site-matrix construction is enabled, `site_metadata` must provide the required site-identity fields.

### Missing-Value Policy

- Output `AnalysisReadyPhosphoDataset.phospho` is required to be complete (missing-value-free).
- Missing values are either forbidden or handled explicitly by configured preprocessing policy.
- Missing-data handling runs before normalisation in preprocessing stage order.
- For non-MinProb missing-data policies, missing-data handling runs before optional `intensity_transform.policy="log2"`.
- When `missing_data.policy="impute_row_median"`, row-median imputation is deterministic.
- Row-median imputation is not left-censored imputation.
- Row-median-imputed values are explicit replacements and are not evidence that the original values were observed.
- `missing_data.policy="impute_minprob"` is explicit opt-in left-censored random imputation.
- `impute_minprob` requires `intensity_transform.policy="log2"` and runs after the log2 stage.
- `impute_minprob` is deterministic for fixed seed and records per-column distribution parameters and dropped-row diagnostics.
- `missing_data.policy="impute_knn"` is explicit opt-in nearest-neighbour imputation.
- `impute_knn` requires explicit `k`, `distance="nan_euclidean"`, and `max_missing_fraction_per_row`.
- `impute_knn` drops rows above the configured row-missingness threshold before imputation and reports those rows as not imputable.
- `impute_knn` is deterministic for fixed input/config and records neighbour count, distance metric, and row-drop diagnostics.
- Output `comparisons` (when produced) must also be numeric and missing-value-free.

### Transformation and Normalisation Policy

- Controlled only by `DatasetPreprocessingConfig`:
  - intensity transform (`identity` or `log2`)
  - normalisation (`none`, `median_center`, `quantile`)
  - optional total-protein correction
  - optional protein-aware preparation
  - optional batch correction with `linear_residualize_batch`
  - optional site-matrix and comparison construction
  - optional `ruv_readiness` reporting for possible future RUV/SPS/RUV-III
    preprocessing
- No hidden transforms are applied outside this configuration.
- `total_protein_correction.policy="subtract_log_total"` is a direct
  phosphosite-minus-total transformation on log-scale values. It changes
  phosphosite matrix values and quantitative meaning. It is not protein-aware
  modelling and is not MSstatsPTM-style inference.
- `protein_aware_preparation.policy="prepare_model_inputs"` is
  preparation-only model-input preparation. It prepares aligned
  phosphosite/protein input contracts and diagnostics only. It does not modify
  phosphosite values, subtract total protein, normalise intensities, run joint
  PTM/protein modelling, or adjust downstream differential models.
- Public builder protein-aware preparation maps explicit protein identifiers
  from `site_metadata` (`protein_accession`, `protein_id`, or
  `protein_group_id`) to `total.index`. Gene-symbol matching is not the public
  default.
- `protein_mapping_policy="require_unambiguous"` excludes missing and ambiguous
  mappings from preparation. `protein_mapping_policy="allow_missing_with_report"`
  reports missing site protein identifiers or missing total-protein rows as
  phospho-only fallback rows. Ambiguous site-to-protein or total-protein-row
  mappings remain excluded.
- Protein-aware preparation requires total-protein input data. Missing
  total-protein rows are per-site diagnostics; missing `total` input is an
  error when `policy="prepare_model_inputs"` is selected.
- Phospho and total sample columns must match in the same order for builder
  protein-aware preparation. Missing, extra, or reordered total-protein sample
  columns are diagnostics and make sites ineligible; the builder does not
  reorder matrices for this preparation stage.
- Phospho and total transformation states must be compatible: same scale kind
  and transformed flag. Incompatible states are diagnostics and make sites
  ineligible for protein-aware preparation.
- `linear_residualize_batch` is fixed-effect residualisation of batch terms
  while preserving condition effects by design. It requires explicit batch and
  condition columns in `sample_metadata`; those columns are used only for this
  preprocessing step and do not define differential-analysis design.
- Confounded batch/condition designs are rejected before correction because the
  method cannot preserve condition effects when batch and condition are
  perfectly confounded.
- `linear_residualize_batch` is not ComBat, not RUV, not limma
  `removeBatchEffect` parity, and not mixed-effects modelling.
- Prefer intent presets (`DatasetPreprocessingConfig.default()`,
  `DatasetPreprocessingConfig.strict()`, and
  `DatasetPreprocessingConfig.from_raw_phosphosite_table()`) for common lanes.

### Scoring Policy

- Not applicable. This workflow prepares analysis-ready quantitative inputs.

### Candidate/Module Selection Policy

- Not a kinase/module-selection workflow.
- If duplicate sites are resolved, selection behavior is controlled by `site_matrix.duplicate_site_policy`.

### Provenance Guarantees

- Dataset provenance records preprocessing-stage execution and table fingerprints.
- `intensity_scale_state` establishment provenance records establishment mode
  (`declared`, `transformed`, `identity`, or `derived`) and establishment source.
- Normalisation provenance is explicit: method, parameters, matrix-shape before/after,
  per-sample summary before/after, and row/column drop diagnostics.
- Each table fingerprint carries both `exact_hash_*` (audit/regression) and `tolerance_hash_*` (cross-platform tolerant comparison) metadata; legacy `hash_*` fields are compatibility aliases of the tolerance hash.
- `provenance.workflow_parameters["preprocessing_plan"]` includes both `stage_order` and
  `resolved_stage_order` with per-stage order index and rationale.
- `intensity_scale_state` and `processing_state` are attached and validated at boundary.
- `preprocessing_report` provides row-level and operation-level preprocessing audit tables.
- `preprocessing_report.batch_correction` is a typed `BatchCorrectionReport`.
  It records the method, status (`"disabled"`, `"applied"`, or `"rejected"`),
  batch and condition columns, observed levels, matrix shapes, condition
  preservation policy, confounding-check status, warnings, and limitations.
- `preprocessing_report.protein_aware_preparation` is a typed
  `ProteinAwarePreparationReport` when preparation is enabled. It records the
  preparation policy, protein-mapping policy, eligibility counts,
  site-eligibility table, missing-total and ambiguous-mapping diagnostics,
  sample-alignment diagnostics, transformation-state diagnostics, and explicit
  limitations: no phospho-matrix modification, no total-protein subtraction, no
  normalisation, no joint modelling, no differential model adjustment, and no
  MSstatsPTM equivalence claim.
- Dataset run provenance includes the active protein-aware preparation summary
  under `provenance.workflow_parameters["protein_aware_preparation"]` when the
  preparation stage runs.
- `processing_state.ruv_readiness` reports whether required controls/groups/batch and
  missingness provenance are present for possible future RUV/SPS/RUV-III
  preprocessing. This is report-only and does not select SPS controls, run
  correction, or block dataset construction.

### Known Limitations

- This workflow does not infer kinase activity or signalome structure.
- Protein identity is not derived automatically from display labels; provide
  protein context for `site_key` derivation and explicit, complete `protein_id`
  only for downstream signalome analysis.
- Protein-aware preparation is preparation-only and is not a supported joint
  PTM/protein differential model. Current differential analysis does not
  consume `ProteinAwarePreparationResult`, does not adjust differential models
  from its covariate matrix, and no MSstatsPTM-style inference or equivalence is
  claimed.
- The only executable batch-correction preprocessing method is
  `linear_residualize_batch`. Broader batch-effect modelling, ComBat, RUV,
  RUV/SPS/RUV-III correction, limma `removeBatchEffect` parity, and
  mixed-effects modelling are not provided by the dataset builder.

### Expected Output Tables

- `dataset.phospho` indexed by `site_key`
- `dataset.site_metadata` indexed by `site_key`, with `site_key`, `display_id`,
  `organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`,
  `site_sequence`, and optional workflow/source metadata such as `protein_id`.
  `protein_id` is signalome grouping metadata, not a `site_key` identity field.
- optional `dataset.sample_metadata`, `dataset.total`, `dataset.comparisons`
- optional `dataset.protein_aware_preparation` containing
  `ProteinAwarePreparationResult` when protein-aware preparation is enabled
- optional preprocessing report tables (`row_counts`, `operations`, `row_audit`, and sidecars)

## EnrichmentWorkflow Contract

### Required Input Tables

- `selected_identifiers` or `input_table` with the configured
  `identifier_column`; exactly one source is required.
- `set_collection`: a caller-supplied `GeneSetCollection`,
  `PtmSetCollection`, or homogeneous `EnrichmentSetCollection`.
- `background_universe`: explicit, non-empty, and in the same identifier
  namespace as the selected identifiers and set collection.
- `config`: `EnrichmentConfig` with `method="over_representation"` and an
  explicit multiple-testing correction policy.

### Identifiers and Alignment Assumptions

- `identifier_kind` is required and must match the collection:
  `gene_symbol` and `protein_id` are gene-level; `site_key`, `display_id`, and
  `phosphosite` are PTM/site-level.
- Gene-level and site-level enrichment are not interchangeable. The workflow
  does not convert gene identifiers to phosphosite identifiers or collapse
  site identifiers to genes.
- Selected identifiers must already be members of `background_universe`.
  Validation fails instead of silently dropping selected identifiers.
- Set members outside the background are excluded from each ORA test and
  reported in summaries/diagnostics.

### Statistical Policy

- Native enrichment is offline over-representation analysis only.
- ORA uses a hypergeometric test over the explicit background universe and the
  supplied set collection.
- Supported multiple-testing correction values are `"benjamini_hochberg"` and
  `"none"`.
- Background universes are not inferred from datasets, reference bundles, or
  set membership.

### Provenance Guarantees

- `EnrichmentWorkflowResult.provenance` is populated by workflow execution.
- Provenance records method, identifier column and kind, collection kind,
  analysis level, explicit background universe size, selected identifier
  count, selected identifier source, set collection source/name/version
  metadata when supplied, multiple-testing correction, and result-table
  fingerprints.
- Provenance records the offline/no-online-resource policy and explicit
  limitations.

### Known Limitations

- GO, KEGG, Reactome, PTM-SEA, and PTMsigDB resources are not bundled by this
  feature unless the caller supplies them as ordinary local collections.
- Online Enrichr, gseapy, clusterProfiler, or other remote-service calls are
  not part of the core workflow.
- This workflow does not implement GSEA, ssGSEA, or PTM-SEA. Kinase activity
  ssGSEA-style scoring is a separate kinase-activity lane, not pathway or
  PTM-SEA enrichment support.

### Expected Output Tables

- `result.table` / `result.to_dataframe()` with one row per tested term,
  including overlap counts, overlap identifiers, p-value, adjusted p-value,
  correction method, and enrichment ratio.
- `result.background_summary`, `result.set_collection_summary`,
  `result.method_metadata`, and `result.diagnostics` sidecars for execution
  audit.

## DifferentialAnalysisWorkflow Contract

### Required Input Tables

- `dataset.phospho`: numeric feature-by-sample table from `AnalysisReadyPhosphoDataset`.
- `design`: typed `ExperimentalDesign` sample records aligned to dataset sample
  IDs (`sample_id`, `condition`, optional replicate/batch/`block_id` metadata).
- `contrasts`: typed `Contrast` condition-vs-condition definitions.

### Experimental-Design Contract Policy

- Differential analysis does not infer conditions from sample names.
- Sample/design alignment is validated before statistical execution.
- `technical_replicate_policy` is explicit and defaults to `reject`.
- Supported explicit technical-replicate aggregation policies are `mean` and
  `median`.
- By default, dataset and design must reference the same sample set.
- `config.allow_design_subset=True` is the only supported path to analyze an explicit
  subset of samples.
- Duplicate design sample IDs, missing condition labels, and unknown contrast
  conditions are hard validation errors.
- Batch can be declared as a fixed-effect covariate and modelled as an ordinary
  design term when the resolved design is full rank and requested contrasts are
  estimable. Categorical and continuous fixed-effect covariates are supported
  under the same validation policy.
- Fixed-effect batch terms are not batch correction and do not implement
  ComBat, RUV, `removeBatchEffect`, `duplicateCorrelation`, or mixed-effects
  modelling.
- Repeated `biological_replicate_id` values within condition groups are treated
  as technical replicates and require an explicit aggregation policy.
- Technical-replicate aggregation requires `biological_replicate_id` for every
  design sample and consistent optional group fields (`batch`, `block_id`) within
  each condition+biological-replicate group.
- Paired/block design intent is explicit via `config.paired_design_policy`,
  which defaults to `"reject"`. `paired_design_policy="fixed_block"` is an
  opt-in fixed-effect block policy. It does not infer `block_id` and does not
  enable mixed-effects modelling. Each sample must have `block_id`; each block
  must contain at least two samples and both sides of every requested contrast;
  the resolved block design must be full rank and contrasts must be estimable.
- Current parity-protected lane is two-condition unpaired simple contrasts.
- Correlated-replicate and mixed-effects modelling are not executable in this
  release.
- Missing values are rejected at `AnalysisReadyPhosphoDataset` boundary before
  differential execution.

### Empirical-Bayes Moderation Policy

- `config.empirical_bayes.method="standard"` applies limma-style moderated variance.
- `config.empirical_bayes.method="robust"` applies winsorized robust hyperparameter estimation and outlier-aware prior-df shrinkage.
- `config.empirical_bayes.trend=True` enables mean-intensity-dependent prior variance (limma-trend style).
- `config.empirical_bayes.winsor_tail_p` controls robust winsor tail clipping; this is used only when `method="robust"`.
- Moderation changes residual variance estimates, moderated t-statistics, and p-values.
- Moderation does not alter `logFC`; fold-change estimates remain OLS contrast estimates.

### Expected Output Tables

- Each contrast table is indexed by encoded protein-scoped `site_key`.
- Each contrast table includes `site_key`, `display_id`, `organism`,
  `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`, `logFC`,
  `t`, `P.Value`, and `adj.P.Val`.
- `site_key` values in the column exactly match the index.
- Required protein context from `dataset.site_metadata` is preserved, and
  optional protein metadata such as `protein_id` is preserved when present.
- Duplicate `display_id` values remain distinct rows when `site_key` differs.
- Display-indexed, `GENE;SITE;`-keyed, arbitrary-keyed, and stat-only public
  `DifferentialAnalysisResult` tables are invalid.

### Diagnostics and Provenance

- Result output includes method/trend flags and prior-parameter diagnostics.
- Trend mode includes mean-intensity vs residual-variance trend diagnostics as typed data objects.
- Core differential output does not require plotting; plotting adapters are out of scope for this core contract.

## KinaseWorkflow Contract

### Required Input Tables

- `dataset.phospho` and `dataset.site_metadata` from an `AnalysisReadyPhosphoDataset`.
- Resolved references with:
  - `kinase_substrate_map` (`kinase`, `substrate_site`)
  - `site_sequences` (`site_sequence` indexed by display site ID)
- For `scoring_config.scoring_mode="kinase_library_motif"` or
  `"combined_profile_motif"`, the normal resolved reference requirements still
  apply. The projected `kinase_substrate_map` must overlap dataset `display_id`
  values and must contain at least one kinase with
  `scoring_config.min_substrates` quantified substrates.
- Kinase Library-style workflow modes additionally require
  `KinaseWorkflowRequest.kinase_library_resource` to be a local
  `KinaseLibraryResource` with compatible organism metadata, a central
  phospho-residue sequence-window definition, matrix lanes, score scale, and
  provenance.
- Site sequences are required for workflow scoring rows. Reference
  `site_sequences` are projected by `display_id`; dataset
  `site_metadata.site_sequence` can supplement missing reference rows. If no
  scoring row has resolved sequence support, interpretation fails before
  scoring.
- The supplied Kinase Library-style resource must include at least one
  residue-class lane matching the resolved scoring-site sequences. For example,
  an all-Ser/Thr scoring set with only Tyr resource lanes fails validation at
  `kinase.interpreter.kinase_library_resource_usability`.

### Identifiers and Alignment Assumptions

- Quantified rows are keyed by `site_key`; reference `substrate_site` and
  `site_sequences` display labels are projected onto those rows through
  dataset `display_id` metadata before scoring.
- This projection is an explicit `display_id` -> `site_key` mapping layer.
  Reference display IDs remain reference/display identifiers and are not
  converted into analysis-ready row identity.
- Sample columns are used as provided by the dataset; scoring compares site profiles across these aligned samples.
- Reference resolution can come from `ReferencePreset` or explicit `ReferenceBundle`.
- Dataset `site_metadata.site_sequence` values can supplement missing reference sequences.
- When a dataset sequence conflicts with a reference sequence for the same site, behavior is controlled by `KinaseWorkflowRequest.site_sequence_conflict_policy`:
  - `"prefer_reference"` (default): keep the reference sequence
  - `"prefer_dataset"`: override with the dataset sequence
  - `"error"`: fail fast with conflict diagnostics and a clear next action

### Missing-Value Policy

- Dataset input is expected to be analysis-ready (complete matrix at boundary).
- Correlation-based scores can still become missing for unsupported cases (for example zero-variance profiles).
- `pred_mat` keeps missing entries where a site-kinase pair is not selected/scored for output.

### Transformation and Normalisation Policy

- No additional intensity transform or normalisation is done in this workflow.
- The workflow consumes dataset values as provided by dataset-building preprocessing.
- Mixed corrected/uncorrected total-protein quantitative meaning is rejected by
  default. Set
  `scoring_config.allow_mixed_total_protein_quantitative_meaning=True` only
  when mixed-state input is intentional and scientifically justified.

### Scoring Policy

- Kinase profiles are built from quantified substrates with `min_substrates >= 2`.
- Site-vs-profile Pearson correlations are shifted to `[0, 1]` support scores.
- Motif and profile evidence are fused by rank-weighted fusion with profile fallback.
- `profile_missing_value_strategy` controls strict vs skip-missing profile medians.
- The default `scoring_mode="phosr_rank_weighted"` uses the PhosR-style
  profile plus motif-frequency rank-fusion lane.
- `scoring_mode="kinase_library_motif"` still builds profile context from the
  resolved kinase-substrate map to validate eligible workflow kinases, but the
  authoritative downstream score matrix is `kinase_library_motif_scores`.
- Kinase Library-style workflow modes do not silently fall back to PhosR-style
  motif scoring. Missing `kinase_library_resource`, organism/window mismatch,
  or missing matching residue-class lanes fail validation.
- Raw science-layer `score_kinase_library_motifs` outputs preserve the
  caller-supplied provider score scale. Workflow
  `kinase_library_motif_scores` are per-kinase min-max normalized support
  scores for within-run ranking and prediction support.
- Substrate-map activity inference is separate from motif-based kinase
  prediction support. Kinase Library motif scores alone do not imply activity;
  activity is produced only by an enabled activity method over workflow
  prediction outputs.
- Prefer `KinaseScoringConfig.default()` or
  `KinaseScoringConfig.strict_missing_values()` before low-level field tuning.

### Candidate/Module Selection Policy

- Candidate substrates per kinase are selected from top-`k` downstream scores with score `> 0`.
- Deterministic mode ranks kinases by mean candidate support and keeps top configured kinases.
- Adaptive mode uses ensemble sampling policy from `prediction_config`.
- Adaptive mode requires explicit `prediction_config.random_state`; missing seed
  is a hard validation error by design.
- Prefer `KinasePredictionConfig.deterministic()` and
  `KinasePredictionConfig.adaptive_reproducible(random_state=...)` presets.
- `adaptive_policy="r_parity"` is a parity-oriented mode, not the default
  production recommendation.

### Provenance Guarantees

- Provenance includes input/output table fingerprints, environment, resolved reference provenance, and active scientific policy records.
- Fingerprints expose explicit exact-vs-tolerance hash metadata; use exact hashes for scientific audit and tolerance hashes only for tolerant comparisons.
- Workflow parameters for scoring, prediction, and optional activity are serialized in provenance.
- `workflow_parameters.scoring_diagnostics` includes motif sequence-validation
  counts and motif site-sequence coverage diagnostics (for example
  `total_sites_considered`, `sites_with_valid_site_sequence`,
  `sites_without_valid_site_sequence`, and `site_sequence_coverage_fraction`).

### Known Limitations

- Scores are relative support values within a run, not calibrated probabilities.
- Kinase Library-style workflow scoring is a PhosPy workflow integration over
  caller-supplied local resources. It is not an official Kinase Library
  implementation and is not parity-tested against an official predictor.
- Activity supports two explicit methods:
  - `simplified_weighted_substrate_activity_v1` (heuristic weighted substrate activity)
  - `ksea_zscore_v1` (KSEA-style z-score substrate-set enrichment)
- KSEA-style activity is not equivalent to PhosR kinase activity inference.
- In this release, `ReferencePreset.AUTO` bundled runtime data is rat-first.

### Expected Output Tables

- Scoring, default mode: `profile_scores`, `rank_weighted_fusion_scores` (plus optional diagnostics)
- Scoring, Kinase Library motif mode: `profile_scores`,
  `kinase_library_motif_scores`, `kinase_library_site_diagnostics`, and
  `kinase_library_kinase_diagnostics`
- Scoring, combined profile/motif mode: `profile_scores`,
  `kinase_library_motif_scores`, `combined_profile_motif_scores`, diagnostics,
  and optional fusion weights
- Prediction: `pred_mat`, `substrate_list`
- Optional activity: primary `activity_scores` matrix (with `weighted_activity` as compatibility alias), `thresholded_substrate_mean_activity`, `thresholded_substrate_counts`, `target_counts`, `target_table`, optional `statistics_table`
- Primary matrices are indexed by `site_key`; site-level tables that materialize
  row identity include both `site_key` and `display_id`.

## SignalomeWorkflow Contract

### Required Input Tables

- `kinase_result.prediction_result.pred_mat` (non-empty, numeric).
- Upstream downstream score matrix from kinase scoring (`rank_weighted_fusion_scores` or `profile_scores` fallback).
- Complete `kinase_result.dataset.site_metadata.protein_id` values for retained
  interpreted sites. This is signalome-specific protein grouping metadata, not
  base dataset row identity.
- A valid upstream `AnalysisReadyPhosphoDataset` with `site_key` row identity;
  signalome does not reinterpret display IDs as row identity.

### Identifiers and Alignment Assumptions

- Signalome runs on the shared intersection of `site_key` values across
  dataset, prediction matrix, and downstream score matrix.
- `display_id` remains display metadata and may repeat.
- Kinase columns must overlap between prediction and downstream score matrices.
- Site and kinase identifiers must be unique within aligned matrices.

### Missing-Value Policy

- All-missing downstream-score rows are handled by `config.validation.score_preconditioning_policy`:
  - `error_on_drop` (default): fail fast
  - `allow_and_report`: explicit opt-in to drop and report
- Non-finite values are rejected.
- Partially missing values are retained for correlation-based steps with explicit diagnostics.
- Prefer `SignalomeConfig.strict()` and
  `SignalomeConfig.permissive_missing_scores()` presets for this choice.

### Transformation and Normalisation Policy

- No additional biological normalisation is applied.
- For clustering internals, missing entries are handled by workflow-defined missing-value policy and imputation diagnostics are recorded.
- Mixed corrected/uncorrected total-protein quantitative meaning is rejected by
  default. Set
  `config.validation.allow_mixed_total_protein_quantitative_meaning=True` only
  when mixed-state input is intentional and scientifically justified.

### Scoring Policy

- Module-count candidate scoring uses within-cluster correlation summaries over downstream score profiles.
- Kinase-network edge scores are pairwise correlations on finite paired observations.
- Edge inclusion is controlled by `network_policy` and `network_correlation_threshold`.

### Candidate/Module Selection Policy

- Kinase-supported substrates are selected by `substrate_support_cutoff`.
- Module count is either explicit (`module_count`) or automatically selected from configured candidate thresholds.
- Module composition uses `assignment_policy` (`cutoff_binary` or `weighted_top`).
- For workflows where automatic module-count candidate evaluation is expensive,
  use `SignalomeConfig.sampled_candidate_scoring()`.
- This selects sampled candidate scoring, preserves exact-tree and
  full-scoring scale-guard checks, and is not a general bypass for large input
  limits.

### Provenance Guarantees

- Provenance includes upstream kinase provenance payload, input/output table fingerprints, resolved configuration, scale-guard diagnostics, and scientific policy records.
- Table fingerprints include explicit exact and tolerance hash metadata; exact hashes are the authoritative audit-trail fingerprints.
- Alignment and score-preconditioning diagnostics are included for interpretability and audit.

### Known Limitations

- Module/network outputs are derived summaries, not direct evidence of causal regulation.
- Signalome requires explicit protein identifiers; gene-symbol prefixes in
  display labels are not a protein-identity substitute.
- `candidate_scoring_policy="sampled"` approximates candidate module-count scoring only; tree generation remains exact-policy governed and is reported in provenance.

### Expected Output Tables

- `module_assignments.table`
- `signalome_modules.table`
- `kinase_network.edges` (plus `nodes` and `candidate_correlations`)
- `expanded_signalome`
- `site_membership`
- `protein_site_context`
