# Workflow Contracts

This page gives scientist-facing contracts for each public workflow lane:

1. `AnalysisReadyDatasetBuilder`
2. `KinaseWorkflow`
3. `SignalomeWorkflow`

For executable usage, see:

- `examples/dataset_builder_demo.py`
- `examples/kinase_workflow_demo.py`
- `examples/signalome_workflow_demo.py`

For explicit performance and scale contracts (input sizes, guardrails,
approximation behavior, and failure modes), see `docs/performance.md`.

## AnalysisReadyDatasetBuilder Contract

### Required Input Tables

- `phospho` (required): numeric phosphosite-by-sample matrix.
- `site_metadata` (required): row metadata aligned to `phospho.index`.
- `sample_metadata` (optional): sample metadata aligned to `phospho.columns`.
- `total` (optional): total-protein matrix aligned to `phospho.columns` when total-protein correction is used.

### Identifiers and Alignment Assumptions

- `site_metadata.index` must exactly match `phospho.index`.
- Site IDs are expected in canonical phosphosite form (for example `TSC2;S939;`).
- `sample_metadata.index` and `total` sample index must match `phospho.columns` when provided.
- If site-matrix construction is enabled, `site_metadata` must provide the required site-identity fields.

### Missing-Value Policy

- Output `AnalysisReadyPhosphoDataset.phospho` is required to be complete (missing-value-free).
- Missing values are either forbidden or handled explicitly by configured preprocessing policy.
- Output `comparisons` (when produced) must also be numeric and missing-value-free.

### Transformation and Normalisation Policy

- Controlled only by `DatasetPreprocessingConfig`:
  - intensity transform (`identity` or `log2`)
  - normalisation (`none`, `median_center`, `quantile`)
  - optional total-protein correction
  - optional site-matrix and comparison construction
- No hidden transforms are applied outside this configuration.
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
- `intensity_scale_state` and `processing_state` are attached and validated at boundary.
- `preprocessing_report` provides row-level and operation-level preprocessing audit tables.

### Known Limitations

- This workflow does not infer kinase activity or signalome structure.
- `protein_id` is not derived automatically from site IDs; provide it explicitly for downstream signalome analysis.

### Expected Output Tables

- `dataset.phospho`
- `dataset.site_metadata`
- optional `dataset.sample_metadata`, `dataset.total`, `dataset.comparisons`
- optional preprocessing report tables (`row_counts`, `operations`, `row_audit`, and sidecars)

## KinaseWorkflow Contract

### Required Input Tables

- `dataset.phospho` and `dataset.site_metadata` from an `AnalysisReadyPhosphoDataset`.
- Resolved references with:
  - `kinase_substrate_map` (`kinase`, `substrate_site`)
  - `site_sequences` (`site_sequence` indexed by site ID)

### Identifiers and Alignment Assumptions

- Quantified phosphosite IDs must align with reference `substrate_site` and `site_sequences` indices.
- Sample columns are used as provided by the dataset; scoring compares site profiles across these aligned samples.
- Reference resolution can come from `ReferencePreset` or explicit `ReferenceBundle`.

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
- Prefer `KinaseScoringConfig.default()` or
  `KinaseScoringConfig.strict_missing_values()` before low-level field tuning.

### Candidate/Module Selection Policy

- Candidate substrates per kinase are selected from top-`k` downstream scores with score `> 0`.
- Deterministic mode ranks kinases by mean candidate support and keeps top configured kinases.
- Adaptive mode uses ensemble sampling policy from `prediction_config`.
- Prefer `KinasePredictionConfig.deterministic()` and
  `KinasePredictionConfig.adaptive_reproducible(random_state=...)` presets.

### Provenance Guarantees

- Provenance includes input/output table fingerprints, environment, resolved reference provenance, and active scientific policy records.
- Workflow parameters for scoring, prediction, and optional activity are serialized in provenance.
- `workflow_parameters.scoring_diagnostics` includes motif sequence-validation
  counts and motif site-sequence coverage diagnostics (for example
  `total_sites_considered`, `sites_with_valid_site_sequence`,
  `sites_without_valid_site_sequence`, and `site_sequence_coverage_fraction`).

### Known Limitations

- Scores are relative support values within a run, not calibrated probabilities.
- Activity outputs are a simplified weighted substrate activity score
  (`simplified_weighted_substrate_activity_v1`), not full KSEA enrichment and
  not PhosR-equivalent kinase activity inference.
- In this release, `ReferencePreset.AUTO` bundled runtime data is rat-first.

### Expected Output Tables

- Scoring: `profile_scores`, `rank_weighted_fusion_scores` (plus optional diagnostics)
- Prediction: `pred_mat`, `substrate_list`
- Optional activity: `weighted_activity`, `thresholded_substrate_mean_activity`, `thresholded_substrate_counts`, `target_counts`, `target_table`

## SignalomeWorkflow Contract

### Required Input Tables

- `kinase_result.prediction_result.pred_mat` (non-empty, numeric).
- Upstream downstream score matrix from kinase scoring (`rank_weighted_fusion_scores` or `profile_scores` fallback).
- `kinase_result.dataset.site_metadata.protein_id` for retained interpreted sites.

### Identifiers and Alignment Assumptions

- Signalome runs on the shared intersection of site IDs across dataset, prediction matrix, and downstream score matrix.
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
- For larger runs, prefer `SignalomeConfig.large_dataset()` to select sampled
  candidate scoring while preserving scale-guard checks.

### Provenance Guarantees

- Provenance includes upstream kinase provenance payload, input/output table fingerprints, resolved configuration, scale-guard diagnostics, and scientific policy records.
- Alignment and score-preconditioning diagnostics are included for interpretability and audit.

### Known Limitations

- Module/network outputs are derived summaries, not direct evidence of causal regulation.
- Signalome requires explicit protein identifiers; gene-symbol prefixes in site IDs are not a protein-identity substitute.
- `candidate_scoring_policy="sampled"` approximates candidate module-count scoring only; tree generation remains exact-policy governed and is reported in provenance.

### Expected Output Tables

- `module_assignments.table`
- `signalome_modules.table`
- `kinase_network.edges` (plus `nodes` and `candidate_correlations`)
- `expanded_signalome`
- `site_membership`
- `protein_site_context`
