# Kinase Workflow

`KinaseWorkflow` scores kinase-substrate evidence, predicts candidate kinase
regulation, and can optionally compute kinase activity tables. Use it after
building an `AnalysisReadyPhosphoDataset`.

## When to Use This Workflow

Use this workflow when your analysis-ready phosphosite dataset has site
sequences, protein-scoped row identity, and suitable kinase reference data.

Good fits:

- rat beginner runs with `ReferencePreset.AUTO`
- custom human or mouse runs with an explicit `ReferenceBundle`
- default PhosR-style rank-weighted scoring
- explicit Kinase Library-style motif scoring with a caller-supplied
  `KinaseLibraryResource`
- optional activity summaries from workflow prediction output

Scores are relative support values within a run, not calibrated probabilities
or proof of causal regulation.

## Inputs

`KinaseWorkflowRequest.dataset` must be an `AnalysisReadyPhosphoDataset` with:

- numeric, missing-value-free phosphosite values
- rows keyed by `site_key`
- `display_id` metadata for reference projection
- `site_sequence` support for scoring rows
- explicit localisation policy applied upstream when site-level localisation
  matters

References are supplied as `ReferencePreset` or `ReferenceBundle`. Reference
substrates may use display IDs; the workflow projects them to dataset `site_key`
rows through `dataset.site_metadata.display_id`.

## Localisation Prerequisite

Site-level kinase interpretation should start from a dataset that failed fast on
missing or low-confidence localisation. A typical upstream policy is:

```python
from phospy.api import DatasetLocalisationConfig, DatasetPreprocessingConfig

preprocessing = DatasetPreprocessingConfig(
    localisation=DatasetLocalisationConfig(
        mode="require_threshold",
        confidence_column="localisation_confidence",
        min_confidence=0.75,
    )
)
```

With this policy, dataset build fails when the localisation column is missing,
invalid, missing per row, or below threshold.

## Request Object

Use `KinaseWorkflowRequest`.

Important fields:

| Field | Meaning |
| --- | --- |
| `dataset` | The `AnalysisReadyPhosphoDataset` to score. |
| `references` | `ReferencePreset` or explicit `ReferenceBundle`. |
| `scoring_config` | `KinaseScoringConfig` for scoring mode, substrate floors, diagnostics, localisation, and mixed total-protein guardrails. |
| `prediction_config` | `KinasePredictionConfig` for deterministic or adaptive prediction. |
| `activity_config` | `KinaseActivityConfig`, `None`, or disabled config for optional activity output. |
| `site_sequence_conflict_policy` | Handles dataset/reference sequence conflicts: `"prefer_reference"`, `"prefer_dataset"`, or `"error"`. |
| `reference_display_ambiguity_policy` | Handles one-display-to-many-`site_key` reference projection: `"error"` or `"allow_with_diagnostics"`. |
| `kinase_library_resource` | Required only for Kinase Library-style workflow scoring modes. |

Constructing the request records intent only. `KinaseWorkflow.run(...)` validates
dataset, reference, config, localisation, sequence, and projection compatibility
before interpretation and execution.

## Request Configuration

Use these config objects:

- `KinaseScoringConfig`
- `KinasePredictionConfig`
- `KinaseActivityConfig`
- `LocalisationRequirement`

Important `KinaseScoringConfig` fields:

| Field | Default | Notes |
| --- | --- | --- |
| `scoring_mode` | `"phosr_rank_weighted"` | Supported modes: `"phosr_rank_weighted"`, `"kinase_library_motif"`, `"combined_profile_motif"`. |
| `min_substrates` | `2` | Minimum quantified substrates for kinase support. |
| `include_diagnostic_scoring_tables` | `False` | Adds non-primary diagnostic scoring tables. |
| `profile_missing_value_strategy` | `"strict"` | Use `"median_skipna"` only when skipping missing profile values is intended. |
| `localisation_requirement` | `LocalisationRequirement()` | Workflow-level localisation requirement. |
| `allow_mixed_total_protein_quantitative_meaning` | `False` | Keep `False` unless mixed corrected/uncorrected rows are intended. |

Important `KinasePredictionConfig` fields:

| Field | Default | Notes |
| --- | --- | --- |
| `mode` | `"deterministic_ranking"` | Use `"adaptive_ensemble"` only with an explicit `random_state`. |
| `top_k` | `30` | Top substrate sites retained per kinase. |
| `deterministic_max_selected_kinases` | `10` | Breadth of deterministic kinase selection. |
| `adaptive_ensemble_runs` | `10` | Adaptive ensemble executions. |
| `adaptive_policy` | `"stable"` | `"r_parity"` is parity-oriented, not the recommended production default. |
| `n_iterations` | `5` | Adaptive sampling iterations. |
| `random_state` | `None` | Required for adaptive mode. |

Important `KinaseActivityConfig` fields:

| Field | Default | Notes |
| --- | --- | --- |
| `enabled` | `True` | Set `activity_config=None` or `enabled=False` to skip activity. |
| `method` | `"simplified_weighted_substrate_activity"` | Supported methods: `"simplified_weighted_substrate_activity"`, `"ksea_zscore"`, `"ssgsea_substrate_enrichment"`. |
| `threshold` | `0.6` | Prediction support threshold. |
| `min_substrates` | `3` | Weighted activity substrate floor. |
| `top_n_substrates` | `20` | Weighted activity substrate cap. |
| `ksea_min_substrates` | `5` | KSEA-style substrate floor. |
| `ssgsea_min_substrates` | `5` | ssGSEA-style substrate floor. |
| `ssgsea_random_seed` | `0` | Required when permutations are enabled. |

## Running the Workflow

```python
from phospy import KinaseWorkflow
from phospy.api import KinaseWorkflowRequest, ReferencePreset

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        activity_config=None,
    )
)
```

`ReferencePreset.AUTO` is suitable for the bundled rat beginner lane. For human
or mouse analysis, pass an explicit `ReferenceBundle` with local provenance.

## Result Object

`KinaseWorkflow.run(...)` returns `KinaseWorkflowResult`.

Important fields:

| Field | Meaning |
| --- | --- |
| `dataset` | Input `AnalysisReadyPhosphoDataset`. |
| `references` | Resolved `ReferenceBundle`. |
| `scoring_result` | `KinaseScoringResult` stage output. |
| `prediction_result` | `KinasePredictionResult` stage output. |
| `activity_result` | Optional `KinaseActivityResult`. |
| `eligibility_report` | Optional compact eligibility counters. |
| `site_attrition_summary` | Optional preprocessing/scoring attrition counters. |
| `provenance` | Workflow provenance. |
| `input_dataset_preprocessing_report` | Input dataset preprocessing report when available. |

Important nested result fields:

- `KinaseScoringResult.profile_scores`
- `KinaseScoringResult.rank_weighted_fusion_scores`
- `KinaseScoringResult.authoritative_scores`
- `KinaseScoringResult.score_source`
- `KinasePredictionResult.pred_mat`
- `KinasePredictionResult.substrate_list`
- `KinaseActivityResult.activity_matrix`
- `KinaseActivityResult.substrate_count_matrix`
- optional activity `p_value_matrix`, `q_value_matrix`, and
  `statistics_table`

`activity_result.activity_matrix` is the method-neutral primary activity matrix.
Deprecated compatibility aliases such as `activity_scores` and
`weighted_activity` are not preferred for new documentation or code.

## Interpreting the Result

Primary scoring and prediction matrices are indexed by `site_key`. Site-level
tables that materialize row identity include both `site_key` and `display_id`.

`reference_display_ambiguity_policy="error"` is the default. It rejects a
reference row such as `MAPK14;Y182;` when that display label matches multiple
dataset `site_key` rows. Use `"allow_with_diagnostics"` only when projecting
the same display-level evidence to every matched `site_key` row is intended.
Diagnostics include the display ID and matched `site_key` values. This policy
does not collapse duplicate display labels.

Kinase Library-style workflow scores are normalized support scores for
within-run ranking. They are not official Kinase Library predictor parity and
not calibrated probabilities.

KSEA-style and ssGSEA-style activity methods are explicit PhosPy activity
summaries. KSEA-style activity is not equivalent to PhosR kinase activity
inference. ssGSEA-style activity is not PTM-SEA support.

## Provenance and Reproducibility

Workflow provenance records resolved references, scoring/prediction/activity
configuration, scientific policy records, table fingerprints, and workflow
diagnostics. Adaptive prediction requires `random_state` when
`mode="adaptive_ensemble"` so reruns can be audited.

## Limitations

- Bundled runtime references are rat-first in this release.
- Scores are relative support values, not calibrated causal inference.
- Kinase Library-style scoring requires a compatible caller-supplied local
  resource and does not silently fall back when resource lanes are incompatible.
- Activity is optional and method-specific.
- No broad PhosR activity equivalence is claimed.

## Minimal Example

```python
from phospy import KinaseWorkflow
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    ReferencePreset,
)

request = KinaseWorkflowRequest(
    dataset=dataset,
    references=ReferencePreset.AUTO,
    scoring_config=KinaseScoringConfig(
        min_substrates=2,
        include_diagnostic_scoring_tables=False,
    ),
    prediction_config=KinasePredictionConfig.deterministic(),
    activity_config=None,
)

kinase_result = KinaseWorkflow().run(request)

scores = kinase_result.scoring_result.authoritative_scores
predictions = kinase_result.prediction_result.pred_mat
print(scores.shape)
print(predictions.head())
```
