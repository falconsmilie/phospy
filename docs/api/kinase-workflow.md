# Kinase Workflow

This page explains the kinase workflow API. Use it after you have an
`AnalysisReadyPhosphoDataset` from the dataset builder.

## Purpose

`KinaseWorkflow` resolves references, scores kinase-substrate evidence, predicts
candidate kinase regulation, and can optionally compute kinase activity tables.

```python
kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
    )
)
```

## Imports

```python
from phospy import KinaseWorkflow
from phospy.api import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    ReferenceBundle,
    ReferencePreset,
)
```

## Request Parameters

| Parameter | Type | Default | Required | How to Use It |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | None | Yes | Dataset returned by `AnalysisReadyDatasetBuilder.run(...)`. It must already be strict and missing-value-free. |
| `references` | `ReferencePreset` or `ReferenceBundle` | `ReferencePreset.AUTO` | No | Reference source for kinase-substrate and site-sequence support. Use `ReferencePreset.AUTO` with the bundled rat beginner lane, or pass a custom `ReferenceBundle`. |
| `scoring_config` | `KinaseScoringConfig` | `KinaseScoringConfig()` | No | Controls kinase scoring support thresholds, diagnostics, missing-value profile handling, and mixed total-protein guardrails. |
| `prediction_config` | `KinasePredictionConfig` | `KinasePredictionConfig()` | No | Controls deterministic or adaptive kinase prediction. |
| `activity_config` | `KinaseActivityConfig` or `None` | `KinaseActivityConfig()` | No | Controls optional activity output. Use `None` to skip activity entirely. |
| `site_sequence_conflict_policy` | `str` | `"prefer_reference"` | No | Controls dataset-vs-reference site-sequence conflicts. Allowed values are `"prefer_reference"`, `"prefer_dataset"`, and `"error"`. |

Minimal call:

```python
result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        activity_config=None,
    )
)
```

Custom references:

```python
result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferenceBundle(
            organism=dataset.organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
        ),
    )
)
```

## Site-Sequence Conflict Policy

| Value | Behaviour | When to Use It |
| --- | --- | --- |
| `"prefer_reference"` | Uses the reference sequence when dataset and reference disagree. | Good default when the reference bundle is curated for the analysis. |
| `"prefer_dataset"` | Uses the dataset sequence when a conflict is found. | Useful when the input table has experimentally reviewed or project-specific sequence context. |
| `"error"` | Fails on a conflict. | Best for audits where any disagreement should be resolved before analysis. |

```python
request = KinaseWorkflowRequest(
    dataset=dataset,
    references=ReferencePreset.AUTO,
    site_sequence_conflict_policy="error",
)
```

## Scoring Configuration

Use presets first:

```python
scoring = KinaseScoringConfig.default()
strict_missing = KinaseScoringConfig.strict_missing_values()
```

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `min_substrates` | `int` | `2` | Integer `>= 2` | Minimum quantified substrates required for a kinase to receive scoring support. Raise this for more conservative scoring. |
| `include_diagnostic_scoring_tables` | `bool` | `False` | `True`, `False` | Includes non-primary scoring diagnostics such as motif score and score-fusion helper tables. |
| `profile_missing_value_strategy` | `str` | `"strict"` | `"strict"`, `"median_skipna"` | Controls how profile medians handle missing values in multi-substrate kinase profiles. |
| `allow_mixed_total_protein_quantitative_meaning` | `bool` | `False` | `True`, `False` | Allows kinase scoring on datasets that contain both corrected and uncorrected rows after partial total-protein correction. Keep `False` unless that mixed interpretation is intentional. |

Example with diagnostic scoring tables:

```python
scoring = KinaseScoringConfig(
    min_substrates=3,
    include_diagnostic_scoring_tables=True,
    profile_missing_value_strategy="strict",
    allow_mixed_total_protein_quantitative_meaning=False,
)

result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        scoring_config=scoring,
    )
)
```

Missing-value profile handling:

```python
scoring = KinaseScoringConfig(
    profile_missing_value_strategy="median_skipna",
)
```

Use `"median_skipna"` only when skipping missing profile values is scientifically
acceptable for your dataset and reference coverage.

## Prediction Configuration

Use intent presets for the two common prediction lanes:

```python
deterministic = KinasePredictionConfig.deterministic()
adaptive = KinasePredictionConfig.adaptive_reproducible(random_state=1)
```

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `top_k` | `int` | `30` | Integer `>= 1` | Number of top predicted substrate sites retained per kinase. |
| `deterministic_max_selected_kinases` | `int` | `10` | Integer `>= 1` | Breadth of retained kinases in deterministic ranking mode. |
| `adaptive_ensemble_runs` | `int` | `10` | Integer `>= 1` | Number of ensemble executions in adaptive mode. |
| `mode` | `str` | `"deterministic_ranking"` | `"deterministic_ranking"`, `"adaptive_ensemble"` | Selects the prediction lane. |
| `adaptive_policy` | `str` | `"stable"` | `"stable"`, `"r_parity"` | Selects adaptive sampling behaviour. `"stable"` is the recommended production policy. |
| `n_iterations` | `int` | `5` | Integer `>= 1` | Adaptive sampling iterations. |
| `random_state` | `int` or `None` | `None` | Integer `>= 0` when `mode="adaptive_ensemble"`; otherwise optional | Required for reproducible adaptive prediction. Unused by deterministic mode. |

Deterministic example:

```python
prediction = KinasePredictionConfig(
    mode="deterministic_ranking",
    top_k=30,
    deterministic_max_selected_kinases=10,
)
```

Adaptive reproducible example:

```python
prediction = KinasePredictionConfig(
    mode="adaptive_ensemble",
    adaptive_policy="stable",
    adaptive_ensemble_runs=20,
    n_iterations=5,
    top_k=30,
    random_state=12345,
)
```

Parity-oriented adaptive example:

```python
prediction = KinasePredictionConfig(
    mode="adaptive_ensemble",
    adaptive_policy="r_parity",
    adaptive_ensemble_runs=10,
    n_iterations=5,
    random_state=1,
)
```

## Activity Configuration

Activity can be disabled in either of these ways:

```python
request = KinaseWorkflowRequest(dataset=dataset, activity_config=None)
```

```python
request = KinaseWorkflowRequest(
    dataset=dataset,
    activity_config=KinaseActivityConfig(enabled=False),
)
```

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `enabled` | `bool` | `True` | `True`, `False` | Enables or disables activity after construction. |
| `method` | `str` | `"simplified_weighted_substrate_activity"` | `"simplified_weighted_substrate_activity"`, `"ksea_zscore"` | Selects weighted heuristic activity or KSEA-style z-score activity. |
| `threshold` | `float` | `0.6` | `0.0 <= value <= 1.0` | Prediction-score threshold for selected substrates. |
| `min_substrates` | `int` | `3` | Integer `>= 1` | Minimum selected substrates per kinase for weighted activity. |
| `top_n_substrates` | `int` | `20` | Integer `>= 1` | Top predicted substrates used in weighted activity. |
| `ksea_min_substrates` | `int` | `5` | Integer `>= 1` | Minimum substrates per kinase and condition for KSEA-style scoring. |
| `ksea_evidence_threshold` | `float` or `None` | `None` | `0.0 <= value <= 1.0`, or `None` | KSEA membership threshold. Uses `threshold` when `None`. |
| `ksea_p_value_method` | `str` | `"normal_approximation"` | `"normal_approximation"` | P-value method for KSEA-style statistics. |
| `ksea_adjust_p_values` | `bool` | `True` | `True`, `False` | Applies Benjamini-Hochberg q-values per condition for KSEA-style statistics. |

Weighted activity example:

```python
activity = KinaseActivityConfig(
    enabled=True,
    method="simplified_weighted_substrate_activity",
    threshold=0.6,
    min_substrates=3,
    top_n_substrates=20,
)
```

KSEA-style activity example:

```python
activity = KinaseActivityConfig(
    enabled=True,
    method="ksea_zscore",
    threshold=0.6,
    ksea_min_substrates=5,
    ksea_evidence_threshold=None,
    ksea_p_value_method="normal_approximation",
    ksea_adjust_p_values=True,
)
```

`method="ksea_zscore"` enables KSEA-style z-score activity inference with
per-condition computability statuses and statistics output. It is not equivalent
to PhosR kinase activity inference.

`activity_result.activity_scores` is the method-neutral primary activity matrix:

- for `simplified_weighted_substrate_activity_v1`, values are weighted substrate activity scores
- for `ksea_zscore_v1`, values are KSEA z-scores

## Full Kinase Workflow Example

```python
from phospy import KinaseWorkflow
from phospy.api import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    ReferencePreset,
)

scoring = KinaseScoringConfig(
    min_substrates=2,
    include_diagnostic_scoring_tables=False,
    profile_missing_value_strategy="strict",
)

prediction = KinasePredictionConfig(
    mode="deterministic_ranking",
    top_k=30,
    deterministic_max_selected_kinases=10,
)

activity = KinaseActivityConfig(
    enabled=True,
    method="simplified_weighted_substrate_activity",
    threshold=0.6,
    min_substrates=3,
    top_n_substrates=20,
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=scoring,
        prediction_config=prediction,
        activity_config=activity,
        site_sequence_conflict_policy="prefer_reference",
    )
)

print(kinase_result.scoring_result.profile_scores.shape)
print(kinase_result.prediction_result.pred_mat.round(4))
if kinase_result.activity_result is not None:
    print(kinase_result.activity_result.activity_scores.round(4))
```

## Output Tables

Common result tables include:

```python
profile_scores = kinase_result.scoring_result.profile_scores
ranked_scores = kinase_result.scoring_result.rank_weighted_fusion_scores
prediction_matrix = kinase_result.prediction_result.pred_mat
```

Activity tables are present when activity is enabled:

```python
activity_scores = kinase_result.activity_result.activity_scores
activity_method = kinase_result.activity_result.activity_method
```

The result provenance records active scientific policies and resolved workflow
parameters:

```python
for policy in kinase_result.provenance.scientific_policies:
    print(policy.policy_id)

print(kinase_result.provenance.workflow_parameters["activity_config"])
```
