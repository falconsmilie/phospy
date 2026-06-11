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

## Localisation Prerequisite

This workflow expects a site-level dataset that already passed explicit
localisation-confidence policy at dataset build time. A typical strict policy:

- metadata column: `site_metadata["localisation_confidence"]`
- minimum threshold: `0.75`
- failure behaviour: dataset build fails if the column is missing, has
  invalid/missing values, or contains below-threshold values
- why it matters: low-confidence localisation can mis-assign kinase-substrate
  interpretation at site level

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

Minimal request snippets below focus on kinase API shape and assume this
upstream dataset-localisation policy is already configured.

## Site Identity and Reference Matching

`KinaseWorkflow` operates on `site_key` rows from the
`AnalysisReadyPhosphoDataset`. Reference tables may still use display IDs:
`kinase_substrate_map.substrate_site` and `site_sequences.index` are matched to
dataset `site_metadata.display_id` and then projected through an explicit
`display_id` -> `site_key` mapping before scoring.

`display_id` is metadata and may repeat. Duplicate `display_id` values remain
valid when the corresponding `site_key` values differ, and the workflow does
not collapse or deduplicate rows by display label.

Display-level kinase-substrate reference matching can be ambiguous. If one
reference substrate display ID matches multiple dataset `site_key` rows, the
workflow rejects the projection by default because the reference row does not
identify which protein-scoped phosphosite row carries the evidence.

Set `reference_display_ambiguity_policy="allow_with_diagnostics"` only when
you intentionally want the same display-level kinase-substrate evidence
projected to every matching `site_key` row. Opt-in diagnostics include:

- `display_id`
- matched `site_key` values
- matched row count
- affected reference row position/index, kinase, and substrate entry

Example ambiguity:

```python
import pandas as pd

from phospy.api import ReferenceBundle

dataset.site_metadata.loc[
    :,
    ["site_key", "display_id", "protein_identifier"],
]
# site_key for MAPK14_A, display_id MAPK14;Y182;
# site_key for MAPK14_B, display_id MAPK14;Y182;

references = ReferenceBundle(
    organism=dataset.organism,
    kinase_substrate_map=pd.DataFrame(
        {
            "kinase": ["MAP2K6"],
            "substrate_site": ["MAPK14;Y182;"],
        }
    ),
    site_sequences=pd.DataFrame(
        {"site_sequence": ["AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    ),
)
```

The reference row says only `MAPK14;Y182;`; it does not say whether the evidence
belongs to the `MAPK14_A` or `MAPK14_B` `site_key`. With the default
`reference_display_ambiguity_policy="error"`, `KinaseWorkflow.run(...)` fails
before scoring. With `"allow_with_diagnostics"`, the workflow projects that one
display-level reference row to both matching `site_key` rows and records the
one-display-to-many-site diagnostic. In both modes, prediction and scoring
outputs remain indexed by `site_key`.

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

`KinaseWorkflowRequest` is a lightweight command payload. Constructing it does
not prove that the dataset, references, configs, localisation metadata, or
reference-projection policy are scientifically valid. `KinaseWorkflow.run(...)`
validates the request before interpretation, scoring, prediction, or activity
execution.

| Parameter | Type | Default | Required | How to Use It |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | None | Yes | Dataset returned by `AnalysisReadyDatasetBuilder.run(...)`. It must already be strict and missing-value-free. |
| `references` | `ReferencePreset` or `ReferenceBundle` | `ReferencePreset.AUTO` | No | Reference source for kinase-substrate and site-sequence support. Use `ReferencePreset.AUTO` with the bundled rat beginner lane, or pass a custom `ReferenceBundle`. |
| `scoring_config` | `KinaseScoringConfig` | `KinaseScoringConfig()` | No | Controls kinase scoring support thresholds, diagnostics, missing-value profile handling, and mixed total-protein guardrails. |
| `prediction_config` | `KinasePredictionConfig` | `KinasePredictionConfig()` | No | Controls deterministic or adaptive kinase prediction. |
| `activity_config` | `KinaseActivityConfig` or `None` | `KinaseActivityConfig()` | No | Controls optional activity output. Use `None` to skip activity entirely. |
| `site_sequence_conflict_policy` | `str` | `"prefer_reference"` | No | Controls dataset-vs-reference site-sequence conflicts. Allowed values are `"prefer_reference"`, `"prefer_dataset"`, and `"error"`. |
| `reference_display_ambiguity_policy` | `str` | `"error"` | No | Controls one-display-to-many-`site_key` kinase-substrate reference projection. Allowed values are `"error"` and `"allow_with_diagnostics"`. |

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

## Reference Display Ambiguity Policy

| Value | Behaviour | When to Use It |
| --- | --- | --- |
| `"error"` | Fails when a display-level kinase-substrate reference row maps to multiple dataset `site_key` rows. | Default. Use when the reference table lacks protein-scoped substrate identity. |
| `"allow_with_diagnostics"` | Preserves one-to-many projection and records structured diagnostics. | Use only after deciding that projecting the same reference evidence to all matching protein-scoped rows is scientifically intended. |

```python
request = KinaseWorkflowRequest(
    dataset=dataset,
    references=custom_references,
    reference_display_ambiguity_policy="allow_with_diagnostics",
)
```

This policy is kinase-reference projection specific. It does not make
`display_id` unique and does not change dataset identity validation. This policy
does not collapse duplicate display labels.

After opt-in projection, audit the ambiguity diagnostics before interpreting
the result:

```python
assert kinase_result.provenance is not None
scoring_diagnostics = kinase_result.provenance.workflow_parameters[
    "scoring_diagnostics"
]
diagnostics = scoring_diagnostics["site_sequence_merge"][
    "display_reference_matching"
]
print(diagnostics["one_to_many_display_reference_matches"])
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

Primary scoring and prediction matrices are indexed by `site_key`. Site-level
tables that materialize row identity, such as `substrate_list` and activity
target tables, include both `site_key` and `display_id` for audit and display.

### Result Construction Contract

`KinaseWorkflowResult` is a workflow-owned container. Direct construction is
intentionally minimal and does not re-run reference compatibility, dataset
alignment, scoring, prediction, activity, attrition, or provenance validation.
Use `KinaseWorkflow.run(...)` when you need a scientifically coherent public
kinase result.

The nested stage results are the directly constructible table contracts:
`KinaseScoringResult`, `KinasePredictionResult`, and `KinaseActivityResult`
validate their own public table schemas and defensive DataFrame ownership.
Cross-table and cross-object identity coherence is guaranteed for
workflow-created results, not arbitrary hand-assembled containers.

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
