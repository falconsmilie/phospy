# Kinase Workflow

`KinaseWorkflow` scores kinase-substrate evidence, ranks candidate kinase
support, and can optionally compute exploratory kinase activity score tables.
Use it after building an `AnalysisReadyPhosphoDataset`.

## When to Use This Workflow

Use this workflow when your analysis-ready phosphosite dataset has site
sequences, protein-scoped row identity, and suitable kinase reference data.

Good fits:

- rat beginner runs with `ReferencePreset.AUTO`
- custom human or mouse runs with an explicit `ReferenceBundle`
- default PhosR-inspired rank-weighted scoring
- explicit Kinase Library-style motif scoring with a caller-supplied
  `KinaseLibraryResource`
- optional activity-like score summaries from workflow prediction output

Scores are relative support values within a run, not calibrated probabilities
or proof of causal regulation.
Kinase activity score outputs are substrate/reference-dependent exploratory
summaries. Sparse substrate coverage weakens interpretation, and causal kinase
activity claims require external validation.
The default `phosr_rank_weighted` value names PhosPy's PhosR-inspired
rank-weighted scoring mode. It is not an exact PhosR implementation and is not
intended to provide numerical parity with PhosR.

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
| `activity_config` | `KinaseActivityConfig`, `None`, or disabled config for optional activity score output. |
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
| `scoring_mode` | `"phosr_rank_weighted"` | Supported modes: `"phosr_rank_weighted"`, `"kinase_library_motif"`, `"combined_profile_motif"`. The default is PhosR-inspired PhosPy scoring, not a PhosR compatibility mode. |
| `min_substrates` | `2` | Minimum unique usable substrates for kinase scoring support. The public floor is `2`. |
| `include_diagnostic_scoring_tables` | `False` | Adds non-primary diagnostic scoring tables. |
| `include_substrate_contributions` | `False` | Adds an optional substrate-level contribution table to `KinaseWorkflowResult`. |
| `profile_missing_value_strategy` | `"strict"` | Use `"median_skipna"` only when skipping missing profile values is intended. |
| `localisation_requirement` | `LocalisationRequirement()` | Workflow-level localisation requirement. |
| `allow_mixed_total_protein_quantitative_meaning` | `False` | Keep `False` unless mixed corrected/uncorrected rows are intended. |

### Minimum Substrate Support

Kinase scoring uses `KinaseScoringConfig.min_substrates` to decide whether a
kinase has enough profile support to enter scoring. The default is `2`, and the
public config cannot be set below `2`.

A usable scoring substrate is a unique kinase-substrate reference entry that:

- resolves through reference projection to a dataset `site_key`
- remains in the workflow scoring phospho matrix
- has usable site-sequence support for kinase scoring

Duplicate map rows for the same kinase/site count once. Reference substrates
that are unmapped, absent from the analysis-ready dataset, missing sequence
support, or filtered before scoring do not count toward `min_substrates`.
Because `AnalysisReadyPhosphoDataset` is the public boundary, quantified
substrates are expected to come from numeric, missing-value-free phosphosite
rows.

Kinases below `min_substrates` are excluded from profile scoring and downstream
score columns. If no kinase reaches the floor, the workflow fails early with
`seam=kinase.interpreter.eligible_kinases` and reports overlap/support counts.
Kinases exactly at the floor are included, but they are minimally supported and
should be interpreted cautiously. PhosPy does not add a separate
`weakly_supported` result column; use the support threshold, eligibility counts,
and provenance to audit these cases.

Single-substrate profiles are not part of the public scoring contract. A profile
based on one quantified substrate can be dominated by that substrate's own
sample pattern, making profile-correlation support fragile.

Diagnostics to check:

- `kinase_result.eligibility_report.eligible_kinases`
- `kinase_result.eligibility_report.excluded_kinases_below_min_substrates`
- `kinase_result.site_attrition_summary.scoring`
- `kinase_result.provenance.workflow_parameters["scoring_config"]`

This differs from broad PhosR-style expectations: `phosr_rank_weighted` is a
PhosPy scoring mode with an explicit public support floor, not an exact PhosR
compatibility or numerical parity mode.

For `scoring_mode="phosr_rank_weighted"`, PhosPy uses available
substrate/reference evidence to build kinase profiles, scores profile support,
uses motif support when sequence/reference evidence allow, and combines
available profile and motif evidence with rank-derived weights. The configured
`min_substrates` floor controls which kinases have enough support to enter this
lane.

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
| `min_substrates` | `3` | Weighted activity-like score substrate floor. |
| `top_n_substrates` | `20` | Weighted activity substrate cap. |
| `ksea_min_substrates` | `5` | KSEA-style substrate floor. |
| `ssgsea_min_substrates` | `5` | ssGSEA-style substrate floor. |
| `ssgsea_random_seed` | `0` | Required when permutations are enabled. |

Activity score substrate support is counted from prediction/activity inputs,
not from the scoring profile count. For the simplified weighted and KSEA-style
methods, finite prediction support at or above the configured threshold is
included. KSEA-style and ssGSEA-style results keep not-computable
kinase-condition pairs in `statistics_table` with insufficient-substrate status
when the selected method can report that detail. Activity diagnostics also
expose `method_summary`, `substrate_count_matrix`,
`thresholded_substrate_counts`, and threshold-membership metadata where
supported.

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
| `substrate_contributions` | Optional substrate-level contribution table when enabled. |
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
- `KinaseEligibilityReport.eligible_kinases`
- `KinaseEligibilityReport.excluded_kinases_below_min_substrates`

`activity_result.activity_matrix` is the method-neutral primary activity matrix
for kinase activity scores. Deprecated compatibility aliases such as
`activity_scores` and `weighted_activity` are not preferred for new
documentation or code.

### Substrate Contribution Table

Set `KinaseScoringConfig(include_substrate_contributions=True)` to attach
`kinase_result.substrate_contributions`.

```python
request = KinaseWorkflowRequest(
    dataset=dataset,
    references=references,
    scoring_config=KinaseScoringConfig(
        min_substrates=2,
        include_substrate_contributions=True,
    ),
    activity_config=None,
)

kinase_result = KinaseWorkflow().run(request)
contributions = kinase_result.substrate_contributions
```

The table contains one row per projected kinase-substrate evidence row. It is
off by default because it can be large.

Stable columns:

| Column | Meaning |
| --- | --- |
| `kinase` | Kinase identifier used in scoring. |
| `substrate_site` | Dataset `site_key` used by the workflow. |
| `substrate_identifier` | Display identifier when available. |
| `value_used_in_scoring` | Score value used for the selected score component, or missing when excluded or unavailable. |
| `score_component` | Score lane summarized by the row, such as `rank_weighted_fusion_scores`. |
| `score_source` | More specific evidence source when available. |
| `reference_source_name`, `reference_source_version`, `reference_bundle_id`, `reference_identifier_namespace` | Reference metadata when supplied. |
| `status` | `included` or `excluded`. |
| `exclusion_reason` | Reason an excluded row was not used in the score. |
| `ambiguous` | `True` when display-level reference mapping was one-to-many. |

Use `status`, `exclusion_reason`, and `ambiguous` to audit substrate support.
`value_used_in_scoring` is an evidence-summary value for the selected score
component. It is not a calibrated effect size and does not prove direct
regulation.

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

Default `phosr_rank_weighted` scores are PhosR-inspired rank-weighted support
scores implemented by PhosPy. They are not exact PhosR scores, and small or
large numerical differences from PhosR should not be interpreted as bugs by
themselves.

KSEA-style and ssGSEA-style activity methods are explicit PhosPy activity
score summaries. KSEA-style activity scores are not equivalent to PhosR kinase
activity inference. ssGSEA-style activity-like scores are not PTM-SEA support.

### Activity Interpretation

Keep these result concepts separate:

- Substrate association means a kinase-substrate edge is present in
  `target_table` or a reference/prediction support table.
- Enrichment means a substrate set is over-represented or concentrated by the
  selected KSEA-style or ssGSEA-style scoring rule.
- A kinase activity score or activity-like score is the numeric
  `activity_matrix` output from the selected method.
- Causal kinase activity means a biological activation claim. PhosPy activity
  score outputs do not prove this by themselves.

Interpret activity scores as exploratory unless the study design and external
validation support stronger claims. Scores depend on substrate coverage,
reference evidence, threshold choices, finite phosphosite values, and the
selected method. Missing or sparse substrate support weakens interpretation.

## Provenance and Reproducibility

Workflow provenance records resolved references, scoring/prediction/activity
configuration, scientific policy records, table fingerprints, and workflow
diagnostics. Adaptive prediction requires `random_state` when
`mode="adaptive_ensemble"` so reruns can be audited.

## Limitations

- Bundled runtime references are rat-first in this release.
- Scores are relative support values, not calibrated causal inference.
- `phosr_rank_weighted` is PhosR-inspired PhosPy scoring, not exact PhosR
  numerical compatibility.
- Kinase Library-style scoring requires a compatible caller-supplied local
  resource and does not silently fall back when resource lanes are incompatible.
- Activity score output is optional, method-specific, and exploratory unless
  the study design supports stronger claims.
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
