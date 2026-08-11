# Kinase Analysis

`KinaseWorkflow` ranks kinase–substrate support for each phosphosite in an
`AnalysisReadyPhosphoDataset`.

Use this workflow when you have phosphosite intensities, centered site
sequences, and suitable reference data, and you want to explore which kinases
have the strongest support for each site. The result includes score and
prediction matrices, optional activity-like summaries, diagnostics, provenance,
and caveats.

!!! note "Interpret scores carefully"
    Kinase scores show **relative support within a run**. They are not
    calibrated probabilities, not calibrated causal inference, and not direct
    proof of kinase activation.

## At a Glance

| Question | Answer |
| --- | --- |
| **Input** | An analysis-ready phosphosite dataset, reference data, and explicit scoring settings. |
| **Main request** | `KinaseWorkflowRequest` |
| **Run** | `KinaseWorkflow().run(request)` |
| **Main result** | `KinaseWorkflowResult` |
| **Primary score output** | `result.scoring_result.authoritative_scores` |
| **Primary prediction output** | `result.prediction_result.pred_mat` |
| **Optional activity output** | `result.activity_result.activity_matrix` |

## Before You Start

The dataset must provide:

- unique, protein-scoped `site_key` row identities;
- numeric phosphosite values on a scale accepted by the selected method;
- required identity metadata, including `display_id`, `organism`,
  `protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`;
- a usable, centered `site_sequence` for each scored site;
- localisation evidence when required by the scoring policy; and
- a reference context that is compatible with the dataset.

Configure localisation while building the dataset. With the policy below, the
build fails when localisation is missing, invalid, absent for a row, or below
`min_confidence`.

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

Bundled runtime references are rat-only for `ReferencePreset.AUTO` in this
release. Human, mouse, and custom analyses require an explicit
`ReferenceBundle`.

Reference rows are matched through `display_id` and projected back to dataset
`site_key` identity. The default
`reference_display_ambiguity_policy="error"` rejects a display label that
matches more than one site. With `"allow_with_diagnostics"`, diagnostics record
which reference display labels matched multiple `site_key` values. PhosPy does
not collapse duplicate display labels.

## Example

The following example builds a small rat dataset and runs deterministic kinase
prediction. The values are illustrative rather than a biological benchmark.

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.advanced import (
    KinasePredictionConfig,
    KinaseReliabilityProfile,
    KinaseScoringConfig,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    IntensityScaleKind,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
)

phospho = pd.DataFrame(
    {
        "control_1": [8200.0, 9100.0, 6000.0],
        "control_2": [8000.0, 9000.0, 5900.0],
        "treated_1": [16200.0, 9150.0, 13000.0],
    },
    index=["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3A", "TSC2"],
        "site": ["Y182", "S21", "S939"],
        "site_sequence": [
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
            "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
        ],
        "protein_identifier": ["MAPK14", "GSK3A", "TSC2"],
        "protein_group_id": ["MAPK14", "GSK3A", "TSC2"],
        "localisation_confidence": [0.95, 0.94, 0.96],
    },
    index=phospho.index,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale=IntensityScaleKind.LINEAR,
        preprocessing_config=DatasetPreprocessingConfig(
            localisation=DatasetLocalisationConfig(
                mode="require_threshold",
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            )
        ),
    )
)

request = KinaseWorkflowRequest(
    dataset=dataset,
    references=ReferencePreset.AUTO,
    scoring_config=KinaseScoringConfig(
        reliability_profile=KinaseReliabilityProfile.CUSTOM,
        min_substrates=2,
    ),
    prediction_config=KinasePredictionConfig.deterministic(),
    activity_config=None,
)

result = KinaseWorkflow().run(request)

scores = result.scoring_result.authoritative_scores
predictions = result.prediction_result.pred_mat
print(scores.head())
print(predictions.head())
```

## Request

Create a `KinaseWorkflowRequest`.

| Parameter | Type | Required or Default | Description | Main Constraint |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | Required | Dataset to score. | Must use `site_key` rows and provide usable `site_sequence` values. |
| `references` | `ReferencePreset` or `ReferenceBundle` | `ReferencePreset.AUTO` | Kinase–substrate and sequence references. | `AUTO` resolves rat-only bundled data. Other organisms require an explicit bundle. |
| `scoring_config` | `KinaseScoringConfig` or `None` | `None` | Scoring policy. | Execution requires an explicit config with `reliability_profile` set. |
| `prediction_config` | `KinasePredictionConfig` | `KinasePredictionConfig()` | Candidate-ranking policy. | Adaptive mode requires a non-negative `random_state`. |
| `activity_config` | `KinaseActivityConfig` or `None` | `None` | Optional activity-like substrate summary. | `None` disables activity. |
| `site_sequence_conflict_policy` | `"prefer_reference"`, `"prefer_dataset"`, or `"error"` | `"prefer_reference"` | Resolves dataset/reference sequence conflicts. | `"error"` stops the run when sequences conflict. |
| `reference_display_ambiguity_policy` | `"error"` or `"allow_with_diagnostics"` | `"error"` | Handles a reference display label that maps to multiple dataset rows. | The diagnostic policy preserves separate `site_key` rows. |
| `kinase_library_resource` | `KinaseLibraryResource` or `None` | `None` | Optional local motif resource. | Required by Kinase Library-style scoring modes; official data are not bundled. |

<details markdown="1">
<summary><strong>Scoring Parameters</strong></summary>

### `KinaseScoringConfig`

| Parameter | Default | Description | Main Constraint |
| --- | --- | --- | --- |
| `reliability_profile` | Required | Declares `"exploratory"`, `"production"`, or `"custom"` reliability intent. | A bare `KinaseScoringConfig()` is rejected. Prefer `.exploratory()` or `.production(...)` when appropriate. |
| `min_substrates` | `2` | Minimum substrate support for profile scoring. | Must be at least 2. |
| `scoring_mode` | `"phosr_rank_weighted"` | Selects profile, motif, or combined scoring. | Also supports `"kinase_library_contextual_motif"`, `"kinase_library_motif_only"`, and `"combined_profile_motif"`. |
| `include_diagnostic_scoring_tables` | `False` | Adds method diagnostics. | Increases result size. |
| `include_substrate_contributions` | `False` | Adds per-substrate contribution records. | Useful for audit and interpretation; increases result size. |
| `profile_missing_value_strategy` | `"strict"` | Controls profile construction when values are missing. | `"median_skipna"` is the alternative. Scoring methods do not impute values. |
| `profile_self_inclusion_policy` | `"allow"` | Controls whether a site contributes to its own profile score. | Production presets use `"leave_one_out"`. |
| `attrition_policy` | `KinaseAttritionPolicy()` | Sets minimum retained fractions. | Violations warn or fail according to `on_violation`. |
| `localisation_requirement` | `LocalisationRequirement()` | Sets workflow-level localisation requirements. | Production site-level settings require present evidence with probability at least 0.75. |
| `reference_context_compatibility_policy` | `"require_known_match"` | Handles unknown dataset/reference context. | `"allow_unknown_with_caveat"` records a caveat and should be deliberate. |
| `allow_mixed_total_protein_quantitative_meaning` | `False` | Allows mixed corrected and uncorrected quantitative meaning. | Keep `False` unless the mixture is intentional and documented. |

### `KinaseAttritionPolicy`

| Parameter | Default | Description |
| --- | --- | --- |
| `minimum_reference_overlap_fraction` | `0.0` | Minimum fraction of dataset sites with reference overlap. |
| `minimum_sequence_supported_fraction` | `0.0` | Minimum fraction with usable sequence support. |
| `minimum_scored_fraction` | `0.0` | Minimum fraction contributing to the final score matrix. |
| `on_violation` | `"warn"` | Use `"error"` to fail instead of recording a warning. |

</details>

<details markdown="1">
<summary><strong>Prediction and Activity Parameters</strong></summary>

### `KinasePredictionConfig`

| Parameter | Default | Description | Main Constraint |
| --- | --- | --- | --- |
| `top_k` | `30` | Number of top substrate candidates retained per kinase. | Must be at least 1. |
| `deterministic_max_selected_kinases` | `10` | Maximum kinases retained by deterministic ranking. | Must be at least 1. |
| `adaptive_ensemble_runs` | `10` | Ensemble runs in adaptive mode. | Must be at least 1. |
| `mode` | `"deterministic_ranking"` | Selects deterministic or `"adaptive_ensemble"` prediction. | Adaptive mode requires `random_state`. |
| `adaptive_policy` | `"stable"` | Adaptive policy; `"r_parity"` is also supported. | Used only in adaptive mode. |
| `n_iterations` | `5` | Iteration count for adaptive logic. | Must be at least 1. |
| `random_state` | `None` | Seed for adaptive reproducibility. | Must be non-negative when supplied. |

### `KinaseActivityConfig`

| Parameter | Default | Description | Main Constraint |
| --- | --- | --- | --- |
| `enabled` | `True` | Enables activity when the config is supplied. | `activity_config=None` disables the stage first. |
| `method` | `"simplified_weighted_substrate_activity"` | Selects the activity-like method. | Also supports `"ksea_zscore"` and `"ssgsea_substrate_enrichment"`. |
| `threshold` | `0.6` | Prediction threshold for simplified weighted activity. | Must be between 0 and 1. |
| `min_substrates` | `3` | Minimum substrates for simplified weighted activity. | Must meet the method floor. |
| `top_n_substrates` | `20` | Maximum substrates used by simplified weighted activity. | Must meet the method floor. |
| `ksea_min_substrates` | `5` | Minimum substrates for KSEA-style output. | Must meet the method floor. |
| `ksea_evidence_threshold` | `None` | Optional KSEA membership threshold. | When set, must be between 0 and 1. |
| `ksea_p_value_method` | `"normal_approximation"` | KSEA-style *p*-value method. | Ordinary *p* and *q* values require independent membership provenance. |
| `ksea_adjust_p_values` | `True` | Applies Benjamini–Hochberg adjustment when eligible. | Adjustment is performed per profile. |
| `ssgsea_min_substrates` | `5` | Minimum substrates for ssGSEA-style output. | Must meet the method floor. |
| `ssgsea_ranking_direction` | `"descending"` | Ranking direction for contrast/effect input. | `"ascending"` is also supported. |
| `ssgsea_permutations` | `0` | Substrate-label permutations for empirical *p* values. | Must be non-negative. |
| `ssgsea_random_seed` | `0` | Seed for permutations. | Required when permutations are greater than 0. |
| `ssgsea_adjust_p_values` | `True` | Applies Benjamini–Hochberg adjustment to permutation *p* values. | Adjustment is performed per profile. |

For ssGSEA-style permutations, `profile_id` identifies the quantitative profile.
The seed policy `stable_by_method_profile_kinase` keeps the permutation stream
stable for the method, profile, and kinase.

</details>

<details markdown="1">
<summary><strong>Quantitative Method Contract</strong></summary>

The selected method determines the accepted intensity scale, quantitative
meaning, profile axis, missing-value behavior, and statistical interpretation.

| Method | Accepted scale | Accepted meaning | Required centring/standardisation | Missing values | Profile axis | Statistical interpretation | P-value interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| kinase_scoring.combined_profile_motif | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| kinase_scoring.kinase_library_contextual_motif | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| kinase_scoring.kinase_library_motif_only | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio, unknown | Requires centered phosphosite sequence context; quantitative centring is not applicable because phospho values are not consumed by motif-only scoring.; No quantitative standardisation is required or performed for motif-only scoring. | Phospho missing values are not read by motif-only scoring; no missing-value transformation or imputation is performed. | Rows are protein-scoped site_key phosphosites with centered sequence context; quantitative columns are not used for motif-only score calculation. | Scores are sequence-motif support scores from the supplied Kinase Library-style resource, not abundance-profile statistics. | none |
| kinase_scoring.phosr_rank_weighted | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| ksea_zscore_v1 | log2 | phosphosite_log_abundance, phospho_total_log_ratio, contrast_log2_fold_change, differential_effect_size | Uses log2 sample, total-corrected, contrast, or effect profiles as declared by the dataset; no centring is performed in the method.; Requires log2 abundance, log2 total-corrected ratio, log2 contrast fold-change, or pre-standardised effect semantics; linear raw abundance is rejected. | Finite values define per-profile substrate and background sets; missing values are omitted from those calculations without imputation. | Columns must represent log-scale sample profiles, contrasts, or standardised effect profiles; linear raw samples are rejected. | Unweighted substrate-set z-score enrichment over declared log-scale sample, contrast, or effect values with background variance checks. | Two-sided normal-approximation p-values for computed z-scores; available only when typed substrate-membership provenance declares the membership independent of the tested quantitative matrix. Eligible p-values use Benjamini-Hochberg q-value adjustment per profile when enabled; adaptive membership reports descriptive z-scores with p/q unavailable. |
| simplified_weighted_substrate_activity_v1 | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | No method-level centring; activity values are weighted means on the declared input scale.; No automatic standardisation; linear and log2 abundance summaries have different meanings. | Missing substrate values are ignored per profile when computing weighted and thresholded means; no imputation is performed. | Columns must represent sample-level abundance or explicit condition-summary abundance profiles. | Heuristic substrate-supported weighted mean; not a statistical enrichment test and not causal kinase activity proof. | none |
| ssgsea_substrate_enrichment_activity_v1 | log2 | contrast_log2_fold_change, differential_effect_size | Uses ranked contrast/effect values supplied by the caller; no centring is performed inside the method.; Requires log2 contrast fold-change or pre-standardised effect semantics; raw abundance is rejected. | Only finite effect values enter the ranked background; missing values are omitted without imputation. | Columns must represent contrasts or standardised effect profiles, not raw samples. | Rank-walk substrate-set enrichment summary over ordered effect values. Equal-valued finite sites are handled inside the method as tie blocks using the documented block-expectation policy, not row order or lexical site labels. Not PTM-SEA parity and not causal kinase activity proof. | No p-values are produced unless seeded permutations are requested; permutation p-values are two-sided empirical substrate-label permutation p-values, with Benjamini-Hochberg q-values per profile when enabled. |

</details>

## Run the Workflow

```python
result = KinaseWorkflow().run(request)
```

Deterministic prediction is reproducible for the same dataset, references,
configuration, and package version. Adaptive ensemble prediction also requires
an explicit `random_state`.

The workflow raises `WorkflowValidationError` before scoring when required
metadata are missing, references are incompatible, sequence or localisation
requirements fail, quantitative semantics are unsupported, reference projection
is ambiguous under the selected policy, or activity settings do not match the
input.

## Response

`KinaseWorkflow.run(...)` returns a `KinaseWorkflowResult`.

### Main Outputs

| Attribute | Format | Meaning |
| --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | The analysis-ready input carried into the result. |
| `references` | `ReferenceBundle` | The resolved reference bundle used by the run. |
| `scoring_result` | `KinaseScoringResult` | Score matrices, source labels, and scoring diagnostics. |
| `scoring_result.authoritative_scores` | `pandas.DataFrame` | Site-by-kinase score matrix used by downstream prediction and signalome analysis. |
| `prediction_result` | `KinasePredictionResult` | Prediction matrix and optional substrate list. |
| `prediction_result.pred_mat` | `pandas.DataFrame` | Site-by-kinase prediction support matrix. |
| `activity_result` | `KinaseActivityResult` or `None` | Optional activity-like summaries. |
| `eligibility_report` | `KinaseEligibilityReport` or `None` | Counts retained or excluded before scoring. |
| `site_attrition_summary` | Typed attrition summary or `None` | Site counts across preprocessing and scoring. |
| `attrition_provenance` | Typed provenance or `None` | Attrition policy, outcome, warnings, and violations. |
| `substrate_contributions` | `pandas.DataFrame` or `None` | Per-substrate contributions when requested. |
| `input_dataset_preprocessing_report` | `DatasetPreprocessingReport` or `None` | Preprocessing report carried from `dataset`. |
| `provenance` | `RunProvenance` or `None` | Resolved settings, references, method contracts, and input fingerprints. |
| `caveats` | `tuple[ResultCaveat, ...]` | Structured scientific limitations and warnings. |

Score and prediction matrices use unique `site_key` rows and kinase columns.
A missing cell means that the method could not produce a score for that pair; it
does not establish biological absence.

### Activity Output

When activity is enabled, `activity_result.activity_matrix` is the primary
activity matrix for user code.
<!-- phospy-deprecation-compat: activities.result.activity_scores activities.result.weighted_activity -->
`activity_scores` and `weighted_activity` are compatibility aliases and are not
preferred.

| Attribute | Availability | Meaning |
| --- | --- | --- |
| `activity_matrix` | Always when activity succeeds | Method-neutral kinase-by-profile activity score matrix. |
| `substrate_count_matrix` | Always when activity succeeds | Substrate count used for each kinase and profile. |
| `p_value_matrix`, `q_value_matrix` | Method-dependent | Inferential output only when the selected method and provenance permit it. |
| `confidence_interval_low`, `confidence_interval_high` | Method-dependent | Confidence bounds when produced. |
| `statistics_table` | Optional | Long-format method statistics and computability status. |
| `method_diagnostics`, `policy_provenance` | Optional | Method decisions, eligibility, and scientific policy records. |
| `input_semantics`, `profile_metadata`, `membership_selection` | Optional | Quantitative axis and membership evidence used by the method. |

<details markdown="1">
<summary><strong>Detailed Result Contract</strong></summary>

### `KinaseScoringResult`

| Attribute | Availability | Meaning |
| --- | --- | --- |
| `profile_scores` | Always | Profile-derived support scores. |
| `motif_scores` | Optional | Motif-derived scores. |
| `rank_weighted_fusion_scores` | Optional | Fused profile/motif scores. |
| `kinase_library_motif_scores` | Optional | Kinase Library-style motif scores. |
| `combined_profile_motif_scores` | Optional | Combined profile and motif scores. |
| `score_fusion_weights` | Optional | Weights used for score fusion. |
| `score_source_matrix`, `score_source_summary` | Optional | Source category by score and summary counts. |
| `profile_score_diagnostics` | Optional | Leave-one-out and substrate-support diagnostics. |
| `kinase_library_site_diagnostics`, `kinase_library_kinase_diagnostics` | Optional | Resource-specific site and kinase diagnostics. |
| `score_source` | Always | Label for the selected downstream score source. |
| `authoritative_scores` | Always | Score matrix consumed by prediction and signalome analysis. |

### `KinasePredictionResult`

| Attribute or Helper | Availability | Meaning |
| --- | --- | --- |
| `pred_mat`, `to_dataframe()` | Always | Prediction matrix and a defensive snapshot. |
| `substrate_list`, `substrate_list_dataframe()` | Optional | Long-format predicted substrates when produced. |

Typical substrate-list fields are `kinase`, `substrate_site`, `score`, `rank`,
`site_key`, and `display_id`.

### Contribution and Diagnostic Tables

When requested, `substrate_contributions` records the kinase, substrate identity,
value used in scoring, score component and source, reference identity, inclusion
status, exclusion reason, and ambiguity status.

`score_source_summary` separates fused evidence, profile-only evidence, missing
motif overlap, and unavailable scores. `profile_score_diagnostics` reports the
substrate counts before and after leave-one-out handling, the configured floor,
and the reason a pair was not scored.

The eligibility and attrition records report dataset input rows, sequence and
localisation eligibility, reference overlap, final scored sites, activity-stage
contribution, configured minimum fractions, and whether the policy passed,
warned, or failed.

</details>

## Interpret the Result

Higher values indicate stronger support under the selected reference data,
sequence context, quantitative scale, and method. Compare values within the same
run and method; do not treat them as transferable probabilities.

`"phosr_rank_weighted"` is a PhosR-inspired PhosPy implementation, not an exact
PhosR numerical compatibility mode. Kinase Library-style modes require a
compatible caller-supplied resource and do not silently fall back when that
resource is missing or incompatible.

Activity output depends on substrate coverage and membership provenance.
KSEA-style *p* and *q* values are available only when typed provenance shows
that substrate membership was selected independently of the tested matrix.
ssGSEA-style output is not PTM-SEA parity. Causal claims require an appropriate
study design and external validation.

## Common Issues

| Issue | What to Check |
| --- | --- |
| `ReferencePreset.AUTO` fails. | Bundled references are rat-only. Use `Organism.RAT` or pass an explicit `ReferenceBundle`. |
| Sequence validation fails. | Confirm that `site_sequence` is present, centered, and compatible with the selected scoring mode. |
| Localisation validation fails. | Add localisation metadata and configure `DatasetLocalisationConfig`; required localisation fails when missing. |
| A bare scoring config fails. | Set `reliability_profile` or use an exploratory or production preset. |
| Reference projection is ambiguous. | Keep `"error"`, or use `"allow_with_diagnostics"` only when repeated display labels intentionally represent separate `site_key` rows. |
| Scores are sparse or empty. | Review reference overlap, sequence eligibility, attrition, score-source diagnostics, and `min_substrates`. |
| Activity is absent. | Confirm that `activity_config` is supplied and enabled, and that the selected method has sufficient eligible substrates. |

## Related Guides

- [Prepare a dataset](dataset-build-workflow.md)
- [Reference data](../reference_bundles.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)
- [Signalome analysis](signalome.md)
