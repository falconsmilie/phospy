# Kinase analysis workflow

## Plain-language introduction

`KinaseWorkflow` scores kinase-substrate support for an
`AnalysisReadyPhosphoDataset`.

Use it when you have phosphosite intensities with site sequences and want to
rank candidate kinases for each site. The workflow expects a strict dataset,
reference data, explicit scoring settings, and prediction settings. It returns a
`KinaseWorkflowResult` with scoring matrices, a prediction matrix, optional
substrate lists, optional activity summaries, diagnostics, provenance, caveats,
and attrition reports.

Kinase scores are relative support values within a run, not calibrated
probabilities or proof of causal regulation.

## Input and dataset requirements

Start from an `AnalysisReadyPhosphoDataset`; see
[Preparing a dataset](dataset-build-workflow.md).

For kinase analysis, the dataset must provide:

- phosphosite rows keyed by `site_key`;
- numeric phosphosite values on a declared or established scale accepted by the
  selected scoring/activity methods;
- required site metadata, including `site_key`, `display_id`, `organism`,
  `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`, and
  required `site_sequence`;
- centered phosphosite sequence context for current kinase scoring modes;
- site-level localisation evidence when required by `scoring_config`;
- reference context compatible with the dataset.

Configure localisation during dataset building. With this policy, dataset build
fails when localisation metadata is missing, invalid, missing per row, or below
threshold.

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
release. For human, mouse, or custom contexts, pass an explicit
`ReferenceBundle`.

Reference rows are projected through dataset `display_id` metadata into dataset
`site_key` identity. The default
`reference_display_ambiguity_policy="error"` rejects one display label matching
multiple `site_key` rows. With `"allow_with_diagnostics"`, diagnostics record
which reference display labels matched more than one `site_key`. PhosPy does
not collapse duplicate display labels.

## Minimal end-to-end example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.advanced import (
    KinasePredictionConfig,
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
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
        "control_rep1": [8200.0, 9100.0, 6000.0],
        "control_rep2": [8000.0, 9000.0, 5900.0],
        "treatment_rep1": [16200.0, 9150.0, 13000.0],
        "treatment_rep2": [15800.0, 9050.0, 12800.0],
    },
    index=["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3A", "TSC2"],
        "site": ["Y182", "S21", "S939"],
        "site_sequence": [
            ("A" * 15) + "Y" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
        ],
        "display_id": ["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
        "organism": ["rat", "rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id", "protein_id"],
        "protein_identifier": ["MAPK14", "GSK3A", "TSC2"],
        "protein_group_id": ["MAPK14", "GSK3A", "TSC2"],
        "localisation_confidence": [0.95, 0.94, 0.96],
    },
    index=phospho.index.copy(),
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
        include_diagnostic_scoring_tables=False,
        reference_context_compatibility_policy=(
            ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
        ),
    ),
    prediction_config=KinasePredictionConfig.deterministic(),
    activity_config=None,
    site_sequence_conflict_policy="prefer_reference",
)

kinase_result = KinaseWorkflow().run(request)

scores = kinase_result.scoring_result.authoritative_scores
predictions = kinase_result.prediction_result.pred_mat
print(scores.shape)
print(predictions.head())
```

## Request model

Use `KinaseWorkflowRequest`.

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | Required | Analysis-ready phosphosite dataset to score. | Must be `site_key` indexed and carry usable `site_sequence`. |
| `references` | `ReferencePreset | ReferenceBundle` | Default: `ReferencePreset.AUTO` | Reference data used for kinase-substrate support. | `ReferencePreset.AUTO` resolves bundled rat-only references. Other organisms require an explicit `ReferenceBundle`. |
| `scoring_config` | `KinaseScoringConfig | None` | Default: `None` | Scoring policy. | Although the request default is `None`, execution requires an explicit scoring config with `reliability_profile` set. |
| `prediction_config` | `KinasePredictionConfig` | Default: `KinasePredictionConfig()` | Candidate kinase prediction policy. | Adaptive ensemble mode requires a non-negative `random_state`. |
| `activity_config` | `KinaseActivityConfig | None` | Default: `None` | Optional activity-like substrate summary policy. | `None` disables activity. `KinaseActivityConfig(enabled=False)` also disables activity. |
| `site_sequence_conflict_policy` | `"prefer_reference" | "prefer_dataset" | "error"` | Default: `"prefer_reference"` | Resolves dataset/reference sequence conflicts. | `"error"` fails on conflicts; the two preference policies record the selected source. |
| `reference_display_ambiguity_policy` | `"error" | "allow_with_diagnostics"` | Default: `"error"` | Handles reference display labels that map to multiple dataset `site_key` rows. | `"allow_with_diagnostics"` keeps matched rows and records diagnostics; it does not collapse duplicate display labels. |
| `kinase_library_resource` | `KinaseLibraryResource | None` | Default: `None` | Optional local resource for Kinase Library-style motif scoring modes. | Required for Kinase Library-style modes; official Kinase Library data is not bundled. |

`KinaseScoringConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `reliability_profile` | `"exploratory" | "production" | "custom"` | Required keyword | Declares scoring reliability intent. | Use `KinaseScoringConfig.exploratory()` or `.production(...)` presets when possible. Bare `KinaseScoringConfig()` is rejected. |
| `min_substrates` | `int` | Default: `2` | Minimum substrate support for profile scoring. | Must be at least `2`; production preset uses a higher default. |
| `scoring_mode` | string enum | Default: `"phosr_rank_weighted"` | Selects the scoring source/fusion mode. | Supported values: `"phosr_rank_weighted"`, `"kinase_library_contextual_motif"`, `"kinase_library_motif_only"`, and `"combined_profile_motif"`. |
| `include_diagnostic_scoring_tables` | `bool` | Default: `False` | Includes diagnostic scoring tables when available. | Larger outputs. |
| `include_substrate_contributions` | `bool` | Default: `False` | Includes per-substrate contribution table. | Larger outputs; useful for auditing why a kinase score was supported or excluded. |
| `profile_missing_value_strategy` | `"strict" | "median_skipna"` | Default: `"strict"` | Missing-value handling during profile construction. | Missing values are not imputed by scoring methods. |
| `profile_self_inclusion_policy` | `"allow" | "leave_one_out"` | Default: `"allow"` | Controls whether a site can contribute to its own kinase profile score. | Production preset uses leave-one-out. |
| `attrition_policy` | `KinaseAttritionPolicy` | Default: `KinaseAttritionPolicy()` | Minimum retained-fraction policy for scoring attrition. | Violations warn or fail depending on `on_violation`. |
| `localisation_requirement` | `LocalisationRequirement` | Default: `LocalisationRequirement()` | Workflow-level localisation requirement. | Production site-level requirement uses `require_present=True` and `minimum_probability=0.75`. |
| `reference_context_compatibility_policy` | `"require_known_match" | "allow_unknown_with_caveat"` | Default: `"require_known_match"` | Controls unknown dataset/reference context. | Allowing unknown context records a caveat and should be used deliberately. |
| `allow_mixed_total_protein_quantitative_meaning` | `bool` | Default: `False` | Allows mixed total-protein quantitative meaning. | Keep `False` unless the mixture is intentional and documented. |

`KinaseAttritionPolicy`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `minimum_reference_overlap_fraction` | `float` | Default: `0.0` | Minimum fraction of dataset sites with reference overlap. | Must be in `[0, 1]`. |
| `minimum_sequence_supported_fraction` | `float` | Default: `0.0` | Minimum fraction with usable site sequence. | Must be in `[0, 1]`. |
| `minimum_scored_fraction` | `float` | Default: `0.0` | Minimum fraction contributing to final scoring output. | Must be in `[0, 1]`. |
| `on_violation` | `"warn" | "error"` | Default: `"warn"` | Policy action when thresholds are not met. | `"error"` fails before returning a result. |

`LocalisationRequirement`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `require_present` | `bool` | Default: `False` | Requires localisation evidence to be present. | Missing localisation fails when `True`. |
| `minimum_probability` | `float | None` | Default: `None` | Minimum acceptable localisation probability. | If set, must be in `[0, 1]`. Production site-level policy uses `0.75`. |

`KinasePredictionConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `top_k` | `int` | Default: `30` | Maximum top candidate substrates/kinases considered by prediction. | Must be at least `1`. |
| `deterministic_max_selected_kinases` | `int` | Default: `10` | Maximum selected kinases in deterministic ranking mode. | Must be at least `1`. |
| `adaptive_ensemble_runs` | `int` | Default: `10` | Number of runs in adaptive ensemble mode. | Must be at least `1`. |
| `mode` | `"deterministic_ranking" | "adaptive_ensemble"` | Default: `"deterministic_ranking"` | Prediction mode. | Adaptive mode requires `random_state`. |
| `adaptive_policy` | `"stable" | "r_parity"` | Default: `"stable"` | Adaptive ensemble policy. | Applies only in adaptive mode. |
| `n_iterations` | `int` | Default: `5` | Iteration count for supported adaptive logic. | Must be at least `1`. |
| `random_state` | `int | None` | Default: `None` | Seed for adaptive ensemble reproducibility. | Must be non-negative when supplied; required for adaptive mode. |

`KinaseActivityConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `enabled` | `bool` | Default: `True` | Enables activity stage when the config is supplied. | `activity_config=None` disables the stage before this field is considered. |
| `method` | string enum | Default: `"simplified_weighted_substrate_activity"` | Activity-like summary method. | Supported values: `"simplified_weighted_substrate_activity"`, `"ksea_zscore"`, `"ssgsea_substrate_enrichment"`. |
| `threshold` | `float` | Default: `0.6` | Prediction support threshold used by simplified weighted activity. | Must be in `[0, 1]`. |
| `min_substrates` | `int` | Default: `3` | Minimum substrates for simplified weighted activity. | Must meet the method floor. |
| `top_n_substrates` | `int` | Default: `20` | Maximum substrates used by simplified weighted activity. | Must meet the method floor. |
| `ksea_min_substrates` | `int` | Default: `5` | Minimum substrates for KSEA-style activity. | Must meet the method floor. |
| `ksea_evidence_threshold` | `float | None` | Default: `None` | Optional KSEA membership threshold. | If set, must be in `[0, 1]`. |
| `ksea_p_value_method` | `"normal_approximation"` | Default: `"normal_approximation"` | KSEA-style p-value method. | Ordinary p/q values require independent substrate-membership provenance. |
| `ksea_adjust_p_values` | `bool` | Default: `True` | Applies Benjamini-Hochberg q-value adjustment when KSEA p-values are eligible. | Per activity profile. |
| `ssgsea_min_substrates` | `int` | Default: `5` | Minimum substrates for ssGSEA-style activity. | Must meet the method floor. |
| `ssgsea_ranking_direction` | `"descending" | "ascending"` | Default: `"descending"` | Ranking direction for ssGSEA-style activity. | Applies to contrast/effect inputs. |
| `ssgsea_permutations` | `int` | Default: `0` | Number of substrate-label permutations for empirical p-values. | Must be non-negative. |
| `ssgsea_random_seed` | `int | None` | Default: `0` | Seed for ssGSEA-style permutations. | Must be non-negative when supplied; required when permutations are greater than `0`. |
| `ssgsea_adjust_p_values` | `bool` | Default: `True` | Applies Benjamini-Hochberg q-value adjustment when permutation p-values are requested. | Per activity profile. |

Quantitative method contract:

| Method | Accepted scale | Accepted meaning | Required centring/standardisation | Missing values | Profile axis | Statistical interpretation | P-value interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| kinase_scoring.combined_profile_motif | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| kinase_scoring.kinase_library_contextual_motif | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| kinase_scoring.kinase_library_motif_only | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio, unknown | Requires centered phosphosite sequence context; quantitative centring is not applicable because phospho values are not consumed by motif-only scoring.; No quantitative standardisation is required or performed for motif-only scoring. | Phospho missing values are not read by motif-only scoring; no missing-value transformation or imputation is performed. | Rows are protein-scoped site_key phosphosites with centered sequence context; quantitative columns are not used for motif-only score calculation. | Scores are sequence-motif support scores from the supplied Kinase Library-style resource, not abundance-profile statistics. | none |
| kinase_scoring.phosr_rank_weighted | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| ksea_zscore_v1 | log2 | phosphosite_log_abundance, phospho_total_log_ratio, contrast_log2_fold_change, differential_effect_size | Uses log2 sample, total-corrected, contrast, or effect profiles as declared by the dataset; no centring is performed in the method.; Requires log2 abundance, log2 total-corrected ratio, log2 contrast fold-change, or pre-standardised effect semantics; linear raw abundance is rejected. | Finite values define per-profile substrate and background sets; missing values are omitted from those calculations without imputation. | Columns must represent log-scale sample profiles, contrasts, or standardised effect profiles; linear raw samples are rejected. | Unweighted substrate-set z-score enrichment over declared log-scale sample, contrast, or effect values with background variance checks. | Two-sided normal-approximation p-values for computed z-scores; available only when typed substrate-membership provenance declares the membership independent of the tested quantitative matrix. Eligible p-values use Benjamini-Hochberg q-value adjustment per profile when enabled; adaptive membership reports descriptive z-scores with p/q unavailable. |
| simplified_weighted_substrate_activity_v1 | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | No method-level centring; activity values are weighted means on the declared input scale.; No automatic standardisation; linear and log2 abundance summaries have different meanings. | Missing substrate values are ignored per profile when computing weighted and thresholded means; no imputation is performed. | Columns must represent sample-level abundance or explicit condition-summary abundance profiles. | Heuristic substrate-supported weighted mean; not a statistical enrichment test and not causal kinase activity proof. | none |
| ssgsea_substrate_enrichment_activity_v1 | log2 | contrast_log2_fold_change, differential_effect_size | Uses ranked contrast/effect values supplied by the caller; no centring is performed inside the method.; Requires log2 contrast fold-change or pre-standardised effect semantics; raw abundance is rejected. | Only finite effect values enter the ranked background; missing values are omitted without imputation. | Columns must represent contrasts or standardised effect profiles, not raw samples. | Rank-walk substrate-set enrichment summary over ordered effect values. Equal-valued finite sites are handled inside the method as tie blocks using the documented block-expectation policy, not row order or lexical site labels. Not PTM-SEA parity and not causal kinase activity proof. | No p-values are produced unless seeded permutations are requested; permutation p-values are two-sided empirical substrate-label permutation p-values, with Benjamini-Hochberg q-values per profile when enabled. |

## Running the workflow

Call `KinaseWorkflow().run(request)`.

```python
from phospy import KinaseWorkflow
from phospy.api import KinaseWorkflowRequest

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=scoring_config,
    )
)
```

Deterministic prediction mode is deterministic for the same dataset,
references, config, and package version. Adaptive ensemble prediction requires
an explicit seed through `random_state`.

The workflow can raise `WorkflowValidationError` before scoring when required
metadata is missing, references are incompatible, `site_sequence` is missing or
unusable, localisation policy fails, quantitative scale/meaning is unsupported,
reference display projection is ambiguous under the selected policy, or activity
settings are not compatible with the dataset.

## Response model and output formats

`KinaseWorkflow.run(...)` returns `KinaseWorkflowResult`.

Top-level result:

| Attribute or helper | Python type | Always present? | Meaning |
| --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | Yes | Input dataset used by the workflow. |
| `references` | `ReferenceBundle` | Yes | Resolved reference bundle. |
| `scoring_result` | `KinaseScoringResult` | Yes | Scoring matrices and diagnostics. |
| `prediction_result` | `KinasePredictionResult` | Yes | Prediction matrix and optional substrate list. |
| `eligibility_report` | `KinaseEligibilityReport | None` | Optional | Counts of sites/kinases retained or excluded before scoring. |
| `site_attrition_summary` | `KinaseWorkflowPreprocessingAttritionSummary | KinaseWorkflowScoringAttritionSummary | None` | Optional | Site-level attrition counts from preprocessing/scoring. |
| `attrition_provenance` | `KinaseWorkflowAttritionProvenance | None` | Optional | Attrition metrics, policy, outcome, violations, and warnings. |
| `activity_result` | `KinaseActivityResult | None` | Optional | Activity-like substrate summaries when requested and enabled. |
| `provenance` | object or `None` | Optional | Workflow provenance, including resolved config and method contracts. |
| `caveats` | `tuple` | Yes | Structured caveats and warnings. |
| `substrate_contributions` | `pandas.DataFrame | None` | Optional | Per-substrate scoring contribution table when requested. |
| `substrate_contributions_dataframe()` | `pandas.DataFrame | None` | Optional | Defensive snapshot of substrate contributions. |
| `input_dataset_preprocessing_report` | object or `None` | Optional | Preprocessing report carried from the input dataset. |
| `scientifically_equals(...)` | `bool` | Yes | Comparison helper for scientific/result equivalence checks. |

`KinaseScoringResult`:

| Attribute or helper | Python type | Always present? | Meaning |
| --- | --- | --- | --- |
| `profile_scores` | `pandas.DataFrame` | Yes | Profile-derived kinase support scores. |
| `motif_scores` | `pandas.DataFrame | None` | Optional | Motif-derived scores when available. |
| `rank_weighted_fusion_scores` | `pandas.DataFrame | None` | Optional | Fused profile/motif scores for rank-weighted mode. |
| `kinase_library_motif_scores` | `pandas.DataFrame | None` | Optional | Kinase Library-style motif score matrix when used. |
| `combined_profile_motif_scores` | `pandas.DataFrame | None` | Optional | Combined profile/motif matrix when used. |
| `score_fusion_weights` | `pandas.DataFrame | None` | Optional | Fusion weights used when score fusion is active. |
| `score_source_matrix` | `pandas.DataFrame | None` | Optional | Per-site/per-kinase source category for authoritative scores. |
| `score_source_summary` | `pandas.DataFrame | None` | Optional | Counts by score-source category for each kinase. |
| `profile_score_diagnostics` | `pandas.DataFrame | None` | Optional | Leave-one-out/profile scoring diagnostics. |
| `kinase_library_site_diagnostics` | table-like or `None` | Optional | Site-level Kinase Library-style diagnostics. |
| `kinase_library_kinase_diagnostics` | table-like or `None` | Optional | Kinase-level Kinase Library-style diagnostics. |
| `score_source` | `str` | Yes | Label for the score source selected by the config. |
| `authoritative_scores` | `pandas.DataFrame` | Yes | Score matrix used by downstream prediction and signalome steps. |

Score and prediction matrices:

| Column or index | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| index | Dataset phosphosite identity. | Unique `site_key` string labels | Yes |
| columns | Kinase labels from the resolved reference/scoring context. | strings | Yes |
| cells | Support score for a site/kinase pair. | numeric; missing may appear where a score is unavailable | Yes |

Some scoring modes produce unit-interval scores; profile-derived support values
are relative evidence summaries and should not be compared across runs as
calibrated probabilities.

`score_source_summary`, when present:

| Column or index | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| index | Kinase label. | string | Yes |
| `fused_motif_profile_evidence_count` | Sites with fused motif/profile evidence. | non-negative integer-like | Yes |
| `profile_only_motif_missing_or_constant_count` | Sites scored from profile because motif evidence was missing or constant. | non-negative integer-like | Yes |
| `profile_only_no_motif_overlap_count` | Sites scored from profile because no motif overlap was available. | non-negative integer-like | Yes |
| `unavailable_no_score_count` | Site/kinase pairs without an available score. | non-negative integer-like | Yes |
| `sites_with_score_count` | Site count with usable score. | non-negative integer-like | Yes |
| `total_sites_count` | Total site count considered. | non-negative integer-like | Yes |

`score_source_matrix`, when present, uses the same row and column axes as the
score matrix. Cell values are `fused_motif_profile_evidence`,
`profile_only_motif_missing_or_constant`, `profile_only_no_motif_overlap`, or
`unavailable_no_score`.

`profile_score_diagnostics`, when present:

| Column | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| `kinase` | Kinase label. | string | Yes |
| `substrate_site` | Dataset site identity. | `site_key` string | Yes |
| `status` | Scoring status. | `scored_after_leave_one_out` or `unscored` | Yes |
| `reason` | Explanation for unscored rows. | string; `insufficient_substrates_after_leave_one_out` for that failure | Yes |
| `substrates_before_leave_one_out` | Substrate count before self-exclusion. | non-negative integer-like | Yes |
| `substrates_after_leave_one_out` | Substrate count after self-exclusion. | non-negative integer-like | Yes |
| `min_substrates` | Configured minimum. | integer-like | Yes |

`KinasePredictionResult`:

| Attribute or helper | Python type | Always present? | Meaning |
| --- | --- | --- | --- |
| `pred_mat` | `pandas.DataFrame` | Yes | Site-by-kinase prediction support matrix. |
| `substrate_list` | `pandas.DataFrame | None` | Optional | Long-format predicted substrate list when produced. |
| `to_dataframe()` | `pandas.DataFrame` | Yes | Defensive snapshot of `pred_mat`. |
| `substrate_list_dataframe()` | `pandas.DataFrame | None` | Optional | Defensive snapshot of the substrate list. |

`pred_mat` rows are unique `site_key` values and columns are kinase labels.
Values are unit-interval support scores where present. Missing values indicate
that a score was unavailable; they do not prove biological absence.

Typical substrate-list columns include:

| Column | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| `kinase` | Kinase label. | string | When substrate list is present |
| `substrate_site` | Substrate identifier. | string; often `site_key` or projected site label depending on source | When substrate list is present |
| `score` | Prediction support score. | numeric | When substrate list is present |
| `rank` | Rank within the reported list. | integer-like | When substrate list is present |
| `site_key` | Dataset site identity. | string | Present in workflow-produced lists when available |
| `display_id` | Human-readable site label. | string | Present in workflow-produced lists when available |

`substrate_contributions`, when requested:

| Column | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| `kinase` | Kinase label. | string | Yes |
| `substrate_site` | Substrate site label from the scoring context. | string | Yes |
| `substrate_identifier` | Identifier used by the contribution record. | string | Yes |
| `value_used_in_scoring` | Numeric value consumed by scoring. | numeric or missing | Yes |
| `score_component` | Score component label. | string | Yes |
| `score_source` | Source of the contribution. | string | Yes |
| `reference_source_name` | Reference source name. | string or missing | Yes |
| `reference_source_version` | Reference source version. | string or missing | Yes |
| `reference_bundle_id` | Reference bundle identifier. | string or missing | Yes |
| `reference_identifier_namespace` | Reference identifier namespace. | string or missing | Yes |
| `status` | Inclusion status. | `included` or `excluded` | Yes |
| `exclusion_reason` | Reason for exclusion. | string or missing; missing for included rows | Yes |
| `ambiguous` | Whether the contribution is ambiguous. | boolean-like | Yes |

`KinaseActivityResult`, when present:

| Attribute | Python type | Always present? | Meaning |
| --- | --- | --- | --- |
| `activity_matrix` | `pandas.DataFrame` | Yes | Primary method-neutral kinase activity score matrix. |
| `p_value_matrix` | `pandas.DataFrame | None` | Method-dependent | P-values when the selected method and provenance support them. |
| `q_value_matrix` | `pandas.DataFrame | None` | Method-dependent | Adjusted p-values when produced. |
| `confidence_interval_low` | `pandas.DataFrame | None` | Method-dependent | Lower confidence bound when produced. |
| `confidence_interval_high` | `pandas.DataFrame | None` | Method-dependent | Upper confidence bound when produced. |
| `substrate_count_matrix` | `pandas.DataFrame` | Yes | Substrate counts used by the activity method. |
| `thresholded_substrate_mean_activity` | `pandas.DataFrame | None` | Method-dependent | Thresholded mean activity summary. |
| `thresholded_substrate_counts` | `pandas.DataFrame | None` | Method-dependent | Thresholded substrate counts. |
| `activity_substrate_counts` | `pandas.DataFrame | None` | Optional | Method-specific substrate counts. |
| `target_counts` | table-like or `None` | Optional | Target-count summary. |
| `target_table` | table-like or `None` | Optional | Target-level table. |
| `statistics_table` | `pandas.DataFrame | None` | Optional | Long-format statistics table when produced. |
| `count_field_semantics` | object or mapping | Optional | Count-field meaning. |
| `method_diagnostics` | object or mapping | Optional | Method diagnostics. |
| `policy_provenance` | object or mapping | Optional | Activity policy provenance. |
| `threshold_membership_diagnostics` | object or mapping | Optional | Threshold membership diagnostics. |
| `method_summary` | object or mapping | Optional | Method summary. |
| `input_semantics` | object or mapping | Optional | Activity input semantics. |
| `profile_metadata` | object or mapping | Optional | Profile metadata. |
| `membership_selection` | object or mapping | Optional | Substrate-membership selection evidence. |

`activity_result.activity_matrix` is the primary activity matrix for user code.
<!-- phospy-deprecation-compat: activities.result.activity_scores activities.result.weighted_activity -->
`activity_scores` and `weighted_activity` are compatibility aliases and are not preferred.

`KinaseEligibilityReport`, when present:

| Attribute | Meaning |
| --- | --- |
| `total_dataset_sites` | Dataset site count entering kinase validation. |
| `sequence_complete_sites` | Sites with complete sequence evidence. |
| `localisation_eligible_sites` | Sites passing localisation requirements. |
| `reference_overlap_sites` | Sites overlapping the reference projection. |
| `excluded_no_reference_match` | Sites excluded because no reference match was found. |
| `excluded_low_localisation` | Sites excluded by localisation policy. |
| `eligible_kinases` | Kinases meeting eligibility criteria. |
| `excluded_kinases_below_min_substrates` | Kinases excluded by substrate-count floor. |

Attrition summaries and provenance:

| Attribute | Meaning |
| --- | --- |
| `input_rows` | Rows entering the reported stage. |
| `rows_removed_during_preprocessing` | Rows removed during preprocessing. |
| `rows_removed_invalid_or_missing_site_identifiers` | Rows removed due to invalid or missing site identifiers. |
| `duplicate_sites_merged_or_resolved` | Duplicate-site rows merged or resolved by preprocessing. |
| `output_rows` | Rows leaving preprocessing. |
| `sequence_complete_sites` | Optional sequence-complete count. |
| `final_quantitative_sites_entering_scoring` | Quantitative sites entering scoring. |
| `sites_with_valid_site_sequence` | Sites with usable sequence. |
| `sites_without_usable_site_sequence` | Sites without usable sequence. |
| `sites_eligible_for_motif_scoring` | Sites eligible for motif scoring. |
| `sites_with_kinase_substrate_reference_profile_evidence` | Sites with reference profile evidence. |
| `sites_contributing_to_final_fused_prediction_scoring_output` | Sites contributing to final fused prediction/scoring output. |
| `sites_contributing_to_activity_scoring` | Sites contributing to activity scoring. |
| `attrition_provenance.to_payload()` | JSON-like metrics, policy, outcome (`passed`, `warned`, or `failed`), violations, and warnings. |

## Interpreting the result

Higher kinase scores indicate stronger support under the selected references,
sequence context, quantitative scale, and scoring mode. They are relative
support values, not calibrated causal inference.

`phosr_rank_weighted` is PhosR-inspired PhosPy scoring. It is not exact PhosR
numerical compatibility. Kinase Library-style modes require compatible
caller-supplied local resources and do not silently fall back when a required
resource is incompatible.

Optional activity output depends on substrate coverage, membership provenance,
and the selected method. KSEA-style p/q values are available only when typed
substrate-membership provenance shows that membership was selected independently
of the tested quantitative matrix. ssGSEA-style activity is not PTM-SEA support.
Causal kinase activity claims require study design support and external
validation.

Absence from a prediction, substrate, or contribution table should be read as
absence from the workflow output after filtering and eligibility checks, not as
biological absence.

## Common problems

| Problem | What to check |
| --- | --- |
| `ReferencePreset.AUTO` fails | Bundled references are rat-only. Use `Organism.RAT` or pass an explicit `ReferenceBundle`. |
| Missing sequence error | Ensure `site_metadata.site_sequence` is present, non-empty, centered, and compatible with the requested scoring mode. |
| Localisation failure | Add localisation metadata and configure `DatasetLocalisationConfig`; missing localisation fails when required. |
| Bare scoring config fails | Provide `KinaseScoringConfig(..., reliability_profile=...)` or a preset such as `KinaseScoringConfig.exploratory()`. |
| Reference display ambiguity | Keep the default `"error"` or use `"allow_with_diagnostics"` only when duplicate display labels represent intentional separate `site_key` rows. |
| Incompatible quantitative scale | Check the quantitative method contract table and dataset scale/meaning. |
| Empty predictions or sparse scores | Inspect eligibility, attrition, score-source diagnostics, reference overlap, and `min_substrates`. |
| Activity result is missing | `activity_config=None` or `enabled=False` disables activity; sparse substrates may also prevent method output. |

## Related documentation

- [Preparing a dataset](dataset-build-workflow.md)
- [Reference data](../reference_bundles.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)
- [Signalome analysis](signalome.md)
