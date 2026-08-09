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

For production workflow validation, use the scoring preset that carries the
same site-level threshold requirement into the kinase request:

```python
from phospy.advanced import KinaseScoringConfig

scoring_config = KinaseScoringConfig.production(
    minimum_reference_overlap_fraction=study_reference_overlap_floor,
    minimum_sequence_supported_fraction=study_sequence_support_floor,
    minimum_scored_fraction=study_scored_site_floor,
)
```

`KinaseScoringConfig.exploratory()` names the historical permissive profile
explicitly. Reference-context compatibility remains conservative by default:
unknown context fails unless
`ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT` is set
explicitly.

## Request Object

Use `KinaseWorkflowRequest`.

Important fields:

| Field | Meaning |
| --- | --- |
| `dataset` | The `AnalysisReadyPhosphoDataset` to score. |
| `references` | `ReferencePreset` or explicit `ReferenceBundle`. |
| `scoring_config` | Explicit `KinaseScoringConfig` for exploratory, production, or custom scoring intent plus scoring mode, substrate floors, diagnostics, localisation, and mixed total-protein guardrails. |
| `prediction_config` | `KinasePredictionConfig` for deterministic or adaptive prediction. |
| `activity_config` | `KinaseActivityConfig`, `None`, or disabled config for optional activity score output. The request default is `None`; activity execution is opt-in. |
| `site_sequence_conflict_policy` | Handles dataset/reference sequence conflicts: `"prefer_reference"`, `"prefer_dataset"`, or `"error"`. |
| `reference_display_ambiguity_policy` | Handles one-display-to-many-`site_key` reference projection: `"error"` or `"allow_with_diagnostics"`. |
| `kinase_library_resource` | Required only for Kinase Library-style workflow scoring modes. |

Constructing the request records intent only. `KinaseWorkflow.run(...)` validates
dataset, reference, config, localisation, sequence, and projection compatibility
before interpretation and execution.

## Request Configuration

Use these config objects. `scoring_config` is required at workflow validation:

- `KinaseScoringConfig`
- `KinasePredictionConfig`
- `KinaseActivityConfig`
- `LocalisationRequirement`

Choose scoring intent explicitly with `KinaseScoringConfig.exploratory()`,
`KinaseScoringConfig.production(...)`, or direct
`KinaseScoringConfig(..., reliability_profile=KinaseReliabilityProfile.CUSTOM)`
for custom values. Bare `KinaseScoringConfig()` is rejected.

Important `KinaseScoringConfig` fields:

| Field | Default | Notes |
| --- | --- | --- |
| `scoring_mode` | `"phosr_rank_weighted"` | Supported modes: `"phosr_rank_weighted"`, `"kinase_library_contextual_motif"`, `"kinase_library_motif_only"`, `"combined_profile_motif"`. The default is PhosR-inspired PhosPy scoring, not a PhosR compatibility mode. |
| `reliability_profile` | Required for direct construction | Use `exploratory()`, `production(...)`, or explicit `CUSTOM`; profiles are caller-selected, not inferred from old defaults. |
| `min_substrates` | `2` | Minimum unique usable substrates for kinase scoring support. The public floor is `2`. |
| `include_diagnostic_scoring_tables` | `False` | Adds non-primary diagnostic scoring tables. |
| `include_substrate_contributions` | `False` | Assembles and adds an optional substrate-level contribution table to `KinaseWorkflowResult`. |
| `profile_missing_value_strategy` | `"strict"` | Use `"median_skipna"` only when skipping missing profile values is intended. |
| `localisation_requirement` | `LocalisationRequirement()` | Workflow-level localisation requirement. Use `KinaseScoringConfig.production(...)` for the 0.75 production site-level threshold plus production attrition requirements. |
| `allow_mixed_total_protein_quantitative_meaning` | `False` | Keep `False` unless mixed corrected/uncorrected rows are intended. |

### Method Quantitative Input Contracts

Kinase scoring and activity methods declare their own quantitative input
contracts. The workflow validates the selected method contract before execution,
records the resolved contract in provenance, and does not transform linear,
log2, abundance, contrast, or effect inputs inside a scoring/activity method.

| Method | Accepted scale | Accepted meaning | Required centring/standardisation | Missing values | Profile axis | Statistical interpretation | P-value interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| kinase_scoring.combined_profile_motif | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| kinase_scoring.kinase_library_contextual_motif | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| kinase_scoring.kinase_library_motif_only | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio, unknown | Requires centered phosphosite sequence context; quantitative centring is not applicable because phospho values are not consumed by motif-only scoring.; No quantitative standardisation is required or performed for motif-only scoring. | Phospho missing values are not read by motif-only scoring; no missing-value transformation or imputation is performed. | Rows are protein-scoped site_key phosphosites with centered sequence context; quantitative columns are not used for motif-only score calculation. | Scores are sequence-motif support scores from the supplied Kinase Library-style resource, not abundance-profile statistics. | none |
| kinase_scoring.phosr_rank_weighted | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | Requires centered phosphosite sequence context; does not center quantitative values during scoring.; No automatic quantitative standardisation; values are consumed on the declared input scale. | Profile construction follows scoring_config.profile_missing_value_strategy; missing values are never imputed by the scoring method and no method-level imputation is performed. | Rows are protein-scoped site_key phosphosites; columns are aligned sample/profile abundance or total-corrected quantitative observations used for profile support. | Profile-derived relative support scores are within-run evidence summaries over the declared abundance/profile axis; linear and log2 inputs are scale-sensitive and are not numerically interchangeable. | none |
| ksea_zscore_v1 | log2 | phosphosite_log_abundance, phospho_total_log_ratio, contrast_log2_fold_change, differential_effect_size | Uses log2 sample, total-corrected, contrast, or effect profiles as declared by the dataset; no centring is performed in the method.; Requires log2 abundance, log2 total-corrected ratio, log2 contrast fold-change, or pre-standardised effect semantics; linear raw abundance is rejected. | Finite values define per-profile substrate and background sets; missing values are omitted from those calculations without imputation. | Columns must represent log-scale sample profiles, contrasts, or standardised effect profiles; linear raw samples are rejected. | Unweighted substrate-set z-score enrichment over declared log-scale sample, contrast, or effect values with background variance checks. | Two-sided normal-approximation p-values for computed z-scores; available only when typed substrate-membership provenance declares the membership independent of the tested quantitative matrix. Eligible p-values use Benjamini-Hochberg q-value adjustment per profile when enabled; adaptive membership reports descriptive z-scores with p/q unavailable. |
| simplified_weighted_substrate_activity_v1 | linear, log2 | phosphosite_abundance, phosphosite_log_abundance, phospho_total_log_ratio | No method-level centring; activity values are weighted means on the declared input scale.; No automatic standardisation; linear and log2 abundance summaries have different meanings. | Missing substrate values are ignored per profile when computing weighted and thresholded means; no imputation is performed. | Columns must represent sample-level abundance or explicit condition-summary abundance profiles. | Heuristic substrate-supported weighted mean; not a statistical enrichment test and not causal kinase activity proof. | none |
| ssgsea_substrate_enrichment_activity_v1 | log2 | contrast_log2_fold_change, differential_effect_size | Uses ranked contrast/effect values supplied by the caller; no centring is performed inside the method.; Requires log2 contrast fold-change or pre-standardised effect semantics; raw abundance is rejected. | Only finite effect values enter the ranked background; missing values are omitted without imputation. | Columns must represent contrasts or standardised effect profiles, not raw samples. | Rank-walk substrate-set enrichment summary over ordered effect values. Equal-valued finite sites are handled inside the method as tie blocks using the documented block-expectation policy, not row order or lexical site labels. Not PTM-SEA parity and not causal kinase activity proof. | No p-values are produced unless seeded permutations are requested; permutation p-values are two-sided empirical substrate-label permutation p-values, with Benjamini-Hochberg q-values per profile when enabled. |

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
| request `activity_config` | `None` | Activity execution is opt-in on `KinaseWorkflowRequest`. |
| `enabled` | `True` on `KinaseActivityConfig` | Set `activity_config=None` or `enabled=False` to skip activity. |
| `method` | `"simplified_weighted_substrate_activity"` | Supported methods: `"simplified_weighted_substrate_activity"`, `"ksea_zscore"`, `"ssgsea_substrate_enrichment"`. |
| `threshold` | `0.6` | Prediction support threshold. |
| `min_substrates` | `3` | Weighted activity-like score substrate floor. |
| `top_n_substrates` | `20` | Weighted activity substrate cap. |
| `ksea_min_substrates` | `5` | KSEA-style substrate floor. |
| `ssgsea_min_substrates` | `5` | ssGSEA-style substrate floor. |
| `ssgsea_random_seed` | `0` | Required when permutations are enabled; deterministic child streams are keyed by method, method version, `profile_id`, kinase, stream name, and caller-supplied seed under `stable_by_method_profile_kinase`. |

Activity score substrate support is counted from prediction/activity inputs,
not from the scoring profile count. For the simplified weighted and KSEA-style
methods, finite prediction support at or above the configured threshold is
included. KSEA-style and ssGSEA-style results keep not-computable
kinase-profile pairs in `statistics_table` with insufficient-substrate status
when the selected method can report that detail. Activity diagnostics also
expose `method_summary`, `substrate_count_matrix`,
`thresholded_substrate_counts`, and threshold-membership metadata where
supported.

KSEA-style p/q-value availability is gated by
`KinaseActivityResult.membership_selection`. Current payloads include
`membership_selection_schema_version="2"` and closed `selection_evidence`.
That evidence declares the scientific selection-process kind, score-source
kind, adaptive state, tested-matrix consumption state, source-specific contract
version, and typed independence evidence when applicable. The record also
retains descriptive provider method/version and score-source labels,
threshold/top-k parameters, reference fingerprints, selected kinase/substrate
universes, and two separate matrix fingerprints:
`selection_quantitative_matrix_fingerprint` for the quantitative matrix used
during membership selection when applicable, and
`tested_quantitative_matrix_fingerprint` for the exact KSEA background matrix.
The activity science domain derives `source_category` and
`inferential_decision`; workflow code does not supply final eligibility.
Serialized source category, status, reason, and missing-evidence fields are
cross-validated against the reconstructed typed evidence. Contradictory records
are rejected.

| Membership source category | Ordinary KSEA p/q behavior |
| --- | --- |
| `profile_derived`, including leave-one-out profile scoring | Descriptive z-scores only. |
| `fused_profile_motif` | Descriptive z-scores only. |
| data-adaptive `prediction_selected` | Descriptive z-scores only. |
| `unknown` or incomplete provenance | Descriptive z-scores only. |
| `fixed_external_reference` | p/q values only when fixed-reference independence evidence, source fingerprints, tested-matrix fingerprint, and selected universes are present. |
| `sequence_only_motif` | p/q values only when sequence-only independence evidence, source fingerprints, tested-matrix fingerprint, selected universes, and the motif score source are present. |

Provider-specific method names, provider versions, score-source labels, and
reference descriptions are descriptive only. Arbitrary strings cannot establish
fixed, external, sequence-only, or independent membership. Fixed-reference and
motif-only inference require closed source-specific selection evidence,
explicit `data_adaptive_membership=False`, `consumed_tested_matrix=False`, no
selection quantitative-matrix fingerprint, supported contract versions,
matching typed independence-policy evidence, the exact tested-matrix
fingerprint, source-reference fingerprints, and non-empty selected universes.
`threshold_top_k_policy` may retain numerical configuration such as top-k and
thresholds, but current payloads must not rely on generic policy-map keys for
adaptive state, independence policy, or source classification.

For ineligible membership, KSEA reports `p_value_matrix=None`,
`q_value_matrix=None`, missing `p_value`/`q_value` cells in
`statistics_table`, and explicit inferential status/reason fields. Direct use
of `KseaZScoreActivityMethod` applies the same science-domain policy and checks
tested-matrix and selected-universe provenance before allocating p/q outputs.
Every KSEA result carries explicit membership provenance. Missing or legacy
membership provenance is represented as an explicit missing/ineligible record
for descriptive-only output; legacy open-string evidence fails closed and is
not migrated to ordinary inference unless it matches a reviewed closed rule.
Finite p/q output without eligible membership
provenance is rejected. Legacy membership payloads are loaded only if serialized
eligibility, status, reason, and nested decision fields agree with the decision
recomputed from the preserved facts. Serialized source-category relabelling,
adaptive facts under a fixed-external label, and added favourable independence
tokens on adaptive records are rejected. Missing fingerprints are not fabricated
or upgraded to eligibility.

Ordinary KSEA normal-approximation assumptions remain scientific assumptions
even after eligibility is established. Eligibility means the ordinary p/q
output is permitted by the provenance gate, not that causal kinase activity has
been proven.

`KinaseActivityResult` also carries typed activity profile semantics through
`input_semantics` and `profile_metadata`. These objects define the profile axis,
quantitative semantics, profile identifiers, sample/condition/contrast
identifiers, and condition-summary aggregation metadata when condition-summary
profiles are supplied. Method identity does not determine these semantics:
for example, KSEA can consume sample log-abundance, contrast log-fold-change,
or standardised-effect profiles according to the explicit input semantics.
For activity statistics, `profile_id` is the required row identifier:
sample-axis results use sample IDs, condition-summary results use condition
summary IDs, contrast-axis results use contrast IDs, and effect-axis results
use neutral effect profile IDs. A `condition` column is reserved for genuine
condition-summary results, where it must match `profile_id`; sample, contrast,
and effect statistics tables do not contain `condition`. Use
`legacy_condition_statistics_table_dataframe()` only for deprecated consumers <!-- phospy-deprecation-compat: activities.result.legacy_condition_statistics_table -->
that still need an old condition-shaped table; the adapter adds
`condition = profile_id` and does not create a biological condition contract.

## Running the Workflow

```python
from phospy import KinaseWorkflow
from phospy.advanced import KinaseScoringConfig
from phospy.api import KinaseWorkflowRequest, ReferencePreset

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig.exploratory(),
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
| `attrition_provenance` | Optional immutable attrition metrics, policy, outcome, violations, and warning messages. |
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
`activity_scores` and `weighted_activity` are not preferred for new <!-- phospy-deprecation-compat: activities.result.activity_scores, activities.result.weighted_activity -->
documentation or code.

Reloadable kinase bundles preserve the exact typed activity semantics in
manifest version 3. Loading a current bundle reconstructs `input_semantics` and
`profile_metadata` from their persisted typed payloads and rejects contradictions
between those payloads, the activity matrix columns, condition-summary
aggregation records, or resolved activity semantics recorded in provenance.
KSEA membership-selection provenance and inferential eligibility are also
persisted; loading preserves whether ordinary p/q values were available or
unavailable.
Kinase version-2 bundles did not persist enough semantic metadata for faithful
activity reconstruction and must be regenerated with a current PhosPy version.

### Substrate Contribution Table

Set a custom `KinaseScoringConfig` with `include_substrate_contributions=True`
to assemble and
attach `kinase_result.substrate_contributions`.

```python
from phospy.advanced import KinaseReliabilityProfile, KinaseScoringConfig

request = KinaseWorkflowRequest(
    dataset=dataset,
    references=references,
    scoring_config=KinaseScoringConfig(
        reliability_profile=KinaseReliabilityProfile.CUSTOM,
        min_substrates=2,
        include_substrate_contributions=True,
    ),
    activity_config=None,
)

kinase_result = KinaseWorkflow().run(request)
contributions = kinase_result.substrate_contributions
```

The table contains one row per projected kinase-substrate evidence row. It is
assembled only when explicitly requested and is off by default because it can be
large.

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

Kinase provenance separates causal site-row attrition from compatibility
metrics. `row_attrition`, when present, records only stage-local site-row
removal proven from that stage's input and output site indexes. It does not
encode site/kinase-pair loss as row loss. `row_attrition_metrics` remains
available for legacy compatibility diagnostics, including
`site_kinase_pairs_considered`, `site_kinase_pairs_scored`, and unscored-pair
counts.

`workflow_parameters["reference_projection_summary"]` records the source
kinase-substrate reference projection before unmatched substrate identifiers are
dropped. `workflow_parameters["universe_attrition"]` always contains separate
categories for `reference_attrition`, `sequence_attrition`,
`membership_attrition`, `finite_value_attrition`, and
`activity_background_attrition`. Reference attrition is expressed in the source
`references.kinase_substrate_map.substrate_site` namespace and includes bounded
examples of unmatched identifiers; projected scoring and activity universes
remain in dataset `site_key` identity.

The reference-projection summary is strictly validated schema version 1 provenance.
PhosPy owns and validates its source/output namespaces, identity-semantics text,
identifier-kind vocabulary, projector version, count invariants, one-to-many
diagnostic token, and the bounded unmatched-example policy. Bundle
reconstruction rejects unsupported or contradictory projection-summary payloads
and checks that `universe_attrition["reference_attrition"]` agrees with the
typed summary instead of treating the two fields as independent authority.

`KinaseWorkflowAttritionProvenance` stores the workflow-calculated attrition
metrics, configured policy payload, policy outcome, policy violations, and
warning messages. The result contract does not calculate those values. It
recursively freezes the JSON-like `metrics`, `policy`, and violation details at
construction, rejects invalid JSON state instead of stringifying it, and
returns fresh ordinary `dict`/`list` payloads from `to_payload()`.

## Limitations

- Bundled runtime references are rat-only for `ReferencePreset.AUTO` in this
  release; human or mouse analysis requires an explicit caller-supplied
  `ReferenceBundle`.
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
from phospy.advanced import (
    KinaseReliabilityProfile,
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api import (
    KinaseWorkflowRequest,
    ReferencePreset,
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
)

kinase_result = KinaseWorkflow().run(request)

scores = kinase_result.scoring_result.authoritative_scores
predictions = kinase_result.prediction_result.pred_mat
print(scores.shape)
print(predictions.head())
```
