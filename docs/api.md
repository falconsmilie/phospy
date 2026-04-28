# API Guide

This guide describes the supported public Python API.

If you are new to PhosPy, start with the [Quickstart](getting-started/quickstart-first-workflow.md) first.

PhosPy does not expose HTTP endpoints. Its supported public interfaces are the
Python API documented here and the file-first [`phospy` CLI](cli.md).

## Import contract

`phospy.api` is the primary namespace where public API types are defined and organised in source.

Use the package namespaces like this:

- top-level `phospy` is a curated convenience surface for only:
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
  `KinaseWorkflow`, `SignalomeWorkflow`
- use `phospy.api` for requests, configs, results, enums, references, and
  public exceptions

Examples:

```python
from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
```

```python
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
```

## Public workflow shape

PhosPy has one dataset boundary and two workflow entrypoints:

1. `DatasetBuildRequest -> AnalysisReadyDatasetBuilder.run(...) -> AnalysisReadyPhosphoDataset`
2. `KinaseWorkflowRequest -> KinaseWorkflow.run(...) -> KinaseWorkflowResult`
3. `SignalomeWorkflowRequest -> SignalomeWorkflow.run(...) -> SignalomeWorkflowResult`

All public executors use `run(request)`.

## Beginner lane

The easiest supported lane is:

1. build a dataset with `organism=Organism.RAT`
2. run kinase with `references=ReferencePreset.AUTO`
3. run signalome only when `site_metadata.protein_id` is present

Bundled runtime references are rat-only in this release. Human and mouse work
need an explicit `ReferenceBundle`.

## Main public entrypoints

| Type | Purpose |
| --- | --- |
| `AnalysisReadyDatasetBuilder` | Builds a validated workflow-ready dataset |
| `AnalysisReadyPhosphoDataset` | Strict workflow input model |
| `KinaseWorkflow` | Runs kinase scoring, prediction, and optional activity |
| `SignalomeWorkflow` | Runs signalome analysis from a kinase result |

## Request models

### `DatasetBuildRequest`

Fields:

- `phospho`: DataFrame or supported file path
- `site_metadata`: DataFrame or supported file path
- `sample_metadata`: optional DataFrame or path
- `total`: optional DataFrame or path
- `organism`: optional `Organism`
- `preprocessing_config`: `DatasetPreprocessingConfig`

Builder notes:

- supported inputs are pandas `DataFrame` values or file paths
- supported site-metadata aliases are narrow: `gene_name` may stand in for `gene_symbol`, and `centralized_sequence` may stand in for `site_sequence`
- unsupported historical aliases such as `gene`, `residue`, `phosphosite`, `site_position`, `sequence`, and `protein` are rejected
- the supported public builder lane must end in a missing-value-free `AnalysisReadyPhosphoDataset`

### `KinaseWorkflowRequest`

Fields:

- `dataset`: `AnalysisReadyPhosphoDataset`
- `references`: `ReferencePreset` or `ReferenceBundle`
- `scoring_config`: `KinaseScoringConfig`
- `prediction_config`: `KinasePredictionConfig`
- `activity_config`: `KinaseActivityConfig | None`

### `SignalomeWorkflowRequest`

Fields:

- `kinase_result`: `KinaseWorkflowResult`
- `config`: `SignalomeConfig`

Signalome requires explicit, non-empty `dataset.site_metadata.protein_id` for
all interpreted sites.

## Config models

Public config dataclasses are strict constructors: invalid local policy or
numeric state is rejected when the config object is created.

Request/workflow validators still run and still own cross-object checks (for
example, `subtract_log_total` requiring `total` and
`intensity_transform.policy="log2"`, `sample_metadata_pairs` requiring
`sample_metadata`, reference compatibility, and workflow input object
boundaries).

### Dataset preprocessing

`DatasetPreprocessingConfig` groups six builder-owned areas:

- `intensity_transform: DatasetIntensityTransformConfig`
- `normalisation: DatasetNormalisationConfig`
- `missing_data: DatasetMissingDataConfig`
- `total_protein_correction: DatasetTotalProteinCorrectionConfig`
- `site_matrix: DatasetSiteMatrixConfig`
- `comparisons: DatasetComparisonBuildingConfig`

Performance contracts for preprocessing and scoring lanes are documented in
[Performance Contracts](performance.md).

All preprocessing methods are explicit opt-in choices. Defaults remain
conservative (`intensity_transform="identity"`, `normalisation="none"`,
`missing_data="forbid"`).

#### `DatasetIntensityTransformConfig`

- `policy="identity"` (default)
- `policy="log2"` applies `log2(value + pseudocount)` to quantitative intensities
- `pseudocount` must be non-negative

Use `log2` when intensities are positive (or become positive after adding the
configured pseudocount).

#### `DatasetNormalisationConfig`

- `policy="none"` (default)
- `policy="median_center"` subtracts each sample column median
- `policy="quantile"` forces sample columns to share one empirical distribution

Use quantile normalisation only when matched-distribution assumptions are
scientifically appropriate for your experiment design. It is a dense sort-heavy
operation with additional float64 matrix-copy cost; see
[Performance Contracts](performance.md#quantile-normalisation).

#### `DatasetMissingDataConfig`

- `policy="forbid"` keeps the strict default lane
- `policy="impute_row_median"` drops rows below `min_observed_values`, then imputes remaining missing values with the row median

#### `DatasetTotalProteinCorrectionConfig`

- `policy="none"`
- `policy="subtract_log_total"` (recommended)
- `identity: DatasetTotalProteinCorrectionIdentityConfig`

Subtractive total-protein correction is log-scale correction:

`corrected = log2(phospho + pseudocount) - log2(total + pseudocount)`

In the public builder lane this requires:

- `intensity_transform.policy="log2"`
- a `total` table aligned to phospho sample columns
- an explicit identity mapping policy for phosphosite-to-total matching

`DatasetTotalProteinCorrectionIdentityConfig` makes identity mapping explicit:

- `mode="direct"`:
  - map `site_metadata[phosphosite_key]` directly to total identity
  - total identity is resolved from `total.index` (use `total_protein_key="__index__"` for index-based matching)
- `mode="mapping_table"`:
  - provide `mapping_table` plus `mapping_phosphosite_key` and `mapping_total_protein_key`
  - use when phosphosite and total identifiers are in different namespaces

Supported strictness controls:

- `duplicate_policy="error"` (default)
- `unmatched_policy="error"` (default) or `unmatched_policy="allow_uncorrected"`

Recommended identity preference order:

1. explicit `mapping_table` mode
2. direct accession/protein-group identifiers
3. direct gene-symbol matching only when scientifically appropriate

Gene-symbol matching is a convenience identity policy, not a universal
biological identity guarantee. Isoform-specific, protein-group, or shared-peptide
datasets should prefer accession/protein-group IDs or explicit mapping tables.

By default PhosPy fails loudly on ambiguous identity states:

- duplicate total identity keys
- duplicate or ambiguous mapping-table rows
- null/empty identity keys
- unknown mapping-table references
- unmatched phosphosite rows (unless `unmatched_policy="allow_uncorrected"`)

When correction runs, diagnostics/provenance records include:

- identity mode and keys
- duplicate/unmatched policies
- mapping-table fingerprint (mapping-table mode)
- corrected/uncorrected/unmatched counts
- unused total-protein rows
- whether gene-symbol matching was used and warning text when applicable

Important state semantics:

- numeric scale and scientific meaning are tracked separately at dataset
  boundary
- subtractive correction keeps numeric scale at `log2`, but changes quantitative
  meaning to `phospho_total_log_ratio`

#### `DatasetSiteMatrixConfig`

- `policy="as_input"`
- `policy="build_from_metadata"`
- `duplicate_site_policy`: `max_mean_signal`, `first`, `aggregate_mean`, `aggregate_median`, `error`
- public missing-data handling for this stage: `missing_data_policy="drop_any_missing"`

This public builder lane is intentionally strict and still ends in a missing-value-free `AnalysisReadyPhosphoDataset`.
`minimum_observed_values` is internal-only state and must stay
`None` in the public config lane.

Duplicate-site policy trade-offs:

`duplicate_site_policy` controls a scientific row-resolution choice. Some
policies drop peptide context; aggregate policies preserve all rows numerically
but collapse distinct peptide/site contexts into one site-level row.

- `error`: cautious mode; fail on duplicate constructed sites instead of silently choosing one row.
- `first`: simple and convenient, but later duplicate rows are discarded by input order.
- `max_mean_signal`: keeps the strongest observed row, but can bias toward high-abundance / strong-signal rows.
- `aggregate_mean` / `aggregate_median`: preserve duplicate rows numerically via aggregation, but can blur distinct phosphosite or peptide contexts.

When site-matrix duplicate handling runs, `dataset.preprocessing_report` includes:

- `row_audit`: unified row-level preprocessing audit trail across stages (`missing_data`, `site_matrix`, and future row-excluding stages).
- `duplicate_site_resolution`: one row per source duplicate row, including retained/dropped or aggregated contribution details.
- `metadata_conflicts`: duplicate-site metadata disagreement records (for example conflicting `protein_id` or `site_sequence` values).

#### `DatasetComparisonBuildingConfig`

- `policy="none"`
- `policy="sample_metadata_pairs"`
- `sample_group_column` defaults to `comparison_group`
- `pairs` may be provided explicitly, otherwise observed group pairs are inferred

When comparison building runs:

- `dataset.comparisons` remains the compact site-by-comparison effect-size matrix used by workflows.
- `dataset.preprocessing_report.comparison_group_stats` provides replicate-group summary context (for example `n`, `mean`, `sd`, `sem`) for each site and group.
- `dataset.preprocessing_report.comparison_pair_stats` provides pairwise evidence (left/right group summaries plus `effect_size`) for each site/comparison row.

These sidecar tables improve transparency and auditability. They are not a replacement for full differential phosphoproteomics modelling.

Planned/future lanes (not currently supported in this public contract) include:
`knn` imputation, `min_prob` imputation, and `combat` batch correction.

### Kinase configs

#### `KinaseScoringConfig`

- `min_substrates` default `2`
- `include_diagnostic_scoring_tables` default `False`
- `profile_missing_value_strategy`: `strict` or `median_skipna`

Motif sequence-context validation is strict in the supported lane:

- default scoring semantics expect a centred phosphosite window of length `15`
  (flank size `7`, centre index `7`)
- default scoring does not midpoint-crop longer sequences
- longer centred sequences are only accepted when the caller explicitly opts into
  centred-sequence extraction semantics (`sequence_semantics="centred_sequence"`)
- supported residue alphabet is the 20 canonical amino acids:
  `A C D E F G H I K L M N P Q R S T V W Y`
- scored sequence identifiers must be site-shaped (`<protein>;<residue><position>;`)
  so centre residue compatibility can be validated
- query/target sequences that are missing, short, non-centred, site-residue mismatched,
  non-phospho-centre (`S/T/Y`), or unsupported are excluded from motif scoring
- reference/library sequences used to build motif profiles are validated with the
  same sequence rules before profile construction; invalid reference windows are
  excluded from motif frequency/profile construction (never neutral/partial
  encoded)
- explicit motif-sequence libraries support:
  - bare sequence entries (less-informative), and
  - structured entries with metadata equivalent to
    `reference_id`, `site_id`, `kinase`, and `sequence`
- bare explicit entries validate motif-window quality, supported residue alphabet,
  and phospho-compatible centre residue (`S/T/Y`), but cannot always prove that
  the intended phosphosite identity matches the centre residue
- structured explicit entries are preferred for reproducible motif-library
  construction because supplied `site_id` can be format-validated and checked for
  centre-residue agreement
- site-residue mismatch checks require site-like metadata, such as `site_id` or a
  reference identifier that encodes the phosphosite
- excluded query sites are reported through
  `result.scoring_result.motif_sequence_validation`
- excluded/accepted reference windows are reported through
  `result.scoring_result.motif_library_validation`

Performance notes:

- motif scoring scales with scored sites, eligible kinases, and motif window width,
- the scoring lane filters motif work to eligible profile-overlap kinases,
- enabling diagnostic scoring tables retains `motif_scores` and
  `score_fusion_weights`, which can increase runtime/memory.

See [Performance Contracts](performance.md#kinase-scoring-and-motif-scoring).

#### `KinasePredictionConfig`

- `top_k` default `30`
- `deterministic_max_selected_kinases` default `10`
- `adaptive_ensemble_runs` default `10`
- `mode`: `deterministic_ranking` or `adaptive_ensemble`
- `adaptive_policy`: `stable` or `r_parity`
- `n_iterations` default `5`
- `random_state` optional

`deterministic_max_selected_kinases` controls how many kinases are selected in
deterministic prediction.
`adaptive_ensemble_runs` controls how many ensemble runs are executed in
adaptive prediction.

#### `KinaseActivityConfig`

- `enabled` default `True`
- `threshold` default `0.6`
- `min_substrates` default `3`
- `top_n_substrates` default `20`

Set `activity_config=None` to skip the activity stage.

### `SignalomeConfig`

Fields:

- `substrate_support_cutoff`
- `network_correlation_threshold`
- `network_policy`: `positive_only`, `absolute_threshold`, `signed`
- `assignment_policy`: `cutoff_binary`, `weighted_top`
- `score_preconditioning_policy`: `allow_and_report`, `error_on_drop`
- `module_count`
- `module_selection_primary_correlation_threshold`
- `module_selection_fallback_correlation_threshold`
- `module_selection_max_clusters`
- `cluster_tree_backend`: `exact`
- `candidate_scoring_backend`: `full`, `sampled`
- `max_exact_cluster_tree_sites` (default `2000`)
- `max_full_correlation_sites` (default `2000`)

Signalome stages are now explicit:

- cluster-tree construction: controlled by `cluster_tree_backend`
- candidate scoring for module-count selection: controlled by `candidate_scoring_backend`
- kinase correlation calculation for network edges: controlled by `network_policy` + `network_correlation_threshold`
- module assignment: controlled by `assignment_policy` + `substrate_support_cutoff`
- network generation: built from downstream score correlations after module selection

Important: `candidate_scoring_backend="sampled"` does not imply approximate
signalome clustering. It only changes candidate module-count evaluation. Exact
cluster-tree construction and final module assignment still require the exact
tree and remain hard-guarded by `max_exact_cluster_tree_sites`.
Exact tree construction is an internal, scale-limited implementation detail of
the signalome workflow contract.
Low-level clustering helpers in `phospy.signalomes.clustering` use the same
guarded exact-tree path and resolve missing/`None`
`max_exact_cluster_tree_sites` to the safe default (`2000`), rather than an
unbounded tree build.

Successful runs record scale decisions in provenance under
`result.provenance.workflow_parameters["scale_guard"]` (`site_count`,
`cluster_tree_backend`, `candidate_scoring_backend`,
`candidate_scoring_requested_backend`,
`max_exact_cluster_tree_sites`, `max_full_correlation_sites`,
`exact_cluster_tree_built`, `candidate_scoring_mode`,
`candidate_scoring_evaluated`, `candidate_scoring_skip_reason`,
`candidate_scoring_sampling`, `candidate_scoring_applies_to`,
`final_module_assignment_backend`,
`final_module_assignment_uses_candidate_scoring`, `scale_guard_passed`).

`candidate_scoring_applies_to` is always
`"candidate_module_count_evaluation_only"`, and
`final_module_assignment_uses_candidate_scoring` is always `False`.
`candidate_scoring_requested_backend` records the backend requested in config.
`candidate_scoring_mode` records what was actually evaluated (`"full"`,
`"sampled"`, or `"not_evaluated"`).
When `module_count` is set explicitly, `candidate_scoring_evaluated` is `False`,
`candidate_scoring_skip_reason` is `"explicit_module_count"`, and
`candidate_scoring_sampling` is `None`.
These fields make stage ownership explicit:

- exact cluster-tree construction: `cluster_tree_backend`,
  `max_exact_cluster_tree_sites`, `exact_cluster_tree_built`
- candidate module-count evaluation: `candidate_scoring_*`
- final module assignment: `final_module_assignment_*`

Scientific score semantics are recorded separately under
`result.provenance.workflow_parameters["signalome_score_semantics"]`. This
payload captures:

- upstream downstream score source (`rank_weighted_fusion_scores` when present,
  otherwise `profile_scores`)
- meaning of the selected downstream score source
- module-selection score meaning (within-cluster correlation summaries over
  downstream score profiles)
- candidate scoring mode and scope
- network-correlation meaning and negative-correlation handling by
  `network_policy`
- missing/constant profile handling
- thresholds and limits used for interpretation
- clustering backend name
- scientific interpretation limits (derived summary semantics; not
  probabilities/calibrated confidence/causal proof)

When sampled candidate scoring actually runs,
`candidate_scoring_sampling` records reproducible sampling details for
candidate module scoring:

- `sampling_cap`
- `sampling_method`
- `deterministic_seed_policy`
- `actual_sampled_pair_count`
- `per_cluster_sample_count_summary` (`min`, `max`, `mean`, `total`)

When candidate scoring is skipped (for example explicit `module_count`),
sampled and full-correlation candidate-scoring diagnostics are not emitted.

See [Performance Contracts](performance.md#signalome-clustering-and-module-selection).

## Result models

### `AnalysisReadyPhosphoDataset`

This is the strict workflow-facing dataset boundary.

Key fields include:

- `phospho`
- `site_metadata`
- `sample_metadata` (optional)
- `total` (optional)
- `comparisons` (optional)
- `organism` (optional)
- `intensity_scale_state`
- `processing_state`
- `preprocessing_report` (optional)
- `provenance` (optional `RunProvenance`; machine-readable audit/replay metadata)

State model responsibilities:

- `intensity_scale_state` answers: "Are quantitative values linear or log2?"
- `intensity_scale_state.quantity` answers: "What do phospho matrix values
  represent scientifically?"
- `processing_state` answers: "What preprocessing policy state crossed the analysis-ready boundary?"
- `preprocessing_report` answers: "What happened during preprocessing, including row-level operations and sidecars?"

`intensity_scale_state.label` (for example `log2`) is never a complete
scientific interpretation on its own; always read it together with
`intensity_scale_state.quantity`.

Common combinations:

- linear phosphosite abundance:
  `label="linear"`, `quantity="phosphosite_abundance"`
- log2 phosphosite abundance:
  `label="log2"`, `quantity="phosphosite_log_abundance"`
- log2 phospho/total corrected ratio:
  `label="log2"`, `quantity="phospho_total_log_ratio"`

`preprocessing_report` and `provenance` serve different purposes:

- `preprocessing_report` is human-facing and table-oriented (`row_counts`,
  `operations`, `row_audit`, duplicate-site and comparison sidecars).
- `provenance` is machine-readable and contract-oriented (table fingerprints,
  environment versions, preprocessing execution hashes, workflow parameters,
  random-state metadata, and output fingerprints).

Quick row-audit inspection example:

```python
report = dataset.preprocessing_report
dropped = report.row_audit[report.row_audit["action"] == "dropped"]
print(dropped[["stage", "source_row_id", "site_id", "reason"]])
```

### `KinaseWorkflowResult`

Key fields:

- `dataset`
- `references`
- `scoring_result`
- `prediction_result`
- `activity_result` (optional)
- `provenance` (optional `RunProvenance`)

Common nested outputs:

- `result.scoring_result.profile_scores`
- `result.scoring_result.rank_weighted_fusion_scores`
- `result.scoring_result.motif_sequence_validation`
- `result.scoring_result.motif_library_validation`
- `result.prediction_result.pred_mat`
- `result.activity_result.weighted_activity` when activity is enabled
- `result.activity_result.thresholded_substrate_mean_activity` when activity is enabled

Method notes:

- `thresholded_substrate_mean_activity` is a simple mean phospho signal over
  predicted substrates whose prediction score is greater than
  `activity_config.threshold`. It is not full KSEA enrichment.
- `rank_weighted_fusion_scores` combine profile-correlation scores and
  motif-frequency scores using rank-derived weights from profile substrate
  counts and motif library sizes.
- `motif_sequence_validation.summary()` returns workflow diagnostics with:
  `total_sequences`, `valid_sequences`, `invalid_sequences`,
  `short_sequences`, `off_centre_sequences`, `site_residue_mismatches`,
  `unsupported_residue_characters`, and
  `sequences_excluded_from_motif_scoring`.
- Excluded sequence rows keep traceability through
  `motif_sequence_validation.rows`; excluded sites are not assigned ordinary
  motif scores (`motif_scores` remains missing/`NaN` for those rows).
- `motif_library_validation.summary()` reports motif-library build diagnostics,
  including provided/accepted/excluded reference counts, exclusion reasons
  (`missing`, `short`, `unsupported`, `off_centre`, `site_residue_mismatch`,
  `invalid_site_id`, `non_phospho_centre_residue`),
  and accepted-window/unsupported-residue policies.
- `motif_library_validation.rows` preserves per-reference provenance
  (`reference_id`, `site_id`, `kinase`, `sequence`, `status`, `reason`,
  observed/expected centre residues).

### `SignalomeWorkflowResult`

Key fields:

- `dataset`
- `kinase_result`
- `module_assignments`
- `signalome_modules`
- `kinase_network`
- `module_selection_diagnostics`
- `score_preconditioning_diagnostics`
- `expanded_signalome` (optional by contract, populated in the supported lane)
- `site_membership` (optional; site-level signalome membership provenance)
- `protein_site_context` (optional; protein-level multi-site context summary)
- `provenance` (optional `RunProvenance`)

Common nested outputs:

- `result.module_assignments.table`
- `result.signalome_modules.table`
- `result.kinase_network.edges`
- `result.kinase_network.candidate_correlations`
- `result.kinase_network.correlation_diagnostics`
- `result.expanded_signalome`
- `result.site_membership`
- `result.protein_site_context`

Signalome correlation semantics:

- PhosPy preserves undefined kinase correlations as missing (`NaN`) values.
- `0.0` means correlation was estimated and is near zero.
- Missing correlation means estimation was not possible (for example constant
  profiles, insufficient paired observations, missing values, or non-finite
  inputs).
- `result.kinase_network.candidate_correlations` includes:
  `source_kinase`, `target_kinase`, `correlation`, `correlation_status`,
  `valid_observations`, and `correlation_reason`.
- By default, only rows with `correlation_status == "finite"` are eligible for
  edge creation.
- `result.kinase_network.correlation_diagnostics` reports finite/undefined
  counts and how many candidates were skipped for non-finite correlation.
- Signalome module/network scores are derived summary statistics over upstream
  downstream score profiles; they are not probabilities, calibrated confidence
  values, or direct evidence of causal regulation.
- Use `result.provenance.workflow_parameters["signalome_score_semantics"]` for
  the run-specific interpretation payload (score source, score meaning,
  candidate-scoring scope, missing/constant handling, threshold/limit context,
  and interpretation limits).

Signalome context note:

- `result.signalome_modules.table` remains the compact protein/module summary table.
- `result.site_membership` shows site-level membership context, including excluded rows and reasons.
- `result.protein_site_context` highlights multi-site proteins and flags potentially ambiguous biological interpretation when site clusters or module context disagree.

## Enums and references

### `Organism`

Public organism enum values:

- `human`
- `mouse`
- `rat`

### `ReferencePreset`

Public preset enum values:

- `auto`
- `human`
- `mouse`
- `rat`

Important: enum presence does not guarantee bundled runtime data for every
organism in this release.

### `ReferenceBundle`

Use `ReferenceBundle` when you want to provide your own references.

Fields:

- `organism`
- `kinase_substrate_map`
- `site_sequences`

Large reference maps are supported; workflow runtime depends more on
dataset/reference overlap after filtering than on raw map size. See
[Performance Contracts](performance.md#large-kinase-substrate-references).

## Public exceptions

All user-facing exception types are available from `phospy.api`.

Common ones to catch:

- `PhosPyValidationError`
- `UnsupportedInputFormatError`
- `ReferenceCompatibilityError`
- `ReferenceResolutionError`
- `UnsupportedOrganismError`
- `WorkflowValidationError`
- `WorkflowBoundaryError`
- `SignalomeScaleError`

## Small working example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import DatasetBuildRequest, KinaseWorkflowRequest, Organism, ReferencePreset

phospho = pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"])
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14"],
        "site": ["Y182"],
        "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
        "protein_id": ["MAPK14"],
    },
    index=phospho.index,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
)

result = KinaseWorkflow().run(
    KinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)
```

## Where next

- [Quickstart](getting-started/quickstart-first-workflow.md)
- [Validation Guide](validation.md)
- [CLI Guide](cli.md)
- [Output Bundles](output_bundles.md)
