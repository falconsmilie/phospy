# API Guide

This guide describes the supported public Python API.

If you are new to PhosPy, start with the [Quickstart](getting-started/quickstart-first-workflow.md) first.

PhosPy does not expose HTTP endpoints. Its supported public interfaces are the
Python API documented here and the file-first [`phospy` CLI](cli.md).

## Import contract

`phospy.api` is the canonical namespace where public API types are defined and organised in source.

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
- unsupported legacy aliases such as `gene`, `residue`, `phosphosite`, `site_position`, `sequence`, and `protein` are rejected
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
- `policy="ratio_to_total"` (deprecated alias that resolves to `subtract_log_total`)

Subtractive total-protein correction is log-scale correction:

`corrected = log2(phospho + pseudocount) - log2(total + pseudocount)`

In the public builder lane this requires:

- `intensity_transform.policy="log2"`
- a `total` table aligned to phospho sample columns and site-to-protein mapping

#### `DatasetSiteMatrixConfig`

- `policy="as_input"`
- `policy="build_from_metadata"`
- `duplicate_site_policy`: `max_mean_signal`, `first`, `aggregate_mean`, `aggregate_median`, `error`
- public missing-data handling for this stage: `missing_data_policy="drop_any_missing"`

This public builder lane is intentionally strict and still ends in a missing-value-free `AnalysisReadyPhosphoDataset`.
`minimum_observed_values` is internal-only compatibility state and must stay
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

- scoring expects a centered motif window of length `15` (flank size `7`, centre index `7`)
- supported residue alphabet is the 20 canonical amino acids:
  `A C D E F G H I K L M N P Q R S T V W Y`
- sequences that are missing, short, off-centre, site-residue mismatched, or
  contain unsupported residue characters are excluded from motif scoring
- excluded sites are reported through
  `result.scoring_result.motif_sequence_validation`

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
Legacy `ensemble_size` is accepted as a deprecated constructor alias.

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

Signalome module-selection scoring uses a full site-by-site correlation path up
to `MAX_FULL_CORRELATION_SITE_COUNT` (`2000`), then switches to sampled
within-cluster correlation estimates above that threshold. Approximation use is
reported through `result.module_selection_diagnostics.reason`.

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
- `processing_state` answers: "What preprocessing policy state crossed the analysis-ready boundary?"
- `preprocessing_report` answers: "What happened during preprocessing, including row-level operations and sidecars?"

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
