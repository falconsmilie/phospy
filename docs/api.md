# API Guide

This guide covers the supported public Python API. PhosPy does not expose HTTP
endpoints.

## Import Contract

Use top-level `phospy` for the main entrypoints:

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

Use `phospy.api` for requests, configs, results, enums, references, and public
exceptions:

```python
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
)
```

All public executors use `run(request)`.

## Public Workflow Shape

1. `DatasetBuildRequest` -> `AnalysisReadyDatasetBuilder.run(...)` -> `AnalysisReadyPhosphoDataset`

2. `KinaseWorkflowRequest` -> `KinaseWorkflow.run(...)` -> `KinaseWorkflowResult`

3. `SignalomeWorkflowRequest` -> `SignalomeWorkflow.run(...)` -> `SignalomeWorkflowResult`

The beginner lane is rat-first because bundled runtime references in `1.5.0` are
rat-only. Human and mouse workflows need an explicit `ReferenceBundle`.

## Request Models

### `DatasetBuildRequest`

| Field | Meaning |
| --- | --- |
| `phospho` | pandas `DataFrame` or supported file path for the site-by-sample intensity matrix |
| `site_metadata` | pandas `DataFrame` or path aligned to `phospho.index` |
| `sample_metadata` | optional pandas `DataFrame` or path aligned to sample columns |
| `total` | optional total-protein matrix used by total-protein correction |
| `organism` | optional `Organism` enum value |
| `preprocessing_config` | `DatasetPreprocessingConfig` |

Supported file suffixes are `.csv`, `.tsv`, `.txt` as tab-separated text, and
`.parquet`. CSV/TSV/TXT inputs are read with the first column as the row index.

Supported site-metadata aliases are deliberately narrow:

- `gene_name` may stand in for `gene_symbol`
- `centralized_sequence` may stand in for `site_sequence`

The builder may derive `gene_symbol` and `site` from index values formatted like
`MAPK14;Y182;`. It does not derive `protein_id`.

### `KinaseWorkflowRequest`

| Field | Meaning |
| --- | --- |
| `dataset` | an `AnalysisReadyPhosphoDataset` |
| `references` | `ReferencePreset` or `ReferenceBundle`; default is `ReferencePreset.AUTO` |
| `scoring_config` | `KinaseScoringConfig` |
| `prediction_config` | `KinasePredictionConfig` |
| `activity_config` | `KinaseActivityConfig` or `None`; `None` disables activity output |

### `SignalomeWorkflowRequest`

| Field | Meaning |
| --- | --- |
| `kinase_result` | a `KinaseWorkflowResult` |
| `config` | `SignalomeConfig` |

Signalome requires explicit, non-empty `dataset.site_metadata.protein_id` for all
interpreted sites.

## Dataset Preprocessing Configs

`DatasetPreprocessingConfig` groups six areas.

| Field | Config class | Default |
| --- | --- | --- |
| `intensity_transform` | `DatasetIntensityTransformConfig` | `policy="identity"` |
| `normalisation` | `DatasetNormalisationConfig` | `policy="none"` |
| `missing_data` | `DatasetMissingDataConfig` | `policy="forbid"` |
| `total_protein_correction` | `DatasetTotalProteinCorrectionConfig` | `policy="none"` |
| `site_matrix` | `DatasetSiteMatrixConfig` | `policy="as_input"` |
| `comparisons` | `DatasetComparisonBuildingConfig` | `policy="none"` |

### Intensity Transform

`DatasetIntensityTransformConfig` supports:

- `policy="identity"`
- `policy="log2"`, applying `log2(value + pseudocount)`
- `pseudocount`, which must be non-negative

Use `log2` only when values are valid for the chosen pseudocount.

### Normalisation

`DatasetNormalisationConfig` supports:

- `policy="none"`
- `policy="median_center"`
- `policy="quantile"`

Quantile normalisation is dense and sort-heavy. Use it only when matched sample
distributions are scientifically appropriate.

### Missing Data

`DatasetMissingDataConfig` supports:

- `policy="forbid"`
- `policy="impute_row_median"`
- `min_observed_values`, used with row-median imputation

The dataset that leaves the builder must be missing-value-free.

### Total-Protein Correction

`DatasetTotalProteinCorrectionConfig` supports:

- `policy="none"`
- `policy="subtract_log_total"`
- `identity=DatasetTotalProteinCorrectionIdentityConfig(...)`

`subtract_log_total` requires `intensity_transform.policy="log2"`, a `total`
table aligned to phospho sample columns, and an explicit identity mapping.

Identity mapping supports:

- `mode="direct"`: match `site_metadata[phosphosite_key]` directly to total IDs;
  use `total_protein_key="__index__"` for total index matching
- `mode="mapping_table"`: provide `mapping_table`, `mapping_phosphosite_key`,
  and `mapping_total_protein_key`

Strictness controls are:

- `duplicate_policy="error"`
- `unmatched_policy="error"` or `"allow_uncorrected"`

Gene-symbol matching is convenient, but it is not a universal biological identity
guarantee. Prefer accessions, protein-group IDs, or an explicit mapping table
when that better matches the experiment.

### Site Matrix Construction

`DatasetSiteMatrixConfig` supports:

- `policy="as_input"`
- `policy="build_from_metadata"`
- `duplicate_site_policy`: `"max_mean_signal"`, `"first"`, `"aggregate_mean"`,
  `"aggregate_median"`, or `"error"`
- `missing_data_policy="drop_any_missing"`

Duplicate-site handling is a scientific choice. `error` is cautious. `first` is
simple but input-order dependent. `max_mean_signal` favours the strongest row.
Aggregate policies collapse rows numerically and can blur peptide context.

### Comparison Building

`DatasetComparisonBuildingConfig` supports:

- `policy="none"`
- `policy="sample_metadata_pairs"`
- `sample_group_column`, defaulting to `comparison_group`
- optional explicit `pairs`

When `sample_metadata_pairs` is used, `sample_metadata` is required.

## Kinase Configs

### `KinaseScoringConfig`

| Field | Meaning |
| --- | --- |
| `min_substrates` | minimum quantified substrates per kinase; must be at least `2` |
| `include_diagnostic_scoring_tables` | include extra scoring diagnostics when `True` |
| `profile_missing_value_strategy` | `"strict"` or `"median_skipna"` |

### `KinasePredictionConfig`

| Field | Default | Meaning |
| --- | --- | --- |
| `top_k` | `30` | top predicted substrate sites per kinase |
| `deterministic_max_selected_kinases` | `10` | retained kinases in deterministic mode |
| `adaptive_ensemble_runs` | `10` | ensemble runs in adaptive mode |
| `mode` | `"deterministic_ranking"` | `"deterministic_ranking"` or `"adaptive_ensemble"` |
| `adaptive_policy` | `"stable"` | `"stable"` or `"r_parity"` |
| `n_iterations` | `5` | adaptive sampling iterations |
| `random_state` | `None` | optional random state |

### `KinaseActivityConfig`

| Field | Default | Meaning |
| --- | --- | --- |
| `enabled` | `True` | set `False` to disable activity after construction |
| `threshold` | `0.6` | prediction-score threshold for selected substrates |
| `min_substrates` | `3` | minimum selected substrates per kinase |
| `top_n_substrates` | `20` | top predicted substrates used in weighted activity |

`thresholded_substrate_mean_activity` is a simple mean phospho signal over
predicted substrates above the configured threshold. It is not full KSEA-style
enrichment.

## Signalome Config

`SignalomeConfig` supports these fields:

| Field | Default | Meaning |
| --- | --- | --- |
| `substrate_support_cutoff` | `0.5` | prediction support cutoff for kinase-supported substrates |
| `network_correlation_threshold` | `0.5` | threshold used by the network policy |
| `network_policy` | `"signed"` | `"positive_only"`, `"absolute_threshold"`, or `"signed"` |
| `assignment_policy` | `"cutoff_binary"` | `"cutoff_binary"` or `"weighted_top"` |
| `score_preconditioning_policy` | `"allow_and_report"` | allow or reject all-missing downstream score rows |
| `module_count` | `None` | explicit module count; omit for automatic selection |
| `module_selection_primary_correlation_threshold` | `0.5` | first threshold for automatic module selection |
| `module_selection_fallback_correlation_threshold` | `0.1` | fallback threshold for automatic module selection |
| `module_selection_max_clusters` | `10` | largest candidate module count considered |
| `tree_engine` | `"exact"` | exact tree construction; this is the only public value |
| `candidate_scoring_policy` | `"full"` | `"full"` or `"sampled"` candidate module-count scoring |
| `max_exact_tree_sites` | `2000` | hard guard for exact tree construction |
| `max_full_candidate_scoring_sites` | `2000` | hard guard for full candidate scoring |
| `clustering_engine` | `"scipy_hierarchical"` | `"exact_python"` or `"scipy_hierarchical"` |

`candidate_scoring_policy="sampled"` reduces candidate module-count scoring
cost. It does not remove the exact tree-construction guard.

`clustering_engine="scipy_hierarchical"` is the preferred production default.
Use `clustering_engine="exact_python"` mainly for reference/debug checks.

Signalome module/network scores are derived summaries over upstream kinase score
profiles. They are not probabilities, calibrated confidence values, or direct
proof of causal regulation.

## Result Models

### `AnalysisReadyPhosphoDataset`

Important fields:

- `phospho`
- `site_metadata`
- `sample_metadata`
- `total`
- `comparisons`
- `organism`
- `intensity_scale_state`
- `processing_state`
- `preprocessing_report`
- `provenance`

Read `intensity_scale_state.label` together with
`intensity_scale_state.quantity`. For example, `log2` describes numeric scale;
`phospho_total_log_ratio` describes what the values mean scientifically.

### `KinaseWorkflowResult`

Important fields:

- `dataset`
- `references`
- `scoring_result`
- `prediction_result`
- `activity_result`
- `provenance`

Common tables include `profile_scores`, `rank_weighted_fusion_scores`,
`pred_mat`, and activity tables when activity is enabled.

`result.provenance.scientific_policies` lists the active scientific scoring
policies with stable IDs, assumptions, parameters, and output-scale notes for
auditability.

### `SignalomeWorkflowResult`

Important fields:

- `dataset`
- `kinase_result`
- `module_assignments`
- `signalome_modules`
- `kinase_network`
- `module_selection_diagnostics`
- `score_preconditioning_diagnostics`
- `expanded_signalome`
- `site_membership`
- `protein_site_context`
- `provenance`

Undefined kinase correlations are preserved as missing values. A correlation of
`0.0` means a correlation was estimated and is near zero.

`result.provenance.scientific_policies` includes signalome-specific scientific
policies (for candidate module-count scoring and protein-module derivation) and
their resolved parameters.

## References

`Organism` values are `"human"`, `"mouse"`, and `"rat"`.

`ReferencePreset` values are `"auto"`, `"human"`, `"mouse"`, and `"rat"`.
Enum presence does not mean bundled runtime data exists for every organism in
this release.

Use `ReferenceBundle` for custom references. It requires:

- `organism`
- `kinase_substrate_map` with `kinase` and `substrate_site`
- `site_sequences` indexed by site ID with `site_sequence`

## Public Exceptions

All user-facing exception types are available from `phospy.api`. Common ones are:

- `PhosPyInputError`
- `UnsupportedInputFormatError`
- `PhosPyValidationError`
- `ReferenceResolutionError`
- `ReferenceCompatibilityError`
- `WorkflowValidationError`
- `WorkflowBoundaryError`
- `SignalomeScaleError`

## Small Working Example

```python
import pandas as pd

from phospy import (
  AnalysisReadyDatasetBuilder,
  KinaseWorkflow
)
from phospy.api import (
  DatasetBuildRequest, 
  KinaseWorkflowRequest,
  Organism,
  ReferencePreset
)

phospho = pd.DataFrame(
    {
      "sample_a": [1.00, 0.70], 
      "sample_b": [1.10, 0.80], 
      "sample_c": [0.95, 0.75]
    },
    index=[
      "TSC2;S939;", 
      "GSK3B;S9;"
    ],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B"],
        "site": ["S939", "S9"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "_______MSGRPRTTSFAESCKPVQQPSAFG",
        ],
        "protein_id": ["TSC2", "GSK3B"],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
      phospho=phospho, 
      site_metadata=site_metadata, 
      organism=Organism.RAT
    )
)
result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
      dataset=dataset, 
      references=ReferencePreset.AUTO, 
      activity_config=None
    )
)
print(result.prediction_result.pred_mat)
```
