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

For concise scientist-facing assumptions and interpretation notes, see
[Workflow Contracts](workflow_contracts.md).

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

For most users, start with intent presets and only drop to low-level policy
constants when you need fine-grained control:

```python
from phospy.api import DatasetPreprocessingConfig

strict = DatasetPreprocessingConfig.strict()
raw_table = DatasetPreprocessingConfig.from_raw_phosphosite_table()
```

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

Quantitative meaning is explicit after preprocessing:

- fully corrected log2 datasets: `phospho_total_log_ratio`
- log2 datasets without total-protein correction: `phosphosite_log_abundance`
- partially corrected datasets (`unmatched_policy="allow_uncorrected"` with retained unmatched rows):
  `mixed_phospho_total_log_ratio_and_phosphosite_log_abundance`

For mixed datasets, row-level correction status is available in
`dataset.processing_state.total_protein_correction.diagnostics`, including
`corrected_phosphosite_row_ids`, `corrected_phosphosite_to_total_protein_row_id`,
and `uncorrected_phosphosite_row_reasons`.

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

Use presets first:

```python
from phospy.api import KinaseScoringConfig

default_scoring = KinaseScoringConfig.default()
strict_missing = KinaseScoringConfig.strict_missing_values()
```

| Field | Meaning |
| --- | --- |
| `min_substrates` | minimum quantified substrates per kinase; must be at least `2` |
| `include_diagnostic_scoring_tables` | include extra scoring diagnostics when `True` |
| `profile_missing_value_strategy` | `"strict"` or `"median_skipna"` |
| `allow_mixed_total_protein_quantitative_meaning` | default `False`; require explicit opt-in to run on mixed corrected/uncorrected datasets |

### `KinasePredictionConfig`

Use intent presets for the two common prediction lanes:

```python
from phospy.api import KinasePredictionConfig

deterministic = KinasePredictionConfig.deterministic()
adaptive = KinasePredictionConfig.adaptive_reproducible(random_state=1)
```

| Field | Default | Meaning |
| --- | --- | --- |
| `top_k` | `30` | top predicted substrate sites per kinase |
| `deterministic_max_selected_kinases` | `10` | retained kinases in deterministic mode |
| `adaptive_ensemble_runs` | `10` | ensemble runs in adaptive mode |
| `mode` | `"deterministic_ranking"` | `"deterministic_ranking"` or `"adaptive_ensemble"` |
| `adaptive_policy` | `"stable"` | `"stable"` or `"r_parity"` |
| `n_iterations` | `5` | adaptive sampling iterations |
| `random_state` | `None` | required when `mode="adaptive_ensemble"`; unused in deterministic mode |

`adaptive_policy="r_parity"` exists for parity/replay validation lanes. The
default recommendation for production workflows is `adaptive_policy="stable"`.

### `KinaseActivityConfig`

| Field | Default | Meaning |
| --- | --- | --- |
| `enabled` | `True` | set `False` to disable activity after construction |
| `method` | `"simplified_weighted_substrate_activity"` | activity method selector: weighted heuristic or KSEA-style z-score |
| `threshold` | `0.6` | prediction-score threshold for selected substrates |
| `min_substrates` | `3` | minimum selected substrates per kinase |
| `top_n_substrates` | `20` | top predicted substrates used in weighted activity |
| `ksea_min_substrates` | `5` | minimum substrates per kinase/condition for KSEA scoring |
| `ksea_evidence_threshold` | `None` | KSEA membership threshold; defaults to `threshold` when `None` |
| `ksea_p_value_method` | `"normal_approximation"` | KSEA p-value method |
| `ksea_adjust_p_values` | `True` | apply Benjamini-Hochberg q-values per condition for KSEA |

`method="ksea_zscore"` enables KSEA-style z-score activity inference with
per-condition computability statuses and statistics output. It is not equivalent
to PhosR kinase activity inference.

When activity is enabled, `result.activity_result.activity_method` exposes
explicit method identity metadata, including:

- weighted method:
  - `activity_method_id="simplified_weighted_substrate_activity_v1"`
  - `activity_method_family="heuristic_weighted_substrate_score"`
  - `is_ksea=False`
  - `is_phosr_kinase_activity_equivalent=False`
- KSEA method:
  - `activity_method_id="ksea_zscore_v1"`
  - `activity_method_family="substrate_set_enrichment"`
  - `is_ksea=True`
  - `is_phosr_kinase_activity_equivalent=False`

## Signalome Config

`SignalomeConfig` is grouped by user intent:

- `scientific=SignalomeScientificConfig(...)`
- `clustering=SignalomeClusteringConfig(...)`
- `validation=SignalomeValidationConfig(...)`
- `output=SignalomeOutputConfig(...)`
- `performance=SignalomePerformanceConfig(...)` (advanced guardrails)

Basic usage stays simple:

```python
result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
```

or:

```python
result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig(),
    )
)
```

Intent presets for common signalome behaviour:

```python
from phospy.api import SignalomeConfig

strict = SignalomeConfig.strict()
permissive = SignalomeConfig.permissive_missing_scores()
sampled = SignalomeConfig.sampled_candidate_scoring()
```

`SignalomeConfig.sampled_candidate_scoring()` changes candidate module-count
scoring only and keeps scale guards enabled.

Grouped options:

| Group | Option | Default | Meaning |
| --- | --- | --- | --- |
| `scientific` | `substrate_support_cutoff` | `0.5` | prediction support cutoff for kinase-supported substrates |
| `scientific` | `assignment_policy` | `"cutoff_binary"` | `"cutoff_binary"` or `"weighted_top"` |
| `clustering` | `module_count` | `None` | explicit module count; omit for automatic selection |
| `clustering` | `module_selection_primary_correlation_threshold` | `0.5` | first threshold for automatic module selection |
| `clustering` | `module_selection_fallback_correlation_threshold` | `0.1` | fallback threshold for automatic module selection |
| `clustering` | `module_selection_max_clusters` | `10` | largest candidate module count considered |
| `clustering` | `candidate_scoring_policy` | `"full"` | `"full"` or `"sampled"` for candidate module-count scoring |
| `clustering` | `clustering_engine` | `"scipy_hierarchical"` | `"exact_python"` or `"scipy_hierarchical"` |
| `validation` | `score_preconditioning_policy` | `"error_on_drop"` | strict-by-default handling of all-missing downstream score rows; set `"allow_and_report"` to opt in to row dropping |
| `validation` | `allow_mixed_total_protein_quantitative_meaning` | `False` | require explicit opt-in to run on mixed corrected/uncorrected datasets |
| `output` | `network_correlation_threshold` | `0.5` | threshold used by the network policy |
| `output` | `network_policy` | `"signed"` | `"positive_only"`, `"absolute_threshold"`, or `"signed"` |
| `performance` | `max_exact_tree_sites` | `2000` | advanced hard guard for exact tree construction |
| `performance` | `max_full_candidate_scoring_sites` | `2000` | advanced hard guard for full candidate scoring |

`candidate_scoring_policy="sampled"` reduces candidate module-count scoring
cost only. It does not remove exact tree construction, and it does not make
the full signalome workflow approximate.

`clustering.clustering_engine="scipy_hierarchical"` is the preferred production
default. Use `"exact_python"` mainly for reference/debug checks.

Tree implementation details are recorded in internal diagnostics/provenance as
backend metadata, not as a public configuration choice.

Low-level constants and nested config fields remain part of the advanced API for
specialized scientific workflows.

Signalome module/network scores are derived summaries over upstream kinase score
profiles. They are not probabilities, calibrated confidence values, or direct
proof of causal regulation.

## Result Models

### DataFrame Ownership

PhosPy owns validated datasets and workflow results internally.

When you pass a DataFrame into PhosPy, the accepted dataset/table keeps its own
copy so later changes to your original DataFrame do not silently affect
analysis.

When you retrieve result tables through the public API, use export helpers
(`to_dataframe(...)`, `to_pandas(...)`, or `*_dataframe(...)`). They always
return safe snapshot copies. Mutating an exported DataFrame does not mutate the
owning dataset/result object.

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

Use `dataset.to_dataframe()` for a safe phospho snapshot.

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

Use `to_dataframe()` on scoring/prediction/activity result objects for safe table
snapshots.

`result.provenance.scientific_policies` lists the active scientific scoring
policies with stable IDs, assumptions, parameters, and output-scale notes for
auditability.

`result.provenance.workflow_parameters["activity_config"]["activity_method"]`
mirrors this method identity, and `activity_method_summary` reports method-level
computed/skipped counts.

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

Use `result.to_dataframe()` for `expanded_signalome`, and
`site_membership_dataframe()` / `protein_site_context_dataframe()` for safe
sidecar snapshots.

Inspect runtime clustering diagnostics from provenance:

```python
scale_guard = result.provenance.workflow_parameters["scale_guard"]
print(scale_guard["clustering_engine"])
print(scale_guard["tree_generation_mode"])          # full_exact_tree_construction
print(scale_guard["tree_generation_is_approximate"])  # False
print(scale_guard["candidate_scoring_mode"])        # full / sampled / not_evaluated
print(scale_guard["candidate_scoring_is_approximate"])
print(scale_guard["max_exact_tree_sites"])
print(scale_guard["max_full_candidate_scoring_sites"])
print(scale_guard["candidate_scoring_sampled_site_total"])
print(scale_guard["candidate_scoring_sampled_pair_count"])
```

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
