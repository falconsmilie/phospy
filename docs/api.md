# API Guide

PhosPy does not have an HTTP API. The supported public surface is the Python API plus the `phospy` CLI.

This page starts with the simplest supported path, then links to the lower-level options.

## Which Entry Point Should You Use?

| Goal | Start here |
| --- | --- |
| Clean and prepare phospho data | `PhosphoDataset` |
| Run the common end-to-end workflow | `SimpleKinaseWorkflow` |
| Build a signalome from workflow outputs | `SignalomeWorkflow` |

`PredMatWorkflow`, `KinaseWorkflow`, and `PhosRPipeline` are internal orchestration helpers and are not part of the supported public API.

## Fastest Path for Most Users

```python
from phospy.api import PredictionRunConfig, SimpleKinaseWorkflow

with SimpleKinaseWorkflow().run(
    phospho="study_phospho.tsv",
    total="study_total.tsv",
    species="rat",
    prediction_config=PredictionRunConfig(
        min_substrates=1,
        min_motif_size=1,
        ensemble_size=2,
        top=3,
        inclusion=2,
        n_iterations=2,
        random_state=7,
    ),
) as result:
    pred_mat = result.pred_mat_result.to_frame(copy=False)
    weighted_activity = result.kinase_activity_result.weighted_activity
```

The returned `result` includes:

- `analysis_ready_dataset`
- `reference_bundle`
- `scoring_result`
- `prediction_result`
- `pred_mat_result`
- `kinase_activity_result`

## Common Recipes

### Preprocessing only

```python
from phospy.datasets import PhosphoDataset
from phospy.preprocessing import CorePreprocessingConfig

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
analysis_ready = dataset.preprocessing.run_analysis_ready(
    config=CorePreprocessingConfig(max_unmatched_fraction=0.1)
)
```

### Signalome analysis from workflow outputs

The next example assumes you already have a `result` from `SimpleKinaseWorkflow.run(...)`.

```python
from phospy.api import SignalomeRunConfig, SignalomeWorkflow

signalome = SignalomeWorkflow().run_from_analysis_ready(
    dataset=result.analysis_ready_dataset,
    scoring_result=result.scoring_result,
    prediction_result=result.prediction_result,
    kinases_of_interest=list(result.pred_mat_result.kinase_names[:2]),
    config=SignalomeRunConfig(signalome_cutoff=0.5),
)
```

### Frame mutability contract

`PhosphoDataset` is a mutable workspace, but safe read access is detached by default:

- `dataset.inputs`, `dataset.total_df_copy`, `dataset.phospho_df_copy`, `dataset.copy_inputs()`

Unsafe mutable access is advanced and intentional:

- `dataset.to_mutable_frames_unsafe()`
- Warning: mutating returned frames mutates the dataset's owned internal state

`SignalomeResult` and nested signalome wrappers (`modules`, `assignments`, `network`) also default to detached reads:

- table properties (for example `signalome.site_assignments`, `signalome.signalome_modules`)
- `to_frames()` (always detached)

Unsafe mutable access is advanced and intentional:

- `signalome.to_mutable_frames_unsafe()`
- `signalome.modules.to_mutable_tables_unsafe()`
- `signalome.assignments.to_mutable_tables_unsafe()`
- `signalome.network.to_mutable_state_unsafe()`
- Warning: mutating returned objects mutates the signalome result's owned internal state

## Configuration Objects

### `DatasetLoadOptions` (`phospy.api`)

```python
DatasetLoadOptions(
    phospho_encoding: str | None = None,
    schema: DatasetSchema = DatasetSchema(),
    comparisons: tuple[tuple[str, str], ...] | None = None,
)
```

Use this with `SimpleKinaseWorkflow.run(..., dataset_options=...)` when you need custom input encoding, a custom schema, or explicit comparison pairs.

### `CorePreprocessingConfig` (`phospy.preprocessing`)

```python
CorePreprocessingConfig(
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    max_unmatched_fraction: float = 0.0,
    site_matrix_policy: SiteMatrixPolicy = SiteMatrixPolicy(),
)
```

Use this to control filtering, missing-value sentinels, protein correction tolerance, and site-matrix behaviour.

### `SiteMatrixPolicy` (`phospy.preprocessing`)

```python
SiteMatrixPolicy(
    duplicate_site_strategy: DuplicateSiteStrategy = "max_mean_signal",
    missing_data_policy: SiteMatrixMissingDataPolicy = "drop_any_missing",
    minimum_observed_values: int | None = None,
)
```

Missing-data options:

- `drop_any_missing`: keep only complete corrected rows
- `retain_missing`: keep partially observed rows and preserve `NaN`
- `require_min_observed_values`: keep rows with at least a minimum number of observed values

### `PredictionRunConfig` (`phospy.api`)

```python
PredictionRunConfig(
    min_substrates: int = 1,
    min_motif_size: int = 1,
    allow_profile_only_fallback: bool = False,
    ensemble_size: int = 10,
    top: int = 50,
    score_threshold: float = 0.8,
    inclusion: int = 20,
    n_iterations: int = 5,
    random_state: int | None = None,
    svm_mode: PredictionSvmMode | None = None,
    profile_policy: KinaseProfilePolicy = KinaseProfilePolicy(),
)
```

Use this to control candidate filtering, scoring thresholds, sampling, and prediction mode.

### `KinaseActivityConfig` (`phospy.api`)

```python
KinaseActivityConfig(
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
)
```

Use this to control weighted activity and KSEA-style downstream summaries.

### `SignalomeRunConfig` (`phospy.api`)

```python
SignalomeRunConfig(
    kinase_network_threshold: float = 0.9,
    kinase_network_policy: SignalomeKinaseNetworkPolicy = "positive_only",
    assignment_policy: SignalomeAssignmentPolicy = "cutoff_binary",
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
    module_selection_policy: SignalomeModuleSelectionPolicy = SignalomeModuleSelectionPolicy(),
)
```

Use this to control signalome assignment, kinase network construction, and module selection.

## Workflow Signatures

### `SimpleKinaseWorkflow`

```python
SimpleKinaseWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = "default",
    *,
    reference_provider: ReferenceProvider | None = None,
    activity_analyzer: KinaseActivityAnalyzer | None = None,
    analysis_ready_builder: AnalysisReadyDatasetBuilder | None = None,
)
```

```python
SimpleKinaseWorkflow().run(
    *,
    phospho: pd.DataFrame | str | Path,
    species: str,
    total: pd.DataFrame | str | Path | None = None,
    reference: str = "auto",
    dataset_options: DatasetLoadOptions | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    prediction_config: PredictionRunConfig | None = None,
    activity_config: KinaseActivityConfig | None = None,
) -> SimpleKinaseWorkflowResult
```

Use a context manager when practical:

```python
with SimpleKinaseWorkflow().run(...) as result:
    ...
```

### `SignalomeWorkflow`

```python
SignalomeWorkflow().run(
    *,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult | PredMatResult,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    site_to_protein: Mapping[str, str] | None = None,
    config: SignalomeRunConfig | None = None,
) -> SignalomeResult
```

```python
SignalomeWorkflow().run_from_analysis_ready(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult | PredMatResult,
    kinases_of_interest: Sequence[str],
    site_to_protein: Mapping[str, str] | None = None,
    metadata_protein_columns: Sequence[str] | None = None,
    metadata_fallback_policy: str = "strict",
    allow_gene_symbol_fallback: bool = False,
    allow_ambiguous_metadata_mapping: bool = False,
    config: SignalomeRunConfig | None = None,
) -> SignalomeResult
```

If `site_to_protein` is omitted, `run(...)` falls back to supported `ENTITY;SITE;` phosphosite IDs such as `BTK;Y551;`.

`run_from_analysis_ready(...)` defaults to strict metadata resolution and requires a `protein_id` metadata column. To opt in to metadata fallback columns, set `metadata_fallback_policy="metadata"` and provide `metadata_protein_columns`. Gene-symbol fallback is disabled by default and must be explicitly enabled with `allow_gene_symbol_fallback=True`.

## Preferred Imports

```python
from phospy.api import (
    DatasetLoadOptions,
    KinaseActivityConfig,
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.api.workflow_results import SimpleKinaseWorkflowResult
from phospy.datasets import AnalysisReadyPhosphoDataset, DatasetSchema, PhosphoDataset
from phospy.preprocessing import CorePreprocessingConfig
```

## CLI Reference

Use `phospy --help` for the full option list.

Main options:

- `--total` required TSV path
- `--phospho` required TSV path
- `--outdir` required output directory
- `--pred-mat` optional `predMat` CSV path
- `--phospho-encoding` optional phospho table encoding
- `--localization-threshold` minimum localisation probability
- `--min-observed` minimum observed values per row
- `--total-sentinel` total-table sentinel value treated as missing
- `--phospho-sentinel` phospho-table sentinel value treated as missing
- `--kinase-activity-threshold` threshold for downstream kinase activity summaries
- `--kinase-activity-min-substrates` minimum substrate count for downstream summaries
- `--kinase-activity-top-n-substrates` top-N substrates for weighted activity summaries
- `--max-unmatched-fraction` allowed fraction of phosphosite rows without matching protein rows during correction

## Common Exceptions

- `RequestValidationError`: invalid public input or configuration
- `NoCandidateKinasesError`: thresholds removed all candidate kinases
- `InputCompatibilityError`: inputs are valid on their own but incompatible together

## Related Docs

- [`validation.md`](validation.md)
- [`parity.md`](parity.md)
- [`../examples/simple_workflow_demo.py`](../examples/simple_workflow_demo.py)
- [`../examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py)
