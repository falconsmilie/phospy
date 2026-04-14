# API Guide

PhosPy has no HTTP API. The supported public surface is the Python API plus the `phospy` CLI.

## Supported Entry Points

| Task | Supported public entry point |
| --- | --- |
| Preprocessing | `PhosphoDataset` |
| Kinase scoring and prediction | `SimpleKinaseWorkflow` |
| Signalome analysis | `SignalomeWorkflow` |

`PredMatWorkflow`, `KinaseWorkflow`, and `PhosRPipeline` are internal orchestration helpers and are not part of the supported public API.

## Quick Start

```python
from phospy.api import PredictionRunConfig, SimpleKinaseWorkflow

result = SimpleKinaseWorkflow().run(
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
)
```

`result` gives you:

- `analysis_ready_dataset`
- `reference_bundle`
- `scoring_result`
- `prediction_result`
- `pred_mat_result`
- `kinase_activity_result`

## Task Recipes

### 1. Build an Analysis-Ready Dataset

```python
from phospy.datasets import PhosphoDataset
from phospy.preprocessing import CorePreprocessingConfig

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
analysis_ready = dataset.preprocessing.run_analysis_ready(
    config=CorePreprocessingConfig(max_unmatched_fraction=0.1)
)
```

### 2. Run Signalome Analysis from `SimpleKinaseWorkflow` Outputs

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

## Configuration Reference

### `DatasetLoadOptions` (`phospy.api`)

```python
DatasetLoadOptions(
    phospho_encoding: str | None = None,
    schema: DatasetSchema = DatasetSchema(),
    comparisons: tuple[tuple[str, str], ...] | None = None,
)
```

Used by `SimpleKinaseWorkflow.run(..., dataset_options=...)`.

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

`SiteMatrixPolicy` controls scientifically important phosphosite matrix decisions:

```python
SiteMatrixPolicy(
    duplicate_site_strategy: DuplicateSiteStrategy = "max_mean_signal",
    missing_data_policy: SiteMatrixMissingDataPolicy = "drop_any_missing",
    minimum_observed_values: int | None = None,
)
```

Missing-data policy trade-offs:

- `drop_any_missing`: strict complete-case matrix rows (default; reproducible legacy behavior).
- `retain_missing`: keeps partially observed rows, preserving coverage but allowing `NaN` in the site matrix.
- `require_min_observed_values`: compromise between coverage and completeness; set `minimum_observed_values` (for example `2`) to require at least that many observed corrected columns per row.

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

Used by `SimpleKinaseWorkflow`.

### `KinaseActivityConfig` (`phospy.api`)

```python
KinaseActivityConfig(
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
)
```

Used by `SimpleKinaseWorkflow`.

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

Used by `SignalomeWorkflow`.

## Workflow Reference

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
    config: SignalomeRunConfig | None = None,
) -> SignalomeResult
```

If `site_to_protein` is omitted, signalome grouping falls back to parsing canonical
site IDs in `ENTITY;SITE;` format (for example `BTK;Y551;`). For non-canonical IDs,
pass `site_to_protein` explicitly or use `run_from_analysis_ready(...)` to resolve
from dataset metadata.

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
- `--phospho-encoding` optional phospho table encoding; defaults to UTF-8 when omitted

## Common Exceptions

- `RequestValidationError`: invalid public input or configuration
- `NoCandidateKinasesError`: thresholds removed all candidate kinases
- `InputCompatibilityError`: inputs are valid individually but incompatible together

## Related Docs

- [`validation.md`](validation.md)
- [`parity.md`](parity.md)
- [`../examples/simple_workflow_demo.py`](../examples/simple_workflow_demo.py)
- [`../examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py)
