# Public API Reference

PhosPy has **no HTTP API**. The supported surface is the Python API plus the `phospy` CLI.

## Table of Contents

- [Import Guide](#import-guide)
- [API Index](#api-index)
- [Configuration Contracts (`phospy.api`)](#configuration-contracts-phospyapi)
- [Datasets and Preprocessing](#datasets-and-preprocessing)
- [Reference Resolution](#reference-resolution)
- [Kinase Activity Utility](#kinase-activity-utility)
- [Workflow Entry Points](#workflow-entry-points)
- [File-Based Pipeline](#file-based-pipeline)
- [Result Objects](#result-objects)
- [CLI](#cli)
- [Common Exceptions](#common-exceptions)

## Import Guide

Preferred imports:

```python
from phospy.activities import KinaseActivityAnalyzer
from phospy.api import (
    DatasetLoadOptions,
    KinaseActivityConfig,
    KinaseWorkflow,
    PredMatWorkflow,
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.api.workflow_results import (
    KinaseWorkflowResult,
    PredMatWorkflowResult,
    SimpleKinaseWorkflowResult,
)
from phospy.datasets import AnalysisReadyPhosphoDataset, DatasetSchema, PhosphoDataset
from phospy.pipeline import PhosRPipeline
from phospy.prediction import PredMatResult
from phospy.preprocessing import CorePreprocessingConfig
from phospy.references import (
    BundledReferenceProvider,
    ReferenceBundle,
    ReferenceBundleProvenance,
    ReferenceBundleSourceMetadata,
    ReferenceProvider,
)
from phospy.signalomes import SignalomeMapData, SignalomeNetworkData, SignalomeResult
```

## API Index

| Symbol | Import | Purpose |
| --- | --- | --- |
| `SimpleKinaseWorkflow` | `phospy.api` | Recommended end-to-end kinase lane |
| `PredMatWorkflow` | `phospy.api` | Build `predMat` from native workflow inputs |
| `KinaseWorkflow` | `phospy.api` | Native scoring + prediction with intermediates |
| `SignalomeWorkflow` | `phospy.api` | Build signalomes from aligned kinase outputs |
| `PhosphoDataset` | `phospy.datasets` | Validated dataset workspace and preprocessing facade |
| `KinaseActivityAnalyzer` | `phospy.activities` | Activity analysis from matrix + `predMat` |
| `PhosRPipeline` | `phospy.pipeline` | One-shot file loading + preprocessing + optional activity |
| `PredMatResult` | `phospy.prediction` | Stable in-memory/export contract for generated `predMat` |

## Configuration Contracts (`phospy.api`)

### `DatasetLoadOptions`

**Signature**

```python
DatasetLoadOptions(
    phospho_encoding: str | None = None,
    schema: DatasetSchema = DatasetSchema(),
    comparisons: tuple[tuple[str, str], ...] | None = None,
)
```

Use for dataset-loading options in high-level workflows and pipeline constructors.

### `PredictionRunConfig`

**Signature**

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

Use for prediction and workflow tuning (`PredMatWorkflow`, `KinaseWorkflow`, `SimpleKinaseWorkflow`).

### `KinaseActivityConfig`

**Signature**

```python
KinaseActivityConfig(
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
)
```

Use for downstream kinase activity analysis configuration.

### `SignalomeRunConfig`

**Signature**

```python
SignalomeRunConfig(
    kinase_network_threshold: float = 0.9,
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
    module_selection_policy: SignalomeModuleSelectionPolicy = SignalomeModuleSelectionPolicy(),
)
```

Use for signalome construction tuning in `SignalomeWorkflow`.

## Datasets and Preprocessing

### `DatasetSchema`

**Signature**

```python
DatasetSchema(
    total_cols: tuple[str, ...] = ("group1", ..., "group6"),
    phospho_cols: tuple[str, ...] = ("p_group1", ..., "p_group6"),
    corrected_cols: tuple[str, ...] = ("phospho_corrected_1", ..., "phospho_corrected_6"),
)
```

Use this when your sample columns differ from defaults.

### `PhosphoDataset`

**Constructor**

```python
PhosphoDataset(
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    *,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
)
```

**`from_files(...)`**

```python
PhosphoDataset.from_files(
    total_path: str | Path,
    phospho_path: str | Path,
    phospho_encoding: str | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
    schema: DatasetSchema | None = None,
) -> PhosphoDataset
```

**Useful accessors**

- `dataset.total_df_copy`
- `dataset.phospho_df_copy`
- `dataset.total_df_live`
- `dataset.phospho_df_live`
- `dataset.copy_inputs()`

### `dataset.preprocessing.run(...)`

**Signature**

```python
dataset.preprocessing.run(
    *,
    config: CorePreprocessingConfig,
) -> CoreProcessingResult
```

**Notes**

- `config` is required.
- `config.site_matrix_policy` controls duplicate phosphosite collapse.
- `config.max_unmatched_fraction=0.0` is strict mode for protein-correction row loss.

### `dataset.preprocessing.to_analysis_ready(...)`

**Signature**

```python
dataset.preprocessing.to_analysis_ready(
    result: CoreProcessingResult,
    *,
    source: str = "dataset preprocessing",
) -> AnalysisReadyPhosphoDataset
```

### `dataset.preprocessing.run_analysis_ready(...)`

**Signature**

```python
dataset.preprocessing.run_analysis_ready(
    *,
    config: CorePreprocessingConfig,
    source: str = "dataset preprocessing",
) -> AnalysisReadyPhosphoDataset
```

### `PhosphoDataset.run_analysis_ready(...)`

**Signature**

```python
dataset.run_analysis_ready(
    *,
    config: CorePreprocessingConfig,
    source: str = "dataset preprocessing",
) -> AnalysisReadyPhosphoDataset
```

### `AnalysisReadyPhosphoDataset.from_core_processing_result(...)`

**Signature**

```python
AnalysisReadyPhosphoDataset.from_core_processing_result(
    result: CoreProcessingResult,
    *,
    schema: DatasetSchema,
    comparisons: Sequence[tuple[str, str]] | None = None,
    source: str = "core preprocessing",
) -> AnalysisReadyPhosphoDataset
```

### `CorePreprocessingConfig`

**Signature**

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

### `build_analysis_ready_dataset(...)`

**Signature**

```python
build_analysis_ready_dataset(
    *,
    phospho: pd.DataFrame | str | Path,
    preprocessing_config: CorePreprocessingConfig,
    total: pd.DataFrame | str | Path | None = None,
    phospho_encoding: str | None = None,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
    source: str = "analysis ready dataset builder",
    phospho_only_source: str = "analysis ready dataset builder (phospho only)",
) -> AnalysisReadyPhosphoDataset
```

### `CoreOutputWriter.write(...)`

**Signature**

```python
CoreOutputWriter().write(
    result: CoreProcessingResult,
    outdir: str | Path,
    *,
    format: Literal["csv", "tsv", "parquet"] = "csv",
) -> None
```

Core outputs:

- `df_total_unique`
- `df_total_filtered`
- `df_phospho_filtered`
- `df_phospho_corrected`
- `phosr_input`
- `mat_phospho_corrected`
- `site_sequences`

## Reference Resolution

### `ReferenceBundle`

**Signature**

```python
ReferenceBundle(
    *,
    substrate_map: Mapping[str, Sequence[str]],
    motif_sequences: Mapping[str, Sequence[str]],
    species: str,
    source_metadata: ReferenceBundleSourceMetadata,
    provenance: ReferenceBundleProvenance,
)
```

**Rules**

- `substrate_map` and `motif_sequences` must be non-empty.
- Kinase keys must match exactly across both mappings.
- `species`, metadata, and provenance fields must be non-empty.

### `ReferenceProvider`

**Protocol**

```python
class ReferenceProvider(Protocol):
    def resolve(self, *, species: str, reference: str = "auto") -> ReferenceBundle: ...
```

### `BundledReferenceProvider.resolve(...)`

**Signature**

```python
BundledReferenceProvider().resolve(
    *,
    species: str,
    reference: str = "auto",
) -> ReferenceBundle
```

Current bundled support:

- species: `rat`
- references: `auto`, `l6`, `l6_native`
- `auto` resolves to `l6_native`

## Kinase Activity Utility

### `KinaseActivityAnalyzer`

**`load_pred_mat(...)`**

```python
analyzer.load_pred_mat(pred_mat_path: str | Path) -> pd.DataFrame
```

**`run(...)`**

```python
analyzer.run(
    pred_mat: pd.DataFrame | PredMatResult,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult
```

**`write_outputs(...)`**

```python
analyzer.write_outputs(result: KinaseActivityResult, outdir: str | Path) -> None
```

Written files:

- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `kinase_target_table.csv`

## Workflow Entry Points

### `SimpleKinaseWorkflow`

**Constructor**

```python
SimpleKinaseWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = "default",
    *,
    reference_provider: ReferenceProvider | None = None,
    activity_analyzer: KinaseActivityAnalyzer | None = None,
    analysis_ready_builder: AnalysisReadyDatasetBuilder | None = None,
) -> None
```

**`run(...)`**

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

Returns:

- `analysis_ready_dataset`
- `reference_bundle`
- `workflow_result`
- `kinase_activity_result`

### `PredMatWorkflow`

**Constructor**

```python
PredMatWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = "default",
) -> None
```

**`run(...)`**

```python
PredMatWorkflow().run(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]] | None = None,
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
    reference_bundle: ReferenceBundle | None = None,
    prediction_config: PredictionRunConfig | None = None,
) -> PredMatWorkflowResult
```

### `KinaseWorkflow`

**Constructor**

```python
KinaseWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = "default",
) -> None
```

**`run(...)`**

```python
KinaseWorkflow().run(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]] | None = None,
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
    reference_bundle: ReferenceBundle | None = None,
    prediction_config: PredictionRunConfig | None = None,
) -> KinaseWorkflowResult
```

### `SignalomeWorkflow.run(...)`

**Signature**

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

## File-Based Pipeline

### `PhosRPipeline`

**Constructor**

```python
PhosRPipeline(
    dataset: PhosphoDataset,
    pred_mat: pd.DataFrame | PredMatResult | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    activity_config: KinaseActivityConfig | None = None,
    *,
    manifest_writer: RunManifestWriter | None = None,
    output_publisher: OutputPublisher | None = None,
) -> None
```

**`from_files(...)`**

```python
PhosRPipeline.from_files(
    total_path: str | Path,
    phospho_path: str | Path,
    pred_mat_path: str | Path | None = None,
    dataset_options: DatasetLoadOptions | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    activity_config: KinaseActivityConfig | None = None,
    *,
    manifest_writer: RunManifestWriter | None = None,
    output_publisher: OutputPublisher | None = None,
) -> PhosRPipeline
```

**`run(...)`**

```python
pipeline.run(outdir: str | Path | None = None) -> CoreOutputs
```

Returned object fields:

- `core`
- `kinase_activity`

When `outdir` is set, writes pipeline outputs plus `run_manifest.json`.

## Result Objects

### `PredMatResult`

**Constructor**

```python
PredMatResult(data_frame: pd.DataFrame)
```

**Methods**

```python
result.to_frame(copy: bool = True) -> pd.DataFrame
result.to_csv(path: str | Path, index_label: str = "phosphosite") -> Path
```

### `SignalomeResult`

Common accessors:

- `result.modules`
- `result.assignments`
- `result.network`
- `result.expanded_signalomes`
- `result.to_map_data()`
- `result.to_network_data()`

CSV exports include module tables, assignments, and network outputs.

### `SignalomeMapData`

Common methods:

- `modules(copy: bool = True)`
- `sites(copy: bool = True)`
- `kinases(copy: bool = True)`
- `links(copy: bool = True)`
- `to_frames(copy: bool = True)`
- `to_csv(directory: str | Path)`

### `SignalomeNetworkData`

Common methods:

- `nodes(copy: bool = True)`
- `edges(copy: bool = True)`
- `adjacency(copy: bool = True)`
- `to_frames(copy: bool = True)`
- `to_csv(directory: str | Path)`

## CLI

Use `phospy --help` for full option details.

Main options:

- `--total` required TSV path
- `--phospho` required TSV path
- `--outdir` required output directory
- `--pred-mat` optional `predMat` CSV path
- `--phospho-encoding` optional phospho table encoding
- `--localization-threshold` default `0.75`
- `--min-observed` default `4`
- `--total-sentinel` default `10.0`
- `--phospho-sentinel` default `12.0`
- `--max-unmatched-fraction` default `0.0`
- `--kinase-activity-threshold` default `0.6`
- `--kinase-activity-min-substrates` default `3`
- `--kinase-activity-top-n-substrates` default `20`

## Common Exceptions

Most common public exceptions:

- `RequestValidationError` for invalid user/public inputs
- `NoCandidateKinasesError` when prediction thresholds eliminate all candidate kinases

## Related Pages

- [`validation.md`](validation.md)
- [`parity.md`](parity.md)
- [`../examples/predmat_workflow_demo.py`](../examples/predmat_workflow_demo.py)
- [`../examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py)
