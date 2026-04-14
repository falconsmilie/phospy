# API Guide

PhosPy has **no HTTP API**. The supported public surface is the Python API plus the `phospy` CLI.

This page is layered:

- start with **Quick Start** if you are new
- use **Task Recipes** for common jobs
- use **Reference** sections for signatures and advanced tuning

## Choose Your Entry Point

| You want to... | Use | Why |
| --- | --- | --- |
| Run the common end-to-end path from user inputs | `SimpleKinaseWorkflow` | Fewest moving parts; best default |
| Build only a `predMat` from workflow-shaped inputs | `PredMatWorkflow` | Focused output contract (`PredMatResult`) |
| Get intermediate scoring and prediction internals | `KinaseWorkflow` | Full native result bundle |
| Build downstream signalomes | `SignalomeWorkflow` | Signalome-specific validated boundary |
| Process files and optionally run activity analysis | `PhosRPipeline` or `phospy` CLI | File-first orchestration + output publishing |
| Analyze kinase activity from an existing `predMat` | `KinaseActivityAnalyzer` | Fastest downstream-only lane |

## Quick Start

If you already have phospho data (and optionally total data), this is the shortest supported path:

```python
from phospy.api import SimpleKinaseWorkflow

result = SimpleKinaseWorkflow().run(
    phospho="study_phospho.tsv",
    total="study_total.tsv",
    species="rat",
    reference="auto",
)

pred_mat_df = result.pred_mat_result.to_frame()
weighted_activity = result.kinase_activity_result.weighted_activity
```

`result` includes:

- `analysis_ready_dataset`
- `reference_bundle`
- `workflow_result`
- `kinase_activity_result`

Bundled reference support is intentionally narrow:

- species: `rat`
- references: `auto`, `l6`, `l6_native`
- `auto` resolves to `l6_native`

Runnable demo: [`../examples/simple_workflow_demo.py`](../examples/simple_workflow_demo.py)

## Core Concepts

- `PhosphoDataset`: validated workspace for total + phospho tables
- `AnalysisReadyPhosphoDataset`: post-preprocessing boundary used by inference workflows
- `ReferenceBundle`: validated kinase prior package (`substrate_map` + `motif_sequences`)
- `PredMatResult`: stable prediction matrix contract
- `KinaseActivityResult`: weighted activity, KSEA outputs, and target summaries
- `SignalomeResult`: module, assignment, and network outputs

## Task Recipes

### 1) Build an Analysis-Ready Dataset

```python
from phospy.datasets import PhosphoDataset
from phospy.preprocessing import CorePreprocessingConfig

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
analysis_ready = dataset.run_analysis_ready(
    config=CorePreprocessingConfig(max_unmatched_fraction=0.1)
)
```

Use this when you want to inspect preprocessing outputs before prediction.

### 2) Generate `predMat`

```python
from phospy.api import PredMatWorkflow, PredictionRunConfig

result = PredMatWorkflow().run(
    phospho_matrix=phospho_matrix,
    substrate_map=substrate_map,
    site_sequences=site_sequences,
    motif_sequences=motif_sequences,
    prediction_config=PredictionRunConfig(min_substrates=2, min_motif_size=2),
)

pred_mat = result.pred_mat_result.to_frame()
```

Runnable demo: [`../examples/predmat_workflow_demo.py`](../examples/predmat_workflow_demo.py)

### 3) Run Activity Analysis from Existing `predMat`

```python
from phospy.activities import KinaseActivityAnalyzer

analyzer = KinaseActivityAnalyzer()
pred_mat = analyzer.load_pred_mat("predMat.csv")
activity = analyzer.run(pred_mat=pred_mat, phospho_matrix=phospho_matrix)
```

### 4) Build Signalomes

```python
from phospy.api import SignalomeWorkflow

signalome = SignalomeWorkflow().run(
    scoring_result=predmat_workflow_result.scoring_result,
    prediction_result=predmat_workflow_result.prediction_result,
    expression_matrix=phospho_matrix,
    kinases_of_interest=["KINASE_A", "KINASE_B"],
    site_to_protein=site_to_protein,
)
```

Runnable demo: [`../examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py)

### 5) File-First Pipeline

```python
from phospy.api import DatasetLoadOptions, KinaseActivityConfig
from phospy.pipeline import PhosRPipeline
from phospy.preprocessing import CorePreprocessingConfig

pipeline = PhosRPipeline.from_files(
    total_path="total.tsv",
    phospho_path="phospho.tsv",
    pred_mat_path="predMat.csv",
    dataset_options=DatasetLoadOptions(phospho_encoding="utf-16le"),
    preprocessing_config=CorePreprocessingConfig(max_unmatched_fraction=0.1),
    activity_config=KinaseActivityConfig(min_substrates=3),
)
outputs = pipeline.run(outdir="output")
```

When `outdir` is provided, PhosPy writes core outputs, optional kinase outputs, and `run_manifest.json`.

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

## Configuration Reference

### `DatasetLoadOptions` (`phospy.api`)

```python
DatasetLoadOptions(
    phospho_encoding: str | None = None,
    schema: DatasetSchema = DatasetSchema(),
    comparisons: tuple[tuple[str, str], ...] | None = None,
)
```

- use in `SimpleKinaseWorkflow.run(..., dataset_options=...)` and `PhosRPipeline.from_files(...)`
- `schema` customizes sample column names
- `comparisons` defines pairwise comparison groups

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

- controls filtering, protein-correction strictness, and site-matrix collapse behavior
- `max_unmatched_fraction=0.0` is strict (no silent row loss in protein correction)

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

Used by `PredMatWorkflow`, `KinaseWorkflow`, and `SimpleKinaseWorkflow`.

### `KinaseActivityConfig` (`phospy.api`)

```python
KinaseActivityConfig(
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
)
```

Used by `SimpleKinaseWorkflow` and `PhosRPipeline`.

### `SignalomeRunConfig` (`phospy.api`)

```python
SignalomeRunConfig(
    kinase_network_threshold: float = 0.9,
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

### `PredMatWorkflow`

```python
PredMatWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = "default",
)
```

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

```python
KinaseWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = "default",
)
```

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

### `KinaseActivityAnalyzer`

```python
analyzer.load_pred_mat(pred_mat_path: str | Path) -> pd.DataFrame
```

```python
analyzer.run(
    pred_mat: pd.DataFrame | PredMatResult,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult
```

```python
analyzer.write_outputs(result: KinaseActivityResult, outdir: str | Path) -> None
```

### `PhosRPipeline`

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

```python
pipeline.run(outdir: str | Path | None = None) -> CoreOutputs
```

## Dataset and Preprocessing Reference

### `DatasetSchema`

```python
DatasetSchema(
    total_cols: tuple[str, ...] = ("group1", ..., "group6"),
    phospho_cols: tuple[str, ...] = ("p_group1", ..., "p_group6"),
    corrected_cols: tuple[str, ...] = ("phospho_corrected_1", ..., "phospho_corrected_6"),
)
```

### `PhosphoDataset`

```python
PhosphoDataset(
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    *,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
)
```

```python
PhosphoDataset.from_files(
    total_path: str | Path,
    phospho_path: str | Path,
    phospho_encoding: str | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
    schema: DatasetSchema | None = None,
) -> PhosphoDataset
```

Useful accessors:

- `total_df_copy`, `phospho_df_copy`
- `total_df_live`, `phospho_df_live`
- `copy_inputs()`

`dataset.preprocessing.run(...)` requires `config=CorePreprocessingConfig(...)`.

### Analysis-Ready Helpers

```python
dataset.preprocessing.to_analysis_ready(
    result: CoreProcessingResult,
    *,
    source: str = "dataset preprocessing",
) -> AnalysisReadyPhosphoDataset
```

```python
dataset.preprocessing.run_analysis_ready(
    *,
    config: CorePreprocessingConfig,
    source: str = "dataset preprocessing",
) -> AnalysisReadyPhosphoDataset
```

```python
AnalysisReadyPhosphoDataset.from_core_processing_result(
    result: CoreProcessingResult,
    *,
    schema: DatasetSchema,
    comparisons: Sequence[tuple[str, str]] | None = None,
    source: str = "core preprocessing",
) -> AnalysisReadyPhosphoDataset
```

```python
from phospy.preprocessing.analysis_ready import build_analysis_ready_dataset
```

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

### Core Output Writer

```python
CoreOutputWriter().write(
    result: CoreProcessingResult,
    outdir: str | Path,
    *,
    format: Literal["csv", "tsv", "parquet"] = "csv",
) -> None
```

Core output basenames:

- `df_total_unique`
- `df_total_filtered`
- `df_phospho_filtered`
- `df_phospho_corrected`
- `phosr_input`
- `mat_phospho_corrected`
- `site_sequences`

## References, Results, and Exports

### Reference Resolution

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

- `substrate_map` and `motif_sequences` must both be non-empty
- kinase keys must match exactly across both mappings

```python
class ReferenceProvider(Protocol):
    def resolve(self, *, species: str, reference: str = "auto") -> ReferenceBundle: ...
```

```python
BundledReferenceProvider().resolve(
    *,
    species: str,
    reference: str = "auto",
) -> ReferenceBundle
```

### Result Objects

- `PredMatResult`: `to_frame(copy=True)`, `to_csv(path, index_label="phosphosite")`
- `SimpleKinaseWorkflowResult`: exposes `pred_mat_result`, `scoring_result`, `prediction_result`
- `SignalomeResult`: `modules`, `assignments`, `network`, `expanded_signalomes`, `to_map_data()`, `to_network_data()`, `to_csv(...)`
- `SignalomeMapData`: `modules()`, `sites()`, `kinases()`, `links()`, `to_frames()`, `to_csv(...)`
- `SignalomeNetworkData`: `nodes()`, `edges()`, `adjacency()`, `to_frames()`, `to_csv(...)`

## CLI Reference

Use `phospy --help` for full details.

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

- `RequestValidationError`: invalid user/public inputs
- `NoCandidateKinasesError`: thresholds removed all candidate kinases
- `InputCompatibilityError`: inputs are valid individually but incompatible together

## Related Docs

- [`validation.md`](validation.md)
- [`parity.md`](parity.md)
- [`../examples/predmat_workflow_demo.py`](../examples/predmat_workflow_demo.py)
- [`../examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py)
