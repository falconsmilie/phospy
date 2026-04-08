# Public API Reference

PhosPy does not expose HTTP endpoints. The supported public surface is the Python API below plus the `phospy` CLI.

## Root Imports

```python
from phospy import (
    KinaseActivityAnalyzer,
    KinaseWorkflow,
    PhosphoDataset,
    PhosRPipeline,
    PredMatResult,
    PredMatWorkflow,
)
```

Use:

- `PhosphoDataset` for validated inputs and core preprocessing
- `KinaseActivityAnalyzer` for analysis from an existing `predMat`
- `PhosRPipeline` for file loading plus publishing
- `PredMatWorkflow` for one obvious public predMat generation path
- `PredMatResult` for the canonical in-memory and export contract of a generated `predMat`
- `KinaseWorkflow` for the native end-to-end prediction flow

## Common File Rules

- total input: TSV
- phospho input: TSV
- `predMat`: CSV, with the first column used as the phosphosite index
- file-loaded total and phospho headers are normalised to lowercase snake case
- duplicate cleaned headers are rejected

## `DatasetSchema`

Use `DatasetSchema` when your numeric sample columns do not use the default names.

```python
from phospy.dataset_schema import DatasetSchema

schema = DatasetSchema(
    total_cols=("sample_a", "sample_b"),
    phospho_cols=("p_sample_a", "p_sample_b"),
    corrected_cols=("corrected_a", "corrected_b"),
)
```

Rules:

- all three column groups must have the same length
- the groups must align one-for-one
- structural columns such as `genes`, `gene_names`, `gene_p_site`, `localization_prob`, and `centralized_sequence` are fixed package columns

## `PhosphoDataset`

Use `PhosphoDataset` when you want one validated dataset owner and the bound helpers that go with it.

### Constructor

```python
PhosphoDataset(
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    *,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
)
```

Validation highlights:

- `total_df` must include `genes` and the schema's `total_cols`
- `phospho_df` must include `uid`, `gene_names`, `gene_p_site`, `localization_prob`, `centralized_sequence`, and the schema's `phospho_cols`
- required identifier columns must not contain null values
- numeric sample columns are coerced to numeric and fail if non-numeric values remain
- `localization_prob` must stay in `[0, 1]`
- `gene_p_site` must split cleanly into gene and site parts such as `BTK_Y551`
- comparison pairs must use known schema group names and must not be duplicated

### `from_files(...)`

```python
PhosphoDataset.from_files(
    total_path: str | Path,
    phospho_path: str | Path,
    phospho_encoding: str | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
    schema: DatasetSchema | None = None,
) -> PhosphoDataset
```

Use this for the normal file-based workflow.

### Data Access

```python
dataset.total_df_copy -> pd.DataFrame
dataset.phospho_df_copy -> pd.DataFrame
dataset.total_df_live -> pd.DataFrame
dataset.phospho_df_live -> pd.DataFrame
dataset.copy_inputs() -> tuple[pd.DataFrame, pd.DataFrame]
```

- `*_copy` and `copy_inputs()` return detached copies
- `*_live` returns the dataset-owned DataFrames directly

### Preprocessing

```python
dataset.preprocessing.run(
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float | int = 10.0,
    phospho_sentinel: float | int = 12.0,
    config: CorePreprocessingConfig | None = None,
) -> CoreProcessingResult
```

This runs the full core preprocessing path.

Rules:

- `localization_threshold` and `max_unmatched_fraction` must be in `[0, 1]`
- `min_observed` must be `>= 1`
- `config` can replace the scalar arguments, but cannot be combined with them when you override values
- `max_unmatched_fraction=0.0` is the strict mode and rejects silent row loss during protein correction

Returns `CoreProcessingResult` with:

- `total_unique`
- `total_filtered`
- `phospho_filtered`
- `phospho_corrected`
- `site_matrix`

Example:

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

matrix = core.site_matrix.matrix
corrected = core.phospho_corrected
```

### Site Matrix Builder

```python
dataset.site_matrix.build(
    corrected_df: pd.DataFrame,
    *,
    gene_p_site_col: str = "gene_p_site",
    sequence_col: str = "centralized_sequence",
) -> SiteMatrixResult
```

Use this when you already have a corrected phosphosite table.

Rules:

- `corrected_df` must include `gene_p_site_col`, `sequence_col`, and the schema's `corrected_cols`
- corrected value columns must be numeric
- `gene_p_site_col` must split cleanly into gene and site parts

Returns `SiteMatrixResult` with:

- `phosr_input`
- `matrix`
- `sequences`
- `row_drop_stats`

## `CorePreprocessingConfig`

`CorePreprocessingConfig` lets you pass one preprocessing object instead of individual scalar options.

```python
from phospy.core_processing import CorePreprocessingConfig

config = CorePreprocessingConfig(
    localization_threshold=0.75,
    min_observed=4,
    total_sentinel=10.0,
    phospho_sentinel=12.0,
    max_unmatched_fraction=0.1,
)
```

You can use it with `dataset.preprocessing.run(config=...)` or `PhosRPipeline(..., preprocessing_config=...)`.

## `CoreOutputWriter`

```python
from phospy.writers import CoreOutputWriter

CoreOutputWriter().write(
    result=core,
    outdir="results",
    format="csv",
)
```

Supported formats:

- `csv`
- `tsv`
- `parquet` (requires a parquet engine such as `pyarrow`)

Core outputs written to disk are:

- `df_total_unique`
- `df_total_filtered`
- `df_phospho_filtered`
- `df_phospho_corrected`
- `phosr_input`
- `mat_phospho_corrected`
- `site_sequences`

## Targeted Preprocessing Helpers

These helpers are public, but most users should start with `PhosphoDataset.preprocessing.run()`.

### `filter_localized_sites(...)`

```python
filter_localized_sites(
    df: pd.DataFrame,
    *,
    localization_col: str = "localization_prob",
    threshold: float = 0.75,
    return_summary: bool = False,
)
```

Keeps rows where `localization_col >= threshold`.

### `filter_sites_by_coverage(...)`

```python
filter_sites_by_coverage(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    min_coverage: float = 0.0,
    return_summary: bool = False,
)
```

Keeps rows whose observed-value fraction across `columns` is at least `min_coverage`.

## `KinaseActivityAnalyzer`

Use `KinaseActivityAnalyzer` when you already have a phosphosite matrix and a `predMat`.

### `load_pred_mat(...)`

```python
analyzer.load_pred_mat(pred_mat_path: str | Path) -> pd.DataFrame
```

Rules:

- the file must be a CSV
- the first column becomes the phosphosite index
- kinase score columns must be numeric
- scores must stay in `[0, 1]`
- the index must be unique and non-null

### `run(...)`

```python
analyzer.run(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult
```

Rules:

- `pred_mat` may be a valid `predMat` `DataFrame` or a `PredMatResult`
- the public validation boundary normalizes `PredMatResult` to the internal validated `DataFrame` contract before analysis runs
- `phospho_matrix` must be a numeric phosphosite matrix with a unique, non-null index
- `threshold` must be in `[0, 1]`
- `min_substrates` and `top_n_substrates` must be `>= 1`
- `pred_mat` and `phospho_matrix` must overlap by at least one phosphosite row
- that overlap must cover at least 10% of the phosphosite matrix

Example:

```python
from phospy import KinaseActivityAnalyzer

analyzer = KinaseActivityAnalyzer()
result = analyzer.run(
    pred_mat=analyzer.load_pred_mat("predMat.csv"),
    phospho_matrix=core.site_matrix.matrix,
    threshold=0.6,
    min_substrates=1,
    top_n_substrates=1,
)
```

### `write_outputs(...)`

```python
analyzer.write_outputs(result: KinaseActivityResult, outdir: str | Path) -> None
```

This writes:

- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `kinase_target_table.csv`

`KinaseActivityResult` exposes:

- `weighted_activity`
- `ksea_scores`
- `ksea_counts`
- `target_counts`
- `target_table`

## `PhosRPipeline`

Use `PhosRPipeline` when you want preprocessing and optional kinase analysis in one place.

### Constructor

```python
PhosRPipeline(
    dataset: PhosphoDataset,
    pred_mat: pd.DataFrame | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> PhosRPipeline
```

Use this when your dataset is already in memory.

Rules:

- `dataset` must be a `PhosphoDataset`
- `pred_mat`, when provided, may be a valid in-memory `predMat` `DataFrame` or a `PredMatResult`
- pipeline construction normalizes `PredMatResult` at the validation boundary; downstream pipeline execution only sees the validated `DataFrame`
- `preprocessing_config` cannot be mixed with overridden scalar preprocessing options

### `from_files(...)`

```python
PhosRPipeline.from_files(
    total_path: str | Path,
    phospho_path: str | Path,
    pred_mat_path: str | Path | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
    phospho_encoding: str | None = None,
    schema: DatasetSchema | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> PhosRPipeline
```

Notes:

- `total_path` and `phospho_path` follow the same rules as `PhosphoDataset.from_files(...)`
- `pred_mat_path`, when provided, is loaded and schema-validated during pipeline construction
- overlap between `predMat` and the processed site matrix is checked at runtime after preprocessing

### `run(...)`

```python
pipeline.run(outdir: str | Path | None = None) -> CoreOutputs
```

- when `outdir` is set, files are published to disk
- when `outdir` is `None`, PhosPy returns only the in-memory results

`CoreOutputs` exposes:

- `core`
- `kinase_activity`

When `outdir` is set, the output bundle also includes `run_manifest.json`.

## `PredMatWorkflow`

Use `PredMatWorkflow` when your goal is to generate a `predMat` from one documented public workflow.

### Constructor

```python
PredMatWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: str = "default",
)
```

Supported public `svm_mode` values:

- `"default"`
- `"r_parity"`

### `run(...)`

```python
workflow.run(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
    min_substrates: int = 1,
    min_motif_size: int = 1,
    allow_profile_only_fallback: bool = False,
    ensemble_size: int = 10,
    top: int = 50,
    score_threshold: float = 0.8,
    inclusion: int = 20,
    n_iterations: int = 5,
    random_state: int | None = None,
    svm_mode: str | None = None,
) -> PredMatWorkflowResult
```

Validation highlights:

- `phospho_matrix` must be a numeric phosphosite matrix with a unique, non-null index
- `substrate_map` must not be empty and must overlap with `phospho_matrix`
- `site_sequences`, when supplied, must be a mapping or a pandas `Series` keyed by phosphosite ID
- `motif_sequences` is required unless `allow_profile_only_fallback=True`
- when `motif_sequences` is provided, `site_sequences` is also required
- `min_substrates`, `min_motif_size`, `ensemble_size`, `top`, `inclusion`, and `n_iterations` must be `>= 1`
- `score_threshold` must be in `[0, 1]`

The returned `PredMatWorkflowResult` exposes:

- `scoring_result`
- `prediction_result`
- `pred_mat_result`

`pred_mat_result` is the canonical predMat contract.

Row and column semantics:

- rows are phosphosite identifiers
- columns are kinase identifiers
- values are prediction scores in `[0, 1]`

In-memory access:

- `result.pred_mat_result.data_frame` returns the owned in-memory `DataFrame`
- `result.pred_mat_result.to_frame()` returns a detached copy

CSV export:

- `result.pred_mat_result.to_csv("predMat.csv")`
- the export is UTF-8 CSV with a `phosphosite` index column and deterministic newline / float formatting
- `load_pred_mat(...)` can read the exported file back as a validated `predMat`

Example:

```python
from phospy import PredMatWorkflow
from phospy.io import load_pred_mat

workflow = PredMatWorkflow(svm_mode="default")
result = workflow.run(
    phospho_matrix=phospho_matrix,
    substrate_map=substrate_map,
    site_sequences=site_sequences,
    motif_sequences=motif_sequences,
)

pred_mat_result: PredMatResult = result.pred_mat_result
pred_mat = pred_mat_result.to_frame(copy=False)
pred_mat_result.to_csv("predMat.csv")
reloaded_pred_mat = load_pred_mat("predMat.csv")
```

## `KinaseWorkflow`

Use `KinaseWorkflow` when you want the fuller native Python scoring and prediction result, including profile and motif scoring outputs.

### Constructor

```python
KinaseWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: str = "default",
)
```

Supported public `svm_mode` values:

- `"default"`
- `"r_parity"`

### `run(...)`

```python
workflow.run(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
    min_substrates: int = 1,
    min_motif_size: int = 1,
    allow_profile_only_fallback: bool = False,
    ensemble_size: int = 10,
    top: int = 50,
    score_threshold: float = 0.8,
    inclusion: int = 20,
    n_iterations: int = 5,
    random_state: int | None = None,
    svm_mode: str | None = None,
) -> KinaseWorkflowResult
```

Validation highlights:

- `phospho_matrix` must be a numeric phosphosite matrix with a unique, non-null index
- `substrate_map` must not be empty and must overlap with `phospho_matrix`
- `site_sequences`, when supplied, must be a mapping or a pandas `Series` keyed by phosphosite ID
- `motif_sequences` is required unless `allow_profile_only_fallback=True`
- when `motif_sequences` is provided, `site_sequences` is also required
- `min_substrates`, `min_motif_size`, `ensemble_size`, `top`, `inclusion`, and `n_iterations` must be `>= 1`
- `score_threshold` must be in `[0, 1]`

The returned `KinaseWorkflowResult` exposes:

- `profile_result`
- `motif_result`
- `scoring_result`
- `prediction_result`

The most commonly used fields on `prediction_result` are:

- `pred_matrix`
- `substrate_list`

For the canonical predMat contract, prefer `PredMatWorkflowResult.pred_mat_result` or `KinasePredictionResult.pred_mat_result` over direct DataFrame access.

Example:

```python
from phospy import KinaseWorkflow

workflow = KinaseWorkflow(svm_mode="default")
result = workflow.run(
    phospho_matrix=phospho_matrix,
    substrate_map=substrate_map,
    site_sequences=site_sequences,
    motif_sequences=motif_sequences,
)

pred_matrix = result.prediction_result.pred_matrix
```

## CLI

Use `phospy --help` for the full help text.

Main options:

- `--total` required TSV path
- `--phospho` required TSV path
- `--outdir` required output directory
- `--pred-mat` optional `predMat` CSV path
- `--phospho-encoding` optional encoding for the phospho table
- `--localization-threshold` default `0.75`
- `--min-observed` default `4`
- `--total-sentinel` default `10.0`
- `--phospho-sentinel` default `12.0`
- `--max-unmatched-fraction` default `0.0`
- `--kinase-activity-threshold` default `0.6`
- `--kinase-activity-min-substrates` default `3`
- `--kinase-activity-top-n-substrates` default `20`

Example:

```bash
phospy \
  --total total.tsv \
  --phospho phospho.tsv \
  --pred-mat predMat.csv \
  --outdir output
```

## Result Objects

The main public workflows return small dataclasses that hold output tables. Treat them as standalone results.

```python
CoreProcessingResult(
    total_unique: pd.DataFrame,
    total_filtered: pd.DataFrame,
    phospho_filtered: pd.DataFrame,
    phospho_corrected: pd.DataFrame,
    site_matrix: SiteMatrixResult,
)

SiteMatrixResult(
    phosr_input: pd.DataFrame,
    matrix: pd.DataFrame,
    sequences: pd.Series,
    row_drop_stats: dict[str, int],
)

CoreOutputs(
    core: CoreProcessingResult,
    kinase_activity: KinaseActivityResult | None,
)

PredMatWorkflowResult(
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult,
    pred_mat_result: PredMatResult,
)

# canonical predMat contract
PredMatWorkflowResult.pred_mat_result -> PredMatResult

PredMatResult(
    data_frame: pd.DataFrame,
)

PredMatResult.to_frame(copy: bool = True) -> pd.DataFrame
PredMatResult.to_csv(path: str | Path, index_label: str = "phosphosite") -> Path

`SiteMatrixResult` is a separate preprocessing result for the corrected phosphosite site matrix. It is not the prediction-matrix contract. The canonical prediction-matrix result is `PredMatResult`.

KinasePredictionResult(
    pred_matrix: pd.DataFrame,
    substrate_list: dict[str, list[str]],
)

# canonical predMat contract on the lower-level result
KinasePredictionResult.pred_mat_result -> PredMatResult

KinaseActivityResult(
    weighted_activity: pd.DataFrame,
    ksea_scores: pd.DataFrame,
    ksea_counts: pd.Series,
    target_counts: pd.Series,
    target_table: pd.DataFrame,
)

KinaseWorkflowResult(
    profile_result: KinaseProfileResult,
    motif_result: MotifScoringResult | None,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult,
)
```

## Validation in One Glance

These checks come up most often:

- total and phospho files are TSV; `predMat` is CSV
- `gene_p_site` must look like `BTK_Y551`
- `localization_prob` and `predMat` scores must stay in `[0, 1]`
- by default, protein correction allows no silent row loss
- kinase activity analysis requires overlap between `predMat` and the phosphosite matrix
- native workflow inputs must share phosphosite IDs across matrix, substrate map, and sequence data
