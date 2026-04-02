# Public API Reference

This page documents the supported Python API.

PhosPy does not expose HTTP endpoints. The supported public surface is the Python API described here together with the
small `phospy` CLI described in the README.

The signatures below match the current code.

## Stable Root Imports

```python
from phospy import (
    CoreOutputWriter,
    CoreOutputs,
    CoreProcessingResult,
    DatasetPreprocessing,
    DatasetSchema,
    DatasetSiteMatrix,
    KinaseActivityAnalyzer,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseWorkflow,
    KinaseWorkflowResult,
    PhosphoDataset,
    PhosRPipeline,
    SiteMatrixResult,
)
```

## Pick the Right Entry Point

- Use `PhosphoDataset` for validated in-memory inputs and core preprocessing.
- Use `CoreOutputWriter` to write a `CoreProcessingResult` to disk.
- Use `KinaseActivityAnalyzer` when you already have a phosphosite matrix and a `predMat`.
- Use `PhosRPipeline` when you want file loading, preprocessing, optional kinase analysis, and output publishing in one place.
- Use `KinaseWorkflow` for the native end-to-end scoring and prediction workflow.

## `DatasetSchema`

Use `DatasetSchema` when your sample/value columns do not use the default group names.

### Signature

```python
DatasetSchema(
    total_cols=("group1", "group2", "group3", "group4", "group5", "group6"),
    phospho_cols=("p_group1", "p_group2", "p_group3", "p_group4", "p_group5", "p_group6"),
    corrected_cols=(
        "phospho_corrected_1",
        "phospho_corrected_2",
        "phospho_corrected_3",
        "phospho_corrected_4",
        "phospho_corrected_5",
        "phospho_corrected_6",
    ),
)
```

### What it Does

Defines the aligned numeric column groups used by dataset validation, protein correction, and site-matrix building.

### Parameters and Validation

- `total_cols`: total-proteome sample columns. Must align one-for-one with `phospho_cols` and `corrected_cols`.
- `phospho_cols`: phosphoproteome sample columns. Must align one-for-one with `total_cols` and `corrected_cols`.
- `corrected_cols`: corrected phosphosite output columns. Must align one-for-one with `total_cols` and `phospho_cols`.

Validation rules:

- all three groups must have the same length
- the aligned column groups must be compatible with each other
- structural columns such as `genes`, `gene_names`, `gene_p_site`, `localization_prob`, and `centralized_sequence` are fixed package columns and are not configured here

### Example

```python
from phospy import DatasetSchema

schema = DatasetSchema(
    total_cols=("sample_a", "sample_b"),
    phospho_cols=("p_sample_a", "p_sample_b"),
    corrected_cols=("corrected_a", "corrected_b"),
)
```

## `PhosphoDataset`

Use `PhosphoDataset` when you want one validated in-memory dataset snapshot and the bound preprocessing helpers that go with it.

### `PhosphoDataset(total_df, phospho_df, *, schema=None, comparisons=None)`

```python
PhosphoDataset(
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    *,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
)
```

Creates a dataset from in-memory pandas data frames.

#### Parameters and Validation

- `total_df`: total-proteome input table.
  - must be a pandas `DataFrame`
  - must include `genes` plus the schema's `total_cols`
  - `genes` must not contain null values
  - sample columns are coerced to numeric and fail if non-numeric values remain
- `phospho_df`: phosphoproteome input table.
  - must be a pandas `DataFrame`
  - must include `uid`, `gene_names`, `gene_p_site`, `localization_prob`, `centralized_sequence`, and the schema's `phospho_cols`
  - `uid`, `gene_names`, and `gene_p_site` must not contain null values
  - `localization_prob` and phospho sample columns are coerced to numeric and fail if non-numeric values remain
  - `localization_prob` must stay in `[0, 1]`
  - `gene_p_site` must split cleanly into gene and site parts such as `BTK_Y551`
- `schema`: optional `DatasetSchema`. Defaults to the package standard six-group schema.
- `comparisons`: optional pairwise comparisons such as `[("group1", "group4")]`.
  - each group name must exist in the schema
  - duplicate comparison pairs are rejected

### `PhosphoDataset.from_files(total_path, phospho_path, phospho_encoding=None, comparisons=None, schema=None)`

```python
PhosphoDataset.from_files(
    total_path: str | Path,
    phospho_path: str | Path,
    phospho_encoding: str | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
    schema: DatasetSchema | None = None,
) -> PhosphoDataset
```

Loads the input tables from disk and validates them.

#### Parameters and Validation

- `total_path`: path to the total table.
  - must point to a tab-delimited text file
- `phospho_path`: path to the phospho table.
  - must point to a tab-delimited text file
- `phospho_encoding`: optional file encoding for the phospho table. Defaults to UTF-8 when omitted.
- `comparisons`: same validation as the in-memory constructor.
- `schema`: same validation as the in-memory constructor.

File-loading rules:

- total and phospho files are read as TSV
- file-loaded headers are cleaned to lowercase snake case before validation
- duplicate raw headers that collapse to the same cleaned name are rejected
- the same table rules used by the in-memory constructor then apply

### `PhosphoDataset.from_validated_inputs(validated_inputs, *, comparisons=None)`

```python
PhosphoDataset.from_validated_inputs(
    validated_inputs: ValidatedCoreInputs,
    *,
    comparisons: Sequence[tuple[str, str]] | None = None,
) -> PhosphoDataset
```

Builds a dataset from already validated loader output.

#### Parameters and Validation

- `validated_inputs`: validated loader output from `DatasetLoader.validate(...)`
- `comparisons`: validated against the schema carried by `validated_inputs`

### `dataset.copy_inputs()`

```python
dataset.copy_inputs() -> tuple[pd.DataFrame, pd.DataFrame]
```

Returns deep copies of the validated total and phospho input tables.

Use this when you want caller-owned mutable copies.

### `dataset.preprocessing`

```python
dataset.preprocessing -> DatasetPreprocessing
```

Returns the bound preprocessing facade for the dataset.

This is the preferred public entry point for core preprocessing.

### `dataset.site_matrix`

```python
dataset.site_matrix -> DatasetSiteMatrix
```

Returns the bound site-matrix facade for the dataset.

Use this when you already have a corrected phosphosite table and want to build the site matrix directly.

## `DatasetPreprocessing`

`DatasetPreprocessing` is the bound preprocessing facade returned by `dataset.preprocessing`. It is part of the public API, but you normally access it through `PhosphoDataset` rather than constructing it yourself.

### `dataset.preprocessing.run(...)`

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

Runs the full core preprocessing flow: total cleanup, phosphosite filtering, protein correction, optional pairwise comparisons, and site-matrix construction.

#### Parameters and Validation

- `localization_threshold`:
  - finite numeric value in `[0, 1]`
  - phosphosite rows are retained when `localization_prob >= localization_threshold`
- `min_observed`:
  - integer with value `>= 1`
  - used in the total and phospho preprocessing steps
- `max_unmatched_fraction`:
  - finite numeric value in `[0, 1]`
  - `0.0` is strict and rejects any silent phosphosite row loss during protein correction
- `total_sentinel`:
  - numeric value used as the total-table sentinel before missing-value handling
- `phospho_sentinel`:
  - numeric value used as the phospho-table sentinel before missing-value handling
- `config`:
  - optional `CorePreprocessingConfig`
  - when supplied, it is used instead of the individual scalar arguments above

#### Return Value

Returns `CoreProcessingResult` with:

- `total_unique`
- `total_filtered`
- `phospho_filtered`
- `phospho_corrected`
- `site_matrix`

#### Example

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

matrix = core.site_matrix.matrix
corrected = core.phospho_corrected
```

## `DatasetSiteMatrix`

`DatasetSiteMatrix` is the bound site-matrix facade returned by `dataset.site_matrix`.

### `dataset.site_matrix.build(corrected_df, *, gene_p_site_col="gene_p_site", sequence_col="centralized_sequence")`

```python
dataset.site_matrix.build(
    corrected_df: pd.DataFrame,
    *,
    gene_p_site_col: str = "gene_p_site",
    sequence_col: str = "centralized_sequence",
) -> SiteMatrixResult
```

Builds a PhosR-style phosphosite matrix from corrected phosphosite rows.

#### Parameters and Validation

- `corrected_df`:
  - must be a pandas `DataFrame`
  - must include `gene_p_site_col`, `sequence_col`, and the dataset schema's `corrected_cols`
  - `gene_p_site_col` must not contain null values
  - the corrected value columns must be numeric
  - `gene_p_site_col` must split cleanly into gene and site parts
- `gene_p_site_col`: source phosphosite identifier column. Defaults to `gene_p_site`.
- `sequence_col`: source central sequence column. Defaults to `centralized_sequence`.

#### Return Value

Returns `SiteMatrixResult` with:

- `phosr_input`: intermediate table used to build the matrix
- `matrix`: phosphosite-by-sample numeric matrix
- `sequences`: centralized sequences indexed by phosphosite ID
- `row_drop_stats`: retained and dropped row counts

#### Example

```python
core = dataset.preprocessing.run()
site_matrix = dataset.site_matrix.build(core.phospho_corrected)
```

## `CoreOutputWriter`

Use `CoreOutputWriter` to write a `CoreProcessingResult` to disk.

### `CoreOutputWriter.write(result, outdir, *, format="csv")`

```python
CoreOutputWriter.write(
    result: CoreProcessingResult,
    outdir: str | Path,
    *,
    format: str = "csv",
) -> None
```

Writes the core preprocessing outputs to a directory.

#### Parameters and Validation

- `result`: a `CoreProcessingResult`
- `outdir`: output directory path
- `format`:
  - supported values are `"csv"`, `"tsv"`, and `"parquet"`
  - parquet output requires an installed pandas parquet engine such as `pyarrow`

#### Example

```python
from phospy import CoreOutputWriter

writer = CoreOutputWriter()
writer.write(core, outdir="results", format="parquet")
```

## Standalone Preprocessing Helpers

These are public helpers for targeted preprocessing work. They are useful when you do not want the full dataset-bound preprocessing path.

### `filter_localized_sites(df, *, localization_col="localization_prob", threshold=0.75, return_summary=False)`

```python
filter_localized_sites(
    df: pd.DataFrame,
    *,
    localization_col: str = "localization_prob",
    threshold: float = 0.75,
    return_summary: bool = False,
) -> pd.DataFrame | LocalizationFilterResult
```

Keeps rows where the localisation score is greater than or equal to `threshold`.

#### Parameters and Validation

- `df`: must be a pandas `DataFrame`
- `localization_col`:
  - column must exist in `df`
  - column values must be numeric
- `threshold`:
  - finite numeric value in `[0, 1]`
- `return_summary`:
  - `False` returns only the filtered frame
  - `True` returns `LocalizationFilterResult(filtered=..., summary=...)`

#### Example

```python
from phospy.preprocessing import filter_localized_sites

result = filter_localized_sites(phospho_df, threshold=0.75, return_summary=True)
filtered = result.filtered
removed_rows = result.summary.removed_rows
```

### `filter_sites_by_coverage(df, *, columns, min_coverage=0.0, return_summary=False)`

```python
filter_sites_by_coverage(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    min_coverage: float = 0.0,
    return_summary: bool = False,
) -> pd.DataFrame | CoverageFilterResult
```

Keeps rows whose observed-value proportion across the selected columns is greater than or equal to `min_coverage`.

#### Parameters and Validation

- `df`: must be a pandas `DataFrame`
- `columns`:
  - must contain at least one column name
  - every listed column must exist in `df`
  - selected columns must be numeric
- `min_coverage`:
  - finite numeric value in `[0, 1]`
- `return_summary`:
  - `False` returns only the filtered frame
  - `True` returns `CoverageFilterResult(filtered=..., summary=...)`

#### Example

```python
from phospy.preprocessing import filter_sites_by_coverage

result = filter_sites_by_coverage(
    phospho_df,
    columns=["p_group1", "p_group2", "p_group3", "p_group4", "p_group5", "p_group6"],
    min_coverage=0.5,
    return_summary=True,
)
required_count = result.summary.required_observed_count
```

## `KinaseActivityAnalyzer`

Use `KinaseActivityAnalyzer` when you already have a phosphosite matrix and a `predMat` and want downstream kinase summaries.

### `load_pred_mat(pred_mat_path)`

```python
KinaseActivityAnalyzer.load_pred_mat(pred_mat_path: str | Path) -> pd.DataFrame
```

Loads and validates a prediction matrix from disk.

#### Parameters and Validation

- `pred_mat_path`:
  - must point to a CSV file
  - the first column is used as the phosphosite index
  - kinase score columns must be numeric
  - scores must stay in `[0, 1]`
  - the matrix must contain at least one row and at least one kinase column
  - the phosphosite index must be unique and non-null

### `analyze(pred_mat, phospho_matrix, threshold=0.6, min_substrates=3, top_n_substrates=20)`

```python
KinaseActivityAnalyzer.analyze(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult
```

Computes weighted activity scores, KSEA-style summaries, and kinase target summaries from in-memory inputs.

#### Parameters and Validation

- `pred_mat`:
  - must satisfy the public `predMat` schema
  - must have a unique, non-null phosphosite index
  - scores must be numeric and in `[0, 1]`
- `phospho_matrix`:
  - must be a numeric phosphosite matrix
  - must have a unique, non-null phosphosite index
  - must contain at least one row and at least one numeric column
- `threshold`:
  - finite numeric value in `[0, 1]`
- `min_substrates`:
  - integer with value `>= 1`
- `top_n_substrates`:
  - integer with value `>= 1`

Compatibility validation:

- `pred_mat` and `phospho_matrix` must overlap by at least one phosphosite row
- that overlap must cover at least 10% of the phosphosite matrix

### `load_and_analyze(pred_mat_path, phospho_matrix, threshold=0.6, min_substrates=3, top_n_substrates=20)`

```python
analyzer.load_and_analyze(
    pred_mat_path: str | Path,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult
```

Convenience method that combines `load_pred_mat(...)` and `analyze(...)`.

#### Example

```python
from phospy import KinaseActivityAnalyzer

analyzer = KinaseActivityAnalyzer()
result = analyzer.load_and_analyze(
    pred_mat_path="predMat.csv",
    phospho_matrix=core.site_matrix.matrix,
    threshold=0.6,
    min_substrates=1,
    top_n_substrates=1,
)
```

### `write_outputs(result, outdir)`

```python
analyzer.write_outputs(result: KinaseActivityResult, outdir: str | Path) -> None
```

Writes the downstream kinase-analysis tables to disk.

#### Parameters

- `result`: `KinaseActivityResult`
- `outdir`: output directory path

#### Return Example

`KinaseActivityResult` exposes:

- `weighted_activity`
- `ksea_scores`
- `ksea_counts`
- `target_counts`
- `target_table`

```python
target_counts = result.target_counts
ksea_scores = result.ksea_scores
```

## `PhosRPipeline`

Use `PhosRPipeline` when you want file loading, preprocessing, optional `predMat` analysis, and bundle publishing in one place.

### `PhosRPipeline.from_files(...)`

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
) -> PhosRPipeline
```

Builds a pipeline from file paths.

#### Parameters and Validation

- `total_path`, `phospho_path`, `phospho_encoding`, `comparisons`, and `schema` follow the same rules as `PhosphoDataset.from_files(...)`
- `pred_mat_path` is optional
  - when supplied, it must point to a CSV `predMat`
  - it is validated when downstream kinase analysis runs
- `localization_threshold`, `min_observed`, `max_unmatched_fraction`, `total_sentinel`, and `phospho_sentinel` follow the same rules as `dataset.preprocessing.run(...)`

### `PhosRPipeline.from_request(request)`

```python
PhosRPipeline.from_request(request: CorePipelineRequest) -> PhosRPipeline
```

Builds a pipeline from an already validated request object.

### `pipeline.run(outdir=None)`

```python
pipeline.run(outdir: str | Path | None = None) -> CoreOutputs
```

Runs the pipeline and optionally publishes files.

#### Parameters and Behaviour

- `outdir`:
  - when provided, PhosPy writes the output bundle to disk
  - when `None`, PhosPy returns in-memory results only

If `pred_mat_path` was supplied, the run includes downstream kinase analysis. Otherwise it returns only the core outputs.

#### Return Value

Returns `CoreOutputs` with:

- `core`: `CoreProcessingResult`
- `kinase_activity`: `KinaseActivityResult | None`

When `outdir` is provided, the published bundle includes the core tables, any downstream kinase-analysis tables, and `run_manifest.json`.

#### Example

```python
from phospy import PhosRPipeline

pipeline = PhosRPipeline.from_files(
    total_path="total.tsv",
    phospho_path="phospho.tsv",
    pred_mat_path="predMat.csv",
    max_unmatched_fraction=0.1,
)
outputs = pipeline.run(outdir="output")
```

## `KinaseWorkflow`

Use `KinaseWorkflow` for the native Python end-to-end scoring and prediction workflow.

### `KinaseWorkflow(flank_size=7, kernel="rbf", svm_mode="default")`

```python
KinaseWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: str = "default",
)
```

Creates a workflow runner.

#### Parameters and Validation

- `flank_size`: flank size used by motif scoring
- `kernel`: SVM kernel name passed to the prediction layer
- `svm_mode`:
  - supported public values are `"default"` and `"r_parity"`
  - `"default"` is the normal native mode
  - `"r_parity"` is a narrower learner-seam comparison mode used in PhosR parity work

### `workflow.run(...)`

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

Runs the native workflow end to end.

#### Parameters and Validation

- `phospho_matrix`:
  - must be a numeric phosphosite matrix
  - must contain at least one row and at least one numeric column
  - must have a unique, non-null phosphosite index
- `substrate_map`:
  - must not be empty
  - must overlap with the phosphosite IDs in `phospho_matrix`
- `site_sequences`:
  - optional
  - when supplied, must be a mapping keyed by phosphosite ID or a pandas `Series` indexed by phosphosite ID
  - when supplied, it must contain an entry for every phosphosite row in `phospho_matrix`
- `motif_sequences`:
  - optional only when `allow_profile_only_fallback=True`
  - if provided, it must not be empty
  - if provided, `site_sequences` is also required
- `min_substrates`, `min_motif_size`, `ensemble_size`, `top`, `inclusion`, `n_iterations`:
  - integers with value `>= 1`
- `allow_profile_only_fallback`:
  - set to `True` to allow prediction without motif data
- `score_threshold`:
  - finite numeric value in `[0, 1]`
- `random_state`: optional integer seed
- `svm_mode`:
  - optional per-run override for the workflow's default SVM mode

### `workflow.run_request(request)`

```python
workflow.run_request(request: KinaseWorkflowRequest) -> KinaseWorkflowResult
```

Runs the workflow from an already validated request object.

#### Return Value

Returns `KinaseWorkflowResult` with:

- `profile_result`
- `motif_result`
- `scoring_result`
- `prediction_result`

`prediction_result` is a `KinasePredictionResult`. The most commonly used fields are:

- `pred_matrix`
- `substrate_list`

#### Example

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

## Result Objects

The main result objects exposed at the package root are small dataclasses. The attributes below are the ones you are most likely to use.

### `CoreProcessingResult`

```python
CoreProcessingResult(
    total_unique: pd.DataFrame,
    total_filtered: pd.DataFrame,
    phospho_filtered: pd.DataFrame,
    phospho_corrected: pd.DataFrame,
    site_matrix: SiteMatrixResult,
)
```

### `SiteMatrixResult`

```python
SiteMatrixResult(
    phosr_input: pd.DataFrame,
    matrix: pd.DataFrame,
    sequences: pd.Series,
    row_drop_stats: object,
)
```

### `CoreOutputs`

```python
CoreOutputs(
    core: CoreProcessingResult,
    kinase_activity: KinaseActivityResult | None,
)
```

### `KinaseActivityResult`

```python
KinaseActivityResult(
    weighted_activity: pd.DataFrame,
    ksea_scores: pd.DataFrame,
    ksea_counts: pd.DataFrame,
    target_counts: pd.Series,
    target_table: pd.DataFrame,
)
```

### `KinasePredictionResult`

The most commonly used fields are `pred_matrix` and `substrate_list`.

### `KinaseWorkflowResult`

```python
KinaseWorkflowResult(
    profile_result: object,
    motif_result: object | None,
    scoring_result: object,
    prediction_result: KinasePredictionResult,
)
```

## Validation Summary

These checks come up most often across the public API:

- total and phospho file inputs are TSV; `predMat` file input is CSV with the first column used as the index
- file-loaded total and phospho headers are cleaned to lowercase snake case before validation
- duplicate cleaned column names are rejected
- `gene_p_site` must split into gene and site parts such as `BTK_Y551`
- `localization_prob` and `predMat` scores must stay in `[0, 1]`
- by default, protein correction allows no silent row loss
- downstream kinase analysis requires overlapping phosphosite IDs between `predMat` and the phosphosite matrix
- native workflow inputs must share phosphosite IDs across the matrix, substrate map, and sequence inputs

For the short guide to validation rules and the PhosR parity scope, see [`validation-and-parity.md`](validation-and-parity.md).
