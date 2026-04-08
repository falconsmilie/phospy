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
    SignalomeMapData,
    SignalomeNetworkData,
    SignalomeResult,
    SignalomeWorkflow,
)
```

Use:

- `PhosphoDataset` for validated inputs and core preprocessing
- `KinaseActivityAnalyzer` for analysis from an existing `predMat`
- `PhosRPipeline` for file loading plus publishing
- `PredMatWorkflow` for one obvious public predMat generation path
- `PredMatResult` for the stable in-memory and export contract of a generated `predMat`
- `KinaseWorkflow` for the native end-to-end prediction flow
- `SignalomeWorkflow` for validated signalome construction from scoring and prediction outputs
- `SignalomeResult` for the stable in-memory and export contract of a constructed signalome
- `SignalomeNetworkData` for graph-friendly kinase-network tables derived from a signalome result

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

Use `PredMatWorkflow` when your goal is to generate a `predMat` from one documented public workflow. For the shortest file-backed example, see [`examples/predmat_workflow_demo.py`](../examples/predmat_workflow_demo.py).

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
- when motif-aware prediction is enabled, sequence coverage is checked only for phosphosites that participate in scoring and prediction; full-matrix sequence coverage is not required
- `min_substrates`, `min_motif_size`, `ensemble_size`, `top`, `inclusion`, and `n_iterations` must be `>= 1`
- `score_threshold` must be in `[0, 1]`

The returned `PredMatWorkflowResult` exposes:

- `scoring_result`
- `prediction_result`
- `pred_mat_result`

`pred_mat_result` is the stable predMat result object.

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
import json
from pathlib import Path

import pandas as pd

from phospy import PredMatWorkflow
from phospy.io import load_pred_mat

phospho_matrix = pd.read_csv("predmat_phospho_matrix.csv", index_col=0)
site_sequences = json.loads(Path("predmat_site_sequences.json").read_text())
substrate_map = json.loads(Path("predmat_substrate_map.json").read_text())
motif_sequences = json.loads(Path("predmat_motif_sequences.json").read_text())

workflow = PredMatWorkflow(flank_size=2)
result = workflow.run(
    phospho_matrix=phospho_matrix,
    substrate_map=substrate_map,
    site_sequences=site_sequences,
    motif_sequences=motif_sequences,
    min_substrates=2,
    min_motif_size=2,
    ensemble_size=3,
    top=4,
    score_threshold=0.75,
    inclusion=3,
    n_iterations=2,
    random_state=17,
)

pred_mat_result = result.pred_mat_result
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

For the stable predMat result, prefer `PredMatWorkflowResult.pred_mat_result` or `KinasePredictionResult.pred_mat_result` over direct DataFrame access.

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


## `SignalomeWorkflow`

Use `SignalomeWorkflow` when you already have aligned scoring and prediction outputs and want one validated step that constructs signalomes before you derive map-ready or network-ready tables.

A runnable end-to-end example lives in [`examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py).

### `run(...)`

```python
SignalomeWorkflow().run(
    *,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult | PredMatResult,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    kinase_network_threshold: float = 0.9,
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
) -> SignalomeResult
```

Validation highlights:

- `scoring_result` and `prediction_result` must align to a shared phosphosite index
- `prediction_result` can be either `KinasePredictionResult` or `PredMatResult`
- `expression_matrix` must be a numeric phosphosite matrix with a unique, non-null index
- `kinases_of_interest` must not be empty and must be present in the aligned kinase columns
- `kinase_network_threshold` and `signalome_cutoff` must be in `[0, 1]`
- `module_count`, when supplied, must be `>= 1`
- `min_kinase_module_share_percent` must be `>= 0`

The returned `SignalomeResult` gives you:

- `result.modules` for module summaries and kinase-to-module relationships
- `result.assignments` for site and protein assignments
- `result.to_map_data()` for map-ready plotting tables
- `result.to_network_data()` for graph-friendly kinase-network tables

Example:

```python
signalome_result = SignalomeWorkflow().run(
    scoring_result=pred_mat_result.scoring_result,
    prediction_result=pred_mat_result.prediction_result,
    expression_matrix=phospho_matrix,
    kinases_of_interest=["KINASE_A", "KINASE_B"],
    signalome_cutoff=0.5,
)

map_data = signalome_result.to_map_data()
network_data = signalome_result.to_network_data()
```

## `SignalomeResult`

`SignalomeWorkflow.run(...)` returns `SignalomeResult`. Use the named result views below as the stable contract for downstream work.

### Stable Access Paths

```python
result.modules
result.assignments
result.network
result.expanded_signalomes
```

Use:

- `result.modules` for module-centric outputs
- `result.assignments` for site-level and protein-level assignments
- `result.network` for graph-friendly kinase-network tables
- `result.expanded_signalomes` for kinase-of-interest views built from the module assignments
- `result.to_map_data()` for serialisable map-ready plotting data derived from the signalome result
- `result.to_network_data()` for graph-friendly kinase-network tables derived from the signalome result

The older convenience attributes remain available:

- `result.signalome_modules`
- `result.kinase_module_relationships`
- `result.site_assignments`
- `result.protein_assignments`
- `result.protein_modules`
- `result.kinase_network_nodes`
- `result.kinase_network_edges`
- `result.kinase_correlation_matrix`

### `result.modules`

```python
result.modules.to_frame(copy: bool = True) -> pd.DataFrame
result.modules.to_relationship_table(copy: bool = True) -> pd.DataFrame
```

`result.modules.to_frame()` returns the wide module-by-kinase percentage matrix.

`result.modules.to_relationship_table()` returns the long relationship table with:

- `module_id`
- `kinase`
- `share_percent`

The relationship table contains one row per non-zero kinase-to-module relationship and is sorted by `module_id`, descending `share_percent`, then `kinase`.

### `result.assignments`

```python
result.assignments.sites(copy: bool = True) -> pd.DataFrame
result.assignments.proteins(copy: bool = True) -> pd.DataFrame
```

`result.assignments.sites()` returns the site assignment table with:

- index: `site_id`
- `protein_id`
- `module_id`
- `top_kinase`
- `top_score`

`result.assignments.proteins()` returns the protein assignment table with:

- index: `protein_id`
- `module_id`
- `site_count`

### `result.network`

```python
result.network.adjacency(copy: bool = True) -> pd.DataFrame
result.network.nodes(copy: bool = True) -> pd.DataFrame
result.network.edges(copy: bool = True) -> pd.DataFrame
```

`result.network.adjacency()` returns the kinase correlation matrix used to derive network edges.

`result.network.nodes()` returns the node table with:

- index: `kinase`
- `degree`
- `n_substrates`

`result.network.edges()` returns the edge table with:

- `source_kinase`
- `target_kinase`
- `correlation`

### Export Contract

```python
result.to_frames(copy: bool = True, include_inputs: bool = False) -> dict[str, pd.DataFrame]
result.to_csv(directory: str | Path, include_inputs: bool = False) -> dict[str, Path]
```

By default, `to_frames(...)` and `to_csv(...)` return or write only the stable user-facing outputs:

- `signalome_modules`
- `kinase_module_relationships`
- `site_assignments`
- `protein_assignments`
- `kinase_network_nodes`
- `kinase_network_edges`
- `kinase_correlation_matrix`

Pass `include_inputs=True` to also include the aligned `scoring_matrix`, `pred_mat`, and `expression_matrix` that fed the signalome construction step.

The CSV export uses UTF-8, deterministic newline handling, and stable float formatting.

## `SignalomeNetworkData`

Build graph-friendly kinase-network data from a signalome result:

```python
network_data = result.to_network_data()
```

This step is intentionally separate from signalome construction. `SignalomeNetworkData` contains deterministic node, edge, and adjacency outputs that graph or plotting code can consume without binding the public API to one graph library.

### Stable Access Paths

```python
network_data.nodes(copy: bool = True) -> pd.DataFrame
network_data.edges(copy: bool = True) -> pd.DataFrame
network_data.adjacency(copy: bool = True) -> pd.DataFrame
```

`network_data.nodes()` returns the node table with:

- index: `kinase`
- `degree`
- `n_substrates`
- `module_count`
- `total_share_percent`
- `is_kinase_of_interest`

`network_data.edges()` returns the edge list with:

- `source_kinase`
- `target_kinase`
- `correlation`
- `shared_module_count`
- `shared_modules` (serialized module-id list such as `[1,2]`)
- `source_is_kinase_of_interest`
- `target_is_kinase_of_interest`

`network_data.adjacency()` returns the kinase correlation matrix used to derive the edge list.

### Export Contract

```python
network_data.to_frames(copy: bool = True) -> dict[str, pd.DataFrame]
network_data.to_csv(directory: str | Path) -> dict[str, Path]
```

The named network exports are:

- `signalome_network_nodes`
- `signalome_network_edges`
- `signalome_network_adjacency`

## `SignalomeMapData`

Build map-ready plotting data from a signalome result:

```python
map_data = result.to_map_data()
```

This step is intentionally separate from signalome construction. `SignalomeMapData` contains deterministic coordinate and relationship tables that plotting code can consume without adding matplotlib or graph-library dependencies to the computational core.

### Stable Access Paths

```python
map_data.modules(copy: bool = True) -> pd.DataFrame
map_data.sites(copy: bool = True) -> pd.DataFrame
map_data.kinases(copy: bool = True) -> pd.DataFrame
map_data.links(copy: bool = True) -> pd.DataFrame
```

`map_data.modules()` returns the module-position table with:

- index: `module_id`
- `x`
- `y`
- `n_sites`
- `n_proteins`
- `dominant_kinase`
- `dominant_share_percent`

`map_data.sites()` returns the site-position table with:

- index: `site_id`
- `protein_id`
- `module_id`
- `top_kinase`
- `top_score`
- `x`
- `y`
- `module_x`
- `module_y`
- `position_in_module`
- `expression_mean`
- `expression_std`

`map_data.kinases()` returns the kinase-position table with:

- index: `kinase`
- `x`
- `y`
- `base_x`
- `module_count`
- `total_share_percent`
- `degree`
- `n_substrates`
- `is_kinase_of_interest`

`map_data.links()` returns the kinase-to-module plotting table with:

- `kinase`
- `module_id`
- `share_percent`
- `kinase_x`
- `kinase_y`
- `module_x`
- `module_y`
- `is_kinase_of_interest`

### Export Contract

```python
map_data.to_frames(copy: bool = True) -> dict[str, pd.DataFrame]
map_data.to_csv(directory: str | Path) -> dict[str, Path]
```

The named map exports are:

- `signalome_map_modules`
- `signalome_map_sites`
- `signalome_map_kinases`
- `signalome_map_links`

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

# preferred predMat result object
PredMatWorkflowResult.pred_mat_result -> PredMatResult

PredMatResult(
    data_frame: pd.DataFrame,
)

PredMatResult.to_frame(copy: bool = True) -> pd.DataFrame
PredMatResult.to_csv(path: str | Path, index_label: str = "phosphosite") -> Path

`SiteMatrixResult` is a separate preprocessing result for the corrected phosphosite site matrix. It is not the prediction-matrix contract. The preferred prediction-matrix result is `PredMatResult`.

SignalomeResult(
    scoring_matrix: pd.DataFrame,
    pred_mat: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    modules: SignalomeModules,
    assignments: SignalomeAssignments,
    network: SignalomeKinaseNetwork,
    expanded_signalomes: dict[str, ExpandedSignalome],
)

SignalomeResult.to_map_data() -> SignalomeMapData
SignalomeResult.to_network_data() -> SignalomeNetworkData

SignalomeNetworkData(
    adjacency_matrix: pd.DataFrame,
    node_table: pd.DataFrame,
    edge_table: pd.DataFrame,
)

SignalomeMapData(
    module_positions: pd.DataFrame,
    site_positions: pd.DataFrame,
    kinase_positions: pd.DataFrame,
    kinase_module_links: pd.DataFrame,
)

SignalomeMapData.to_frames(copy: bool = True) -> dict[str, pd.DataFrame]
SignalomeMapData.to_csv(directory: str | Path) -> dict[str, Path]

SignalomeNetworkData.to_frames(copy: bool = True) -> dict[str, pd.DataFrame]
SignalomeNetworkData.to_csv(directory: str | Path) -> dict[str, Path]

# stable signalome result views
SignalomeResult.modules -> SignalomeModules
SignalomeResult.assignments -> SignalomeAssignments
SignalomeResult.network -> SignalomeKinaseNetwork

SignalomeResult.to_frames(
    copy: bool = True,
    include_inputs: bool = False,
) -> dict[str, pd.DataFrame]
SignalomeResult.to_csv(
    directory: str | Path,
    include_inputs: bool = False,
) -> dict[str, Path]

KinasePredictionResult(
    pred_matrix: pd.DataFrame,
    substrate_list: dict[str, list[str]],
)

# preferred predMat result object on the lower-level result
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
