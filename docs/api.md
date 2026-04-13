# Public API Reference

PhosPy has **no HTTP API**. The supported surface is the Python API plus the `phospy` CLI.

## Start With the Right Import

Preferred imports:

```python
from phospy.activities import KinaseActivityAnalyzer
from phospy.api import (
    KinaseWorkflow,
    PredMatWorkflow,
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

The root package still exposes a small convenience surface, but new code should usually import from the owning domain package.

## Pick an Entry Point

Use:

- `SimpleKinaseWorkflow` for the shortest supported end-to-end lane
- `PhosphoDataset` for validated inputs plus preprocessing
- `KinaseActivityAnalyzer` when you already have a `predMat`
- `PhosRPipeline` for one-shot file loading, preprocessing, optional activity analysis, and publishing
- `PredMatWorkflow` for native `predMat` generation
- `KinaseWorkflow` for the fuller native scoring and prediction path
- `SignalomeWorkflow` for downstream signalome construction

## Root Convenience Surface

`import phospy` intentionally re-exports a small set of high-level types:

- `AnalysisReadyPhosphoDataset`
- `BundledReferenceProvider`
- `KinaseActivityAnalyzer`
- `KinaseWorkflow`
- `PhosphoDataset`
- `PhosRPipeline`
- `PredMatResult`
- `PredMatWorkflow`
- `ReferenceBundle`
- `ReferenceBundleProvenance`
- `ReferenceBundleSourceMetadata`
- `ReferenceProvider`
- `SignalomeMapData`
- `SignalomeNetworkData`
- `SignalomeResult`
- `SignalomeWorkflow`
- `SimpleKinaseWorkflow`

That surface is for convenience, not for broad re-exporting of internals.

## Shared Input Rules

### File Formats

- total input: TSV
- phospho input: TSV
- `predMat`: CSV, with the first column used as the phosphosite index

### Default Required Columns

Total table:

- `genes`
- `group1` to `group6`

Phospho table:

- `uid`
- `gene_names`
- `gene_p_site`
- `localization_prob`
- `centralized_sequence`
- `p_group1` to `p_group6`

### Common Validation Rules

- file-loaded headers are cleaned to lowercase snake case
- duplicate cleaned headers are rejected
- required identifier columns must not be null
- numeric sample columns must stay numeric after coercion
- `localization_prob` and `predMat` scores must be in `[0, 1]`
- `gene_p_site` must split cleanly into gene and site parts such as `BTK_Y551`

For the short checklist, see [`validation.md`](validation.md).

## Shared Data and Preprocessing

### `DatasetSchema`

Use this when your sample columns do not match the defaults.

```python
from phospy.datasets import DatasetSchema

schema = DatasetSchema(
    total_cols=("sample_a", "sample_b"),
    phospho_cols=("p_sample_a", "p_sample_b"),
    corrected_cols=("corrected_a", "corrected_b"),
)
```

Rules:

- all three column groups must have the same length
- the groups must align one-for-one
- structural columns stay fixed package columns

### `PhosphoDataset`

Use `PhosphoDataset` when you want one validated dataset owner plus bound helpers.

Constructor:

```python
PhosphoDataset(
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    *,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
) -> None
```

`from_files(...)`:

```python
PhosphoDataset.from_files(
    total_path: str | Path,
    phospho_path: str | Path,
    phospho_encoding: str | None = None,
    comparisons: Sequence[tuple[str, str]] | None = None,
    schema: DatasetSchema | None = None,
) -> PhosphoDataset
```

Key rules:

- `total_df` must include `genes` and the schema's `total_cols`
- `phospho_df` must include required structural columns plus the schema's `phospho_cols`
- comparison pairs must use known schema groups and must not be duplicated

Useful data access:

```python
dataset.total_df_copy -> pd.DataFrame
dataset.phospho_df_copy -> pd.DataFrame
dataset.total_df_live -> pd.DataFrame
dataset.phospho_df_live -> pd.DataFrame
dataset.copy_inputs() -> tuple[pd.DataFrame, pd.DataFrame]
```

- `*_copy` and `copy_inputs()` return detached copies
- `*_live` returns dataset-owned DataFrames directly

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

Key rules:

- `localization_threshold` and `max_unmatched_fraction` must be in `[0, 1]`
- `min_observed` must be `>= 1`
- `config` can replace the scalar arguments, but should not be mixed with overridden scalar values
- `config.site_matrix_policy` controls duplicate phosphosite collapse during site-matrix creation
- `max_unmatched_fraction=0.0` is strict mode and rejects silent row loss during protein correction

Returns `CoreProcessingResult` with:

- `total_unique`
- `total_filtered`
- `phospho_filtered`
- `phospho_corrected`
- `site_matrix`

### `dataset.preprocessing.to_analysis_ready(...)`

```python
dataset.preprocessing.to_analysis_ready(
    result: CoreProcessingResult,
    *,
    source: str = "dataset preprocessing",
) -> AnalysisReadyPhosphoDataset
```

Use this to adapt an existing preprocessing result from the same dataset into the stable analysis-ready boundary.

### `dataset.preprocessing.run_analysis_ready(...)`

```python
dataset.preprocessing.run_analysis_ready(
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float | int = 10.0,
    phospho_sentinel: float | int = 12.0,
    config: CorePreprocessingConfig | None = None,
    source: str = "dataset preprocessing",
) -> AnalysisReadyPhosphoDataset
```

Use this when you want preprocessing and the analysis-ready boundary in one step.

### `AnalysisReadyPhosphoDataset`

Use this as the explicit handoff between preprocessing and kinase inference.

It owns:

- `phospho_matrix`
- `site_metadata`
- `site_sequences`
- `phospho_corrected`
- `provenance`

`from_core_processing_result(...)`:

```python
AnalysisReadyPhosphoDataset.from_core_processing_result(
    result: CoreProcessingResult,
    *,
    schema: DatasetSchema,
    comparisons: Sequence[tuple[str, str]] | None = None,
    source: str = "core preprocessing",
) -> AnalysisReadyPhosphoDataset
```

### `dataset.site_matrix.build(...)`

```python
dataset.site_matrix.build(
    corrected_df: pd.DataFrame,
    *,
    gene_p_site_col: str = "gene_p_site",
    sequence_col: str = "centralized_sequence",
) -> SiteMatrixResult
```

Use this when you already have a corrected phosphosite table.

Returns `SiteMatrixResult` with:

- `phosr_input`
- `matrix`
- `sequences`
- `row_drop_stats`

## Reference Resolution

### `ReferenceBundle`

Use `ReferenceBundle` as the explicit kinase-prior boundary between reference resolution and workflow execution.

Constructor:

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

Key rules:

- `substrate_map` must not be empty
- `motif_sequences` must not be empty
- the kinase sets in both mappings must match exactly
- `species`, source metadata, and provenance fields must not be empty

### `ReferenceProvider`

Protocol:

```python
class ReferenceProvider(Protocol):
    def resolve(
        self,
        *,
        species: str,
        reference: str = "auto",
    ) -> ReferenceBundle: ...
```

### `BundledReferenceProvider`

Use this when you want the supported packaged reference lane.

```python
from phospy.references import BundledReferenceProvider

provider = BundledReferenceProvider()
reference_bundle = provider.resolve(species="rat", reference="auto")
```

Current bundled support is intentionally narrow:

- supported species: `rat`
- supported references for `rat`: `auto`, `l6`, `l6_native`
- `auto` resolves to `l6_native`

Unsupported species or reference selections fail with explicit validation errors.

## Configuration and Output Helpers

### `CorePreprocessingConfig`

Use this when you want one configuration object instead of individual preprocessing arguments.

```python
from phospy.preprocessing import CorePreprocessingConfig

config = CorePreprocessingConfig(
    localization_threshold=0.75,
    min_observed=4,
    total_sentinel=10.0,
    phospho_sentinel=12.0,
    max_unmatched_fraction=0.1,
)
```

### `CoreOutputWriter`

```python
from phospy.io.writers import CoreOutputWriter

CoreOutputWriter().write(
    result=core,
    outdir="results",
    format="csv",
)
```

Supported formats:

- `csv`
- `tsv`
- `parquet` when a parquet engine such as `pyarrow` is installed

Written core outputs:

- `df_total_unique`
- `df_total_filtered`
- `df_phospho_filtered`
- `df_phospho_corrected`
- `phosr_input`
- `mat_phospho_corrected`
- `site_sequences`

## Supporting Utility

### `KinaseActivityAnalyzer`

Use this when you already have a phosphosite matrix and a `predMat`.

`load_pred_mat(...)`:

```python
analyzer.load_pred_mat(pred_mat_path: str | Path) -> pd.DataFrame
```

Rules:

- the file must be a CSV
- the first column becomes the phosphosite index
- kinase score columns must be numeric
- scores must stay in `[0, 1]`
- the index must be unique and non-null

`run(...)`:

```python
analyzer.run(
    pred_mat: pd.DataFrame | PredMatResult,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult
```

Rules:

- `pred_mat` may be a DataFrame or a `PredMatResult`
- `phospho_matrix` must be a numeric phosphosite matrix with a unique, non-null index
- `threshold` must be in `[0, 1]`
- `min_substrates` and `top_n_substrates` must be `>= 1`
- `pred_mat` and `phospho_matrix` must overlap by at least one row
- that overlap must cover at least 50% of the phosphosite matrix

`write_outputs(...)`:

```python
analyzer.write_outputs(result: KinaseActivityResult, outdir: str | Path) -> None
```

Written files:

- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `kinase_target_table.csv`

## Simple Workflow Lane

### `SimpleKinaseWorkflow`

Use this for the recommended public end-to-end lane.

Constructor:

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

`run(...)`:

```python
SimpleKinaseWorkflow().run(
    *,
    phospho: pd.DataFrame | str | Path,
    species: str,
    total: pd.DataFrame | str | Path | None = None,
    reference: str = "auto",
    phospho_encoding: str | None = None,
    comparisons: Sequence[ComparisonSpec] | None = None,
    schema: DatasetSchema | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float | int = 10.0,
    phospho_sentinel: float | int = 12.0,
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
    profile_policy: KinaseProfilePolicy | None = None,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> SimpleKinaseWorkflowResult
```

Notes:

- `phospho` and `species` are required
- `total` is optional
- when `total` is provided, the workflow reuses the dataset preprocessing path
- the default bundled provider currently supports only `species="rat"` and `reference in {"auto", "l6", "l6_native"}`

Returned result bundle:

- `analysis_ready_dataset`
- `reference_bundle`
- `workflow_result`
- `kinase_activity_result`
- convenience accessors for `pred_mat_result`, `scoring_result`, and `prediction_result`

Example:

```python
from phospy.api import SimpleKinaseWorkflow

result = SimpleKinaseWorkflow().run(
    phospho="study_phospho.tsv",
    total="study_total.tsv",
    species="rat",
    reference="auto",
)

pred_mat = result.pred_mat_result.to_frame()
weighted_activity = result.kinase_activity_result.weighted_activity
```

## File-Based Orchestration

### `PhosRPipeline`

Use this when you want preprocessing and optional kinase analysis in one place.

Constructor:

```python
PhosRPipeline(
    dataset: PhosphoDataset,
    pred_mat: pd.DataFrame | PredMatResult | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
    *,
    manifest_writer: RunManifestWriter | None = None,
    output_publisher: OutputPublisher | None = None,
) -> None
```

`from_files(...)`:

```python
PhosRPipeline.from_files(
    total_path: str | Path,
    phospho_path: str | Path,
    pred_mat_path: str | Path | None = None,
    comparisons: Sequence[ComparisonSpec] | None = None,
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
    *,
    manifest_writer: RunManifestWriter | None = None,
    output_publisher: OutputPublisher | None = None,
) -> PhosRPipeline
```

`run(...)`:

```python
pipeline.run(outdir: str | Path | None = None) -> CoreOutputs
```

The returned bundle exposes:

- `core`
- `kinase_activity`

When `outdir` is set, the pipeline also writes `run_manifest.json`.

## Advanced Native Workflow Lane

### `PredMatWorkflow`

Use this when your goal is to generate a `predMat`.

Constructor:

```python
PredMatWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = "default",
) -> None
```

Supported public `svm_mode` values:

- `"default"` for the recommended stable native path
- `"r_parity"` for the parity-oriented preset

`run(...)`:

```python
workflow.run(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]] | None = None,
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
    reference_bundle: ReferenceBundle | None = None,
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
    profile_policy: KinaseProfilePolicy | None = None,
) -> PredMatWorkflowResult
```

Key rules:

- `phospho_matrix` must be numeric and indexed by unique phosphosite IDs
- pass either `reference_bundle` or explicit reference inputs
- `motif_sequences` is required unless `allow_profile_only_fallback=True`
- when `motif_sequences` is provided, `site_sequences` is also required
- `min_substrates`, `min_motif_size`, `ensemble_size`, `top`, `inclusion`, and `n_iterations` must be `>= 1`
- `score_threshold` must be in `[0, 1]`

Returns:

- `scoring_result`
- `prediction_result`
- `pred_mat_result`

When thresholds are too strict and no kinase candidates qualify, PhosPy raises `NoCandidateKinasesError` instead of returning an empty invalid `predMat`.

### `PredMatResult`

This is the stable in-memory and export contract for a generated `predMat`.

```python
PredMatResult(data_frame: pd.DataFrame)
```

Useful accessors:

```python
result.pred_mat_result.data_frame
result.pred_mat_result.to_frame(copy: bool = True) -> pd.DataFrame
result.pred_mat_result.to_csv(path: str | Path, index_label: str = "phosphosite") -> Path
```

### `KinaseWorkflow`

Use this when you want the fuller native scoring and prediction path, including intermediate profile and motif outputs.

Constructor:

```python
KinaseWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = "default",
) -> None
```

`run(...)`:

```python
workflow.run(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]] | None = None,
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
    reference_bundle: ReferenceBundle | None = None,
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
    profile_policy: KinaseProfilePolicy | None = None,
) -> KinaseWorkflowResult
```

Returns:

- `profile_result`
- `motif_result`
- `scoring_result`
- `prediction_result`

## Signalome Workflow

### `SignalomeWorkflow`

Use this when you already have aligned scoring and prediction outputs and want one validated signalome step.

`run(...)`:

```python
SignalomeWorkflow().run(
    *,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult | PredMatResult,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    site_to_protein: Mapping[str, str] | None = None,
    kinase_network_threshold: float = 0.9,
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
    module_selection_policy: SignalomeModuleSelectionPolicy | None = None,
) -> SignalomeResult
```

Key rules:

- `scoring_result` and `prediction_result` must align to a shared phosphosite index
- `prediction_result` may be `KinasePredictionResult` or `PredMatResult`
- `expression_matrix` must be numeric and indexed by unique phosphosite IDs
- aligned `pred_mat` values must be finite because signalome assignment needs a concrete top kinase per row
- `kinases_of_interest` must not be empty and must be present in the aligned kinase columns
- `kinase_network_threshold` and `signalome_cutoff` must be in `[0, 1]`
- `module_count`, when supplied, must be `>= 1`
- `min_kinase_module_share_percent` must be `>= 0`

For larger datasets, set `module_count` explicitly if you want to skip the extra automatic module-selection scoring pass.

## `SignalomeResult`

Useful accessors:

```python
result.modules
result.assignments
result.network
result.expanded_signalomes
result.to_map_data()
result.to_network_data()
```

Main helpers:

```python
result.modules.to_frame(copy: bool = True) -> pd.DataFrame
result.modules.to_relationship_table(copy: bool = True) -> pd.DataFrame
result.assignments.sites(copy: bool = True) -> pd.DataFrame
result.assignments.proteins(copy: bool = True) -> pd.DataFrame
result.network.adjacency(copy: bool = True) -> pd.DataFrame
result.network.nodes(copy: bool = True) -> pd.DataFrame
result.network.edges(copy: bool = True) -> pd.DataFrame
result.to_frames(copy: bool = True, include_inputs: bool = False) -> dict[str, pd.DataFrame]
result.to_csv(directory: str | Path, include_inputs: bool = False) -> dict[str, Path]
```

Default CSV exports:

- `signalome_modules`
- `kinase_module_relationships`
- `site_assignments`
- `protein_assignments`
- `kinase_network_nodes`
- `kinase_network_edges`
- `kinase_correlation_matrix`

Pass `include_inputs=True` to also export the aligned `scoring_matrix`, `pred_mat`, and `expression_matrix`.

## `SignalomeMapData`

```python
map_data = result.to_map_data()
```

Methods:

```python
map_data.modules(copy: bool = True) -> pd.DataFrame
map_data.sites(copy: bool = True) -> pd.DataFrame
map_data.kinases(copy: bool = True) -> pd.DataFrame
map_data.links(copy: bool = True) -> pd.DataFrame
map_data.to_frames(copy: bool = True) -> dict[str, pd.DataFrame]
map_data.to_csv(directory: str | Path) -> dict[str, Path]
```

CSV exports:

- `signalome_map_modules`
- `signalome_map_sites`
- `signalome_map_kinases`
- `signalome_map_links`

## `SignalomeNetworkData`

```python
network_data = result.to_network_data()
```

Methods:

```python
network_data.nodes(copy: bool = True) -> pd.DataFrame
network_data.edges(copy: bool = True) -> pd.DataFrame
network_data.adjacency(copy: bool = True) -> pd.DataFrame
network_data.to_frames(copy: bool = True) -> dict[str, pd.DataFrame]
network_data.to_csv(directory: str | Path) -> dict[str, Path]
```

CSV exports:

- `signalome_network_nodes`
- `signalome_network_edges`
- `signalome_network_adjacency`

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

## Common Exceptions

You will most often see:

- `RequestValidationError` for invalid public inputs or incompatible options
- `NoCandidateKinasesError` when prediction thresholds leave no kinase candidates to score

## Related Pages

- [`validation.md`](validation.md) for the short validation checklist
- [`parity.md`](parity.md) for parity scope and `svm_mode` guidance
- [`../examples/predmat_workflow_demo.py`](../examples/predmat_workflow_demo.py) for a runnable `predMat` example
- [`../examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py) for a runnable signalome example
