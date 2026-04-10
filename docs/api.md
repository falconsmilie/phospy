# Public API Reference

PhosPy has no HTTP API. The supported surface is the Python API below plus the `phospy` CLI.

## Root Imports

```python
from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseActivityAnalyzer,
    KinaseWorkflow,
    PhosphoDataset,
    ReferenceBundle,
    ReferenceBundleProvenance,
    ReferenceBundleSourceMetadata,
    ReferenceProvider,
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

- `PhosphoDataset` for validated total and phospho inputs plus core preprocessing
- `AnalysisReadyPhosphoDataset` for the analysis-ready phosphosite boundary between preprocessing and kinase inference
- `KinaseActivityAnalyzer` for analysis from an existing `predMat`
- `ReferenceBundle` for the kinase-prior boundary between reference resolution and workflow execution
- `ReferenceProvider` for providers that resolve species and reference selections into a `ReferenceBundle`
- `PhosRPipeline` for file loading, preprocessing, optional kinase analysis, and publishing
- `PredMatWorkflow` for the supported `predMat` generation path
- `KinaseWorkflow` for the fuller native scoring and prediction path
- `SignalomeWorkflow` for downstream signalome construction

## Shared Input Rules

### File formats

- total input: TSV
- phospho input: TSV
- `predMat`: CSV, with the first column used as the phosphosite index

### Default required columns

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

### Common validation rules

- file-loaded headers are normalised to lowercase snake case
- duplicate cleaned headers are rejected
- required identifier columns must not be null
- numeric sample columns are coerced to numeric and fail if non-numeric values remain
- `localization_prob` and `predMat` scores must stay in `[0, 1]`
- `gene_p_site` must split cleanly into gene and site parts such as `BTK_Y551`

## `DatasetSchema`

Use `DatasetSchema` when your sample columns do not use the defaults.

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
- structural columns stay fixed package columns

## `PhosphoDataset`

Use `PhosphoDataset` when you want one validated dataset owner and the helpers bound to it.

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
- `phospho_df` must include the required structural columns plus the schema's `phospho_cols`
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

Use this for the normal file-based preprocessing path.

### Data access

```python
dataset.total_df_copy -> pd.DataFrame
dataset.phospho_df_copy -> pd.DataFrame
dataset.total_df_live -> pd.DataFrame
dataset.phospho_df_live -> pd.DataFrame
dataset.copy_inputs() -> tuple[pd.DataFrame, pd.DataFrame]
```

- `*_copy` and `copy_inputs()` return detached copies
- `*_live` returns the dataset-owned DataFrames directly

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

Rules:

- `localization_threshold` and `max_unmatched_fraction` must be in `[0, 1]`
- `min_observed` must be `>= 1`
- `config` can replace the scalar arguments, but cannot be mixed with overridden scalar values
- `max_unmatched_fraction=0.0` is strict mode and rejects silent row loss during protein correction

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

### `dataset.preprocessing.to_analysis_ready(...)`

```python
dataset.preprocessing.to_analysis_ready(
    result: CoreProcessingResult,
    *,
    source: str = "dataset preprocessing",
) -> AnalysisReadyPhosphoDataset
```

Use this as the supported adapter from the bound preprocessing lane into `AnalysisReadyPhosphoDataset` when you already have a `CoreProcessingResult` from the same dataset.

Rules:

- `result` must be a `CoreProcessingResult`
- the adapter reuses the existing site-matrix output from preprocessing rather than rebuilding it differently
- schema and comparisons are taken from the bound dataset preprocessing facade, so callers do not need to pass them back in manually

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

Use this when you want the normal preprocessing path and the analysis-ready boundary in one supported step.

Example:

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
analysis_ready = dataset.preprocessing.run_analysis_ready(
    max_unmatched_fraction=0.1,
)

phospho_matrix = analysis_ready.phospho_matrix
site_sequences = analysis_ready.site_sequences
```

## `AnalysisReadyPhosphoDataset`

Use `AnalysisReadyPhosphoDataset` as the explicit boundary between preprocessing and kinase inference.

It owns:

- `phospho_matrix`: the aligned phosphosite matrix keyed by stable site ID
- `site_metadata`: aligned site metadata keyed by the same site ID
- `site_sequences`: aligned site-centred sequences keyed by the same site ID
- `phospho_corrected`: the corrected phosphosite table the analysis matrix was derived from
- `provenance`: preprocessing provenance including schema, comparisons, row counts, and site-matrix row-drop diagnostics

### `from_core_processing_result(...)`

```python
AnalysisReadyPhosphoDataset.from_core_processing_result(
    result: CoreProcessingResult,
    *,
    schema: DatasetSchema,
    comparisons: Sequence[tuple[str, str]] | None = None,
    source: str = "core preprocessing",
) -> AnalysisReadyPhosphoDataset
```

Use this for advanced cases when preprocessing has already completed outside the bound dataset adapter path and you want one owned, analysis-ready dataset object instead of passing lower-level tables around.

Rules:

- `result` must be a `CoreProcessingResult`
- `schema` defines which corrected value columns belong to the analysis matrix
- the resulting `phospho_matrix`, `site_metadata`, and `site_sequences` must stay exactly aligned on the same unique `site_id` index
- the object owns deep copies of its matrix, metadata, sequences, and corrected source table

Example:

```python
from phospy import AnalysisReadyPhosphoDataset, PhosphoDataset

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)
analysis_ready = AnalysisReadyPhosphoDataset.from_core_processing_result(
    core,
    schema=dataset.schema,
    comparisons=dataset.comparisons,
)

phospho_matrix = analysis_ready.phospho_matrix
site_sequences = analysis_ready.site_sequences
```

This boundary keeps preprocessing outputs reusable without coupling new workflow layers directly to raw input tables or to the full `PhosphoDataset` workspace.

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

Rules:

- `corrected_df` must include `gene_p_site_col`, `sequence_col`, and the schema's `corrected_cols`
- corrected value columns must be numeric
- rows with missing sequence data or incomplete corrected values may be dropped
- if the same phosphosite appears more than once after correction, PhosPy keeps the row with the highest mean corrected signal

Returns `SiteMatrixResult` with:

- `phosr_input`
- `matrix`
- `sequences`
- `row_drop_stats`


## `ReferenceBundle` and `ReferenceProvider`

Use `ReferenceBundle` as the explicit kinase-prior boundary between reference resolution and workflow setup.

A bundle owns:

- `substrate_map`: kinase-to-site mappings used for profile scoring
- `motif_sequences`: kinase-to-sequence-window mappings used for motif scoring
- `species`: the resolved species label for the reference inputs
- `source_metadata`: typed source metadata describing where the reference came from
- `provenance`: typed provenance describing which provider resolved it

### Constructor

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

Validation highlights:

- `substrate_map` must not be empty
- `motif_sequences` must not be empty
- every kinase entry in both mappings must contain at least one value
- the kinase sets in `substrate_map` and `motif_sequences` must match exactly
- `species` must not be empty
- source metadata and provenance fields must not be empty
- motif sequences are validated as a unit through the package motif-library validation path

Example:

```python
from phospy import (
    ReferenceBundle,
    ReferenceBundleProvenance,
    ReferenceBundleSourceMetadata,
)

reference_bundle = ReferenceBundle(
    substrate_map={"KINASE_A": ["SITE_1", "SITE_2"]},
    motif_sequences={"KINASE_A": ["QQSQQ", "QQTQQ"]},
    species="human",
    source_metadata=ReferenceBundleSourceMetadata(
        source="bundled",
        reference="phosr-like",
        version="2026.04",
    ),
    provenance=ReferenceBundleProvenance(
        provider="BundledReferenceProvider",
        notes=("validated",),
    ),
)
```

Use `ReferenceProvider` as the protocol for objects that resolve a species and reference selection into a `ReferenceBundle`:

```python
class ReferenceProvider(Protocol):
    def resolve(
        self,
        *,
        species: str,
        reference: str = "auto",
    ) -> ReferenceBundle: ...
```

This keeps kinase-prior assembly out of user code and gives later bundled or downloaded providers one stable contract.

## `CorePreprocessingConfig`

Use `CorePreprocessingConfig` when you want one configuration object instead of scalar preprocessing options.

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

You can pass it to `dataset.preprocessing.run(config=...)` or `PhosRPipeline(..., preprocessing_config=...)`.

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
- `parquet` with a parquet engine such as `pyarrow`

Written core outputs:

- `df_total_unique`
- `df_total_filtered`
- `df_phospho_filtered`
- `df_phospho_corrected`
- `phosr_input`
- `mat_phospho_corrected`
- `site_sequences`

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
    pred_mat: pd.DataFrame | PredMatResult,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult
```

Rules:

- `pred_mat` may be an in-memory `predMat` `DataFrame` or a `PredMatResult`
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

Writes:

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
) -> PhosRPipeline
```

Rules:

- `dataset` must be a `PhosphoDataset`
- `pred_mat`, when provided, may be a valid `DataFrame` or a `PredMatResult`
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
- `pred_mat_path`, when provided, is loaded and validated during pipeline construction
- overlap between `predMat` and the processed site matrix is checked at runtime after preprocessing

### `run(...)`

```python
pipeline.run(outdir: str | Path | None = None) -> CoreOutputs
```

- when `outdir` is set, files are published to disk
- when `outdir` is `None`, only the in-memory result bundle is returned

`CoreOutputs` exposes:

- `core`
- `kinase_activity`

When `outdir` is set, the pipeline also writes `run_manifest.json`.

## `PredMatWorkflow`

Use `PredMatWorkflow` when your goal is to generate a `predMat` from one supported public workflow.

### Constructor

```python
PredMatWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: str = "default",
)
```

Supported public `svm_mode` values:

- `"default"` for the recommended stable native path
- `"r_parity"` for the supported parity-oriented learner, sampling, and final-scoring preset

### `run(...)`

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
    svm_mode: str | None = None,
) -> PredMatWorkflowResult
```

Validation highlights:

- `phospho_matrix` must be a numeric phosphosite matrix with a unique, non-null index
- pass either explicit `substrate_map` / `motif_sequences` inputs or one `reference_bundle`
- `substrate_map` must not be empty and must overlap with `phospho_matrix`
- `site_sequences`, when supplied, must be a mapping or pandas `Series` keyed by phosphosite ID
- `motif_sequences` is required unless `allow_profile_only_fallback=True`
- when `motif_sequences` is provided, `site_sequences` is also required
- motif-aware validation only requires sequence coverage for phosphosites that actually participate in scoring and prediction
- `min_substrates`, `min_motif_size`, `ensemble_size`, `top`, `inclusion`, and `n_iterations` must be `>= 1`
- `score_threshold` must be in `[0, 1]`

Returns `PredMatWorkflowResult` with:

- `scoring_result`
- `prediction_result`
- `pred_mat_result`

When thresholds are too strict and no kinase candidates qualify, PhosPy raises `NoCandidateKinasesError` instead of returning an empty invalid `predMat`.

## `PredMatResult`

`PredMatResult` is the stable in-memory and export contract for a generated `predMat`.

```python
PredMatResult(data_frame: pd.DataFrame)
```

Access and export:

```python
result.pred_mat_result.data_frame
result.pred_mat_result.to_frame(copy: bool = True) -> pd.DataFrame
result.pred_mat_result.to_csv(path: str | Path, index_label: str = "phosphosite") -> Path
```

Notes:

- rows are phosphosite identifiers
- columns are kinase identifiers
- values are prediction scores in `[0, 1]`
- `data_frame` returns the owned in-memory table
- `to_frame()` returns a detached copy by default

## `KinaseWorkflow`

Use `KinaseWorkflow` when you want the fuller native scoring and prediction path, including intermediate profile and motif outputs.

### Constructor

```python
KinaseWorkflow(
    flank_size: int = 7,
    kernel: str = "rbf",
    svm_mode: str = "default",
)
```

Supported public `svm_mode` values:

- `"default"` for the recommended stable native path
- `"r_parity"` for the supported parity-oriented preset

### `run(...)`

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
    svm_mode: str | None = None,
) -> KinaseWorkflowResult
```

You can supply kinase priors either as explicit `substrate_map` and `motif_sequences` inputs or as one `reference_bundle`. Mixed usage is rejected so the workflow boundary stays unambiguous.

Returns `KinaseWorkflowResult` with:

- `profile_result`
- `motif_result`
- `scoring_result`
- `prediction_result`

The most commonly used fields on `prediction_result` are:

- `pred_matrix`
- `substrate_list`
- `pred_mat_result`

## `SignalomeWorkflow`

Use `SignalomeWorkflow` when you already have aligned scoring and prediction outputs and want one validated signalome step.

### `run(...)`

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

For larger datasets, set `module_count` explicitly if you want to skip the extra automatic module-selection scoring pass.

## `SignalomeResult`

`SignalomeWorkflow.run(...)` returns `SignalomeResult`.

Preferred access paths:

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
```

Export helpers:

```python
result.to_frames(copy: bool = True, include_inputs: bool = False) -> dict[str, pd.DataFrame]
result.to_csv(directory: str | Path, include_inputs: bool = False) -> dict[str, Path]
```

Default exports:

- `signalome_modules`
- `kinase_module_relationships`
- `site_assignments`
- `protein_assignments`
- `kinase_network_nodes`
- `kinase_network_edges`
- `kinase_correlation_matrix`

Pass `include_inputs=True` to also export the aligned `scoring_matrix`, `pred_mat`, and `expression_matrix` used to build the signalome.

Older convenience attributes such as `signalome_modules`, `site_assignments`, and `kinase_network_edges` still exist, but prefer `modules`, `assignments`, and `network` for new code.

## `SignalomeMapData`

Build map-ready plotting data from a signalome result:

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

Named CSV exports:

- `signalome_map_modules`
- `signalome_map_sites`
- `signalome_map_kinases`
- `signalome_map_links`

## `SignalomeNetworkData`

Build graph-friendly kinase-network tables from a signalome result:

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

Named CSV exports:

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
