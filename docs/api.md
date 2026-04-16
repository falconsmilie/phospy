# API Guide

PhosPy does not have an HTTP API. The supported public surface is the Python API plus the `phospy` CLI.

This page starts with the simplest supported path, then links to the lower-level options.

## Which Entry Point Should You Use?

| Goal | Start here |
| --- | --- |
| Clean and prepare phospho data | `PhosphoDataset` |
| Run the common end-to-end workflow | `SimpleKinaseWorkflow` |
| Build a signalome from workflow outputs | `SignalomeWorkflow` |

`KinaseWorkflow` and `PhosRPipeline` are internal orchestration helpers and are not part of the supported public API.

## Supported Result Contracts

Supported public lane:

- `phospy.api.SimpleKinaseWorkflow` and `phospy.api.SignalomeWorkflow`
- `phospy.api.workflow_results.SimpleKinaseWorkflowResult`

Internal helpers such as `phospy.internal.kinase_workflows.KinaseWorkflow` exist for package implementation and contributor workflows, but they are not part of the supported public result contract.

`predMat` is part of the prediction result contract (`prediction_result.pred_mat_result`). You do not need a separate predMat workflow.

### Public common lane result: `SimpleKinaseWorkflowResult`

`SimpleKinaseWorkflow().run(...)` returns one `SimpleKinaseWorkflowResult` object:

```text
SimpleKinaseWorkflowResult
|- analysis_ready_dataset (AnalysisReadyPhosphoDataset)
|- reference_bundle (ReferenceBundle)
|- scoring_result (KinaseScoringResult)
|  |- profile_scores
|  |- combined_scores
|  `- weights
|- prediction_result (KinasePredictionResult)
|  |- pred_matrix
|  |- pred_mat_result
|  |- substrate_list
|  `- optional traces (debug_traces, trace_level, trace_sink)
`- kinase_activity_result (KinaseActivityResult)
   |- weighted_activity
   |- ksea_scores
   |- ksea_counts
   |- target_counts
   |- target_table
   `- overlap_summary
```

Primary result access paths on `SimpleKinaseWorkflowResult`:

- `prediction_result.pred_mat_result` is the canonical predMat contract
- `scoring_result.profile_scores`, `scoring_result.combined_scores`, and `scoring_result.weights` expose scoring outputs
- `prediction_result.substrate_list` exposes predicted substrate memberships

`SimpleKinaseWorkflowResult` convenience properties (`profile_scores`, `combined_scores`, `weights`, `substrate_list`) are detached reads; use `to_owned_*` or `to_mutable_*_unsafe` accessors when you intentionally need shared state.

Common access pattern:

```python
from phospy.api import SimpleKinaseWorkflow

with SimpleKinaseWorkflow().run(...) as result:
    prediction_result = result.prediction_result
    scoring_result = result.scoring_result
    kinase_activity_result = result.kinase_activity_result

    pred_mat_result = prediction_result.pred_mat_result
    pred_mat = pred_mat_result.to_owned_frame()
    pred_matrix = prediction_result.pred_matrix
    combined_scores = scoring_result.combined_scores
    weighted_activity = kinase_activity_result.weighted_activity
```

Use `prediction_result.pred_mat_result` for the stable predMat table contract. Use `prediction_result` when you need full prediction payload details (for example `substrate_list` or optional traces).

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
    prediction_result = result.prediction_result
    scoring_result = result.scoring_result
    kinase_activity_result = result.kinase_activity_result

    pred_mat_result = prediction_result.pred_mat_result
    pred_mat = pred_mat_result.to_owned_frame()
    pred_matrix = prediction_result.pred_matrix
    combined_scores = scoring_result.combined_scores
    weighted_activity = kinase_activity_result.weighted_activity
```

The returned `result` includes:

- `analysis_ready_dataset`
- `reference_bundle`
- `scoring_result`
- `prediction_result`
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
    kinases_of_interest=list(result.prediction_result.pred_mat_result.kinase_names[:2]),
    config=SignalomeRunConfig(signalome_cutoff=0.5),
)
```

### Save and reload workflow output bundles

`SimpleKinaseWorkflowResult` supports an explicit reproducible output-bundle format.

```python
from phospy.api import (
    PredictionRunConfig,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowConfigSnapshot,
)

prediction_config = PredictionRunConfig(
    min_substrates=1,
    min_motif_size=1,
    ensemble_size=2,
    top=3,
    inclusion=2,
    n_iterations=2,
    random_state=7,
)

with SimpleKinaseWorkflow().run(
    phospho="study_phospho.tsv",
    total="study_total.tsv",
    species="rat",
    prediction_config=prediction_config,
) as result:
    snapshot = SimpleKinaseWorkflowConfigSnapshot.from_workflow_inputs(
        prediction_config=prediction_config,
    )
    bundle_dir = result.save_output_bundle(
        "out/workflow_bundle",
        config_snapshot=snapshot,
    )

metadata = result.load_output_bundle_metadata(bundle_dir)
bundle = result.load_output_bundle(
    bundle_dir,
    table_ids=("pred_mat", "kinase_activity_matrix"),
)
```

The bundle manifest includes:

- workflow type (`SimpleKinaseWorkflowResult`)
- config snapshot
- reference identity and provenance
- output inventory (`table_id` to relative file path and value type)

### Ownership and mutability

Result objects in prediction, workflow, and signalome modules follow one naming rule:

- `to_<resource>()`: detached safe copy
- `to_owned_<resource>()`: cheap shared owned state (no copy)
- `to_mutable_<resource>_unsafe()`: explicit mutable shared state

Example (`PredMatResult`):

```python
safe = result.pred_mat_result.to_frame()
owned = result.pred_mat_result.to_owned_frame()
mutable = result.pred_mat_result.to_mutable_frame_unsafe()
```

`PhosphoDataset` is still a mutable workspace, but safe read access is detached by default:

- `dataset.inputs`, `dataset.total_df_copy`, `dataset.phospho_df_copy`, `dataset.copy_inputs()`

Unsafe mutable access is advanced and intentional:

- `dataset.to_mutable_frames_unsafe()`
- Warning: mutating returned frames mutates the dataset's owned internal state

`SignalomeResult` and nested signalome wrappers (`modules`, `assignments`, `network`) also default to detached reads:

- table properties (for example `signalome.site_assignments`, `signalome.signalome_modules`)
- `to_frames()` (always detached)
- `to_owned_frames()` (cheap shared owned access)

Unsafe mutable access is advanced and intentional:

- `signalome.to_mutable_frames_unsafe()`
- `signalome.to_mutable_expanded_signalomes_unsafe()`
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
)
```

```python
SimpleKinaseWorkflow.from_execution_graph(
    execution_graph: SimpleKinaseExecutionGraph,
) -> SimpleKinaseWorkflow
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

Advanced composition seam:

- `phospy.api.simple_workflow_composition.SimpleKinaseExecutionGraph`
- collaborator method contracts:
- `analysis_ready_builder.build(...)`
- `reference_provider.resolve(...)`
- `activity_analyzer.run(...)`
- `workflow_executor.validate_request(...)`
- `workflow_executor.execute_validated_request(...)`

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

`run_from_analysis_ready(...)` accepts only `AnalysisReadyPhosphoDataset` inputs (for example `result.analysis_ready_dataset` from `SimpleKinaseWorkflow.run(...)`).

`run_from_analysis_ready(...)` defaults to strict metadata resolution and requires a `protein_id` metadata column. To opt in to metadata fallback columns, set `metadata_fallback_policy="metadata"` and provide `metadata_protein_columns`. Gene-symbol fallback is disabled by default and must be explicitly enabled with `allow_gene_symbol_fallback=True`.

Site-to-protein resolution now exposes structured diagnostics:

- `analysis_ready.resolve_site_to_protein_mapping_with_diagnostics(...)` returns `SiteToProteinResolutionResult(mapping, diagnostics)`.
- `analysis_ready.resolve_site_to_protein_mapping(...)` still returns a `pd.Series`, with diagnostics attached at `series.attrs["site_to_protein_resolution_diagnostics"]`.
- `SignalomeWorkflow.run_from_analysis_ready(...)` attaches the same diagnostics to `signalome_result.site_to_protein_resolution_diagnostics`.

Diagnostics include the chosen identifier column, fallback mode, checked/incomplete metadata columns, and ambiguous-identifier counts.

## Prediction Package Surface (`phospy.prediction`)

`phospy.prediction` exports only a stable, intentionally small default surface:

- Stable default exports (`from phospy.prediction import ...`): `KinasePredictor`, `KinasePredictionResult`, `PredMatResult`, `KinaseScorer`, `KinaseScoringResult`, `KinaseMotifScorer`, `MotifScoringResult`, `KinaseProfilePolicy`

Everything else must be imported from concrete modules. Common advanced imports:

- `phospy.prediction.contracts.EnsemblePredictorContract`
- `phospy.prediction.traces.PredictionSamplingTrace`
- `phospy.prediction.scoring.combine_profile_and_motif_scores`
- `phospy.prediction.candidates.build_candidate_substrate_list`

Internal-only helper types in low-level modules (for example
`phospy.prediction.execution` and `phospy.prediction.trace_runtime`) are
contributor seams and may change between releases unless promoted to this
default package surface.

Advanced prediction result/debug structures:

- `KinasePredictionResult.debug_traces` values are `KinasePredictionDebugTrace`
- `KinasePredictionDebugTrace.ensemble_traces` contain `AdaptiveSamplingEnsembleTrace`
- `AdaptiveSamplingEnsembleTrace.iterations` contain `AdaptiveSamplingIterationTrace`
- these advanced trace types are for diagnostics and replay workflows; production integrations should prefer `pred_mat_result`, `substrate_list`, and `prediction_debug_trace_tables(...)`

## Advanced Extension Contract: Custom Ensemble Predictors

If you inject a custom `ensemble_predictor` into `phospy.prediction.KinasePredictor`,
implement `phospy.prediction.contracts.EnsemblePredictorContract.predict_kinase(...)` and return
`KinasePredictionBatch` objects that satisfy all output constraints:

- `batch.kinase` must match the requested kinase argument.
- `batch.score_values` must be a 1D numeric vector with exactly one value per
  row in `feature_mat`.
- score values must be finite (`NaN`/`+inf`/`-inf` are rejected).
- `batch.score_index` uses the `label_aligned` contract:
  it must contain exactly the same phosphosite labels as `feature_mat.index`
  (reordered labels are allowed and are aligned by label).

These checks run at the predictor boundary before aggregation; invalid outputs
raise `CustomPredictorOutputError` with the failing kinase and reason.

## Preferred Imports

```python
from phospy.api import (
    DatasetLoadOptions,
    KinaseActivityConfig,
    PredictionRunConfig,
    SimpleKinaseWorkflowConfigSnapshot,
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
- `CustomPredictorOutputError`: custom ensemble predictor returned invalid batch output

## Related Docs

- [`validation.md`](validation.md)
- [`parity.md`](parity.md)
- [`../examples/simple_workflow_demo.py`](../examples/simple_workflow_demo.py)
- [`../examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py)
