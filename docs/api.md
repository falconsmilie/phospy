# Public API Guide

This is the compact guide to the supported Python API.

PhosPy does not expose HTTP endpoints. The supported public surface is the Python API on this page together with
the small `phospy` command-line interface described in the README.

If you stay within the classes on this page, you are using the public surface that PhosPy means to support.

## Stable Root Imports

```python
from phospy import (
    CoreOutputs,
    CoreProcessingResult,
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

## `PhosphoDataset`

Use this when you want the core preprocessing steps without the higher-level pipeline wrapper.

### Create it

From files:

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files(
    "total.tsv",
    "phospho.tsv",
    comparisons=[("group1", "group4")],
    phospho_encoding="utf-8",  # default is utf-8 when omitted
)
```

From in-memory data frames:

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset(total_df=total_df, phospho_df=phospho_df)
```

From already-validated frames:

```python
from phospy import PhosphoDataset

validated = PhosphoDataset.from_validated_inputs(
    total_df=total_df,
    phospho_df=phospho_df,
)
```

A practical detail: `from_files(...)` cleans file headers to lowercase snake case before validation. In-memory data
frames are validated as provided.

### Main methods

#### `process_core(...) -> CoreProcessingResult`

The usual one-call entry point for:

- total-table cleanup
- phospho filtering
- protein correction
- optional pairwise comparison columns
- site-matrix construction

Key parameters:

- `localization_threshold=0.75`
- `min_observed=4`
- `max_unmatched_fraction=0.0`
- `total_sentinel=10.0`
- `phospho_sentinel=12.0`

A practical detail: `max_unmatched_fraction=0.0` is strict. Protein correction fails if the inner join would silently
drop any phosphosite rows.

#### Other supported methods

These are useful when you want the steps separately:

- `prepare_total(...)`
- `prepare_phospho(...)`
- `correct_to_protein(...)`
- `add_pairwise_comparisons(...)`
- `build_site_matrix(...)`
- `write_core_outputs(...)`

You can also customise `total_cols`, `phospho_cols`, and `corrected_cols` at construction time. If you override those
column groups, they must still align by length.

## `CoreProcessingResult`

Returned by `PhosphoDataset.process_core(...)`.

Attributes:

- `total_unique`
- `total_filtered`
- `phospho_filtered`
- `phospho_corrected`
- `site_matrix`

## `SiteMatrixResult`

Available as `core.site_matrix`.

Attributes:

- `phosr_input`: the intermediate table used to build the site matrix
- `matrix`: the phosphosite-by-sample numeric matrix
- `sequences`: centralized sequences indexed by phosphosite ID
- `row_drop_stats`: counts for dropped and retained rows

A useful detail: site-matrix construction drops rows with missing sequences or incomplete corrected values, then
collapses duplicate phosphosites by keeping the row with the highest mean corrected signal.

## `KinaseActivityAnalyzer`

Use this when you already have a phosphosite matrix and a `predMat` and want downstream kinase summaries.
`KinaseActivityAnalyzer` is the sole public orchestration entrypoint for standalone activity analysis. Lower-level
computational helpers remain in `phospy.activities`.

### Preferred usage

```python
from phospy import KinaseActivityAnalyzer

analyzer = KinaseActivityAnalyzer()
result = analyzer.load_and_analyze(
    pred_mat_path="predMat.csv",
    phospho_matrix=phospho_matrix,
)
analyzer.write_outputs(result, outdir="output")
```

### Supported methods

#### `load_pred_mat(pred_mat_path) -> pd.DataFrame`

Loads a prediction matrix from disk and validates it against the public `predMat` schema.

#### `analyze(pred_mat, phospho_matrix, threshold=0.6, min_substrates=3, top_n_substrates=20) -> KinaseActivityResult`

Computes downstream kinase summaries from in-memory `predMat` and phosphosite matrices.

#### `load_and_analyze(pred_mat_path, phospho_matrix, threshold=0.6, min_substrates=3, top_n_substrates=20) -> KinaseActivityResult`

Convenience entry point that combines `load_pred_mat(...)` and `analyze(...)`. This is the preferred entrypoint when your `predMat` lives on disk.

#### `write_outputs(result, outdir) -> None`

Writes the same downstream kinase-analysis tables that `PhosRPipeline` emits when `pred_mat_path` is supplied. Use this to materialize the class-based workflow outputs without going through the pipeline.

Returns:

- weighted activity scores
- KSEA-style scores
- KSEA substrate counts
- kinase target counts
- a long-form kinase target table

The input `predMat` and phosphosite matrix must overlap by at least one row and at least 10% of the phosphosite
matrix.

### Related lower-level functions

These live in `phospy.activities`:

- `build_kinase_target_table(...)`
- `count_predicted_targets(...)`
- `compute_weighted_kinase_activity(...)`
- `compute_ksea_scores(...)`

## `KinaseActivityResult`

Returned by `KinaseActivityAnalyzer.analyze(...)`.

Attributes:

- `weighted_activity`
- `ksea_scores`
- `ksea_counts`
- `target_counts`
- `target_table`

## `PhosRPipeline`

Use this when you want the core preprocessing path and optional `predMat` analysis in one place.

### Create it

```python
from phospy import PhosRPipeline

pipeline = PhosRPipeline.from_files(
    total_path="total.tsv",
    phospho_path="phospho.tsv",
    pred_mat_path="predMat.csv",
    comparisons=[("group1", "group4")],
    phospho_encoding="utf-8",  # default is utf-8 when omitted
    max_unmatched_fraction=0.1,
)
```

### Run it

```python
outputs = pipeline.run(outdir="output")
```

If `pred_mat_path` was supplied, `run()` writes both core outputs and downstream kinase-analysis outputs. If not, it
still writes the core outputs. If you call `run(outdir=None)`, PhosPy still returns the in-memory results and simply
skips file writing.

If you already have a validated `CorePipelineRequest`, you can build the pipeline with `PhosRPipeline.from_request(...)`.
That is the supported entry point when you want path validation separated from pipeline construction.

## `CoreOutputs`

Returned by `PhosRPipeline.run(...)`.

Attributes:

- `core`: a `CoreProcessingResult`
- `kinase_activity`: a `KinaseActivityResult` or `None`

## `KinaseWorkflow`

Use this for the native end-to-end kinase scoring and prediction workflow.

### Create it

```python
from phospy import KinaseWorkflow

workflow = KinaseWorkflow(flank_size=7, svm_mode="default")
```

Constructor parameters:

- `flank_size=7`
- `kernel="rbf"`
- `svm_mode="default"`

### Run it

```python
result = workflow.run(
    phospho_matrix=phospho_matrix,
    substrate_map=substrate_map,
    site_sequences=site_sequences,
    motif_sequences=motif_sequences,
)
```

Key parameters:

- `min_substrates=1`
- `min_motif_size=1`
- `allow_profile_only_fallback=False`
- `ensemble_size=10`
- `top=50`
- `score_threshold=0.8`
- `inclusion=20`
- `n_iterations=5`
- `random_state=None`
- optional per-run `svm_mode`

Important behaviour:

- `site_sequences` must be keyed by phosphosite ID, or supplied as a pandas Series with a phosphosite index.
- `motif_sequences` require matching `site_sequences`.
- If you omit `motif_sequences`, you must set `allow_profile_only_fallback=True`.
- `svm_mode="default"` is the preferred native mode.
- `svm_mode="r_parity"` exists for a narrower learner-seam comparison against committed parity references.

If you already have a validated `KinaseWorkflowRequest`, use `run_request(...)`.

## `KinaseWorkflowResult`

Returned by `KinaseWorkflow.run(...)`.

Attributes:

- `profile_result`
- `motif_result` (`None` in profile-only mode)
- `scoring_result` (`combined_scores` and `weights` are `None` in profile-only mode)
- `prediction_result`

## `KinasePredictionResult`

This is the final prediction bundle from the native workflow.

The two fields most people reach for first are:

- `pred_matrix`
- `substrate_list`

Despite the name, `substrate_list` is a dictionary keyed by kinase name, with each value holding the ranked phosphosite
IDs selected for that kinase.

The full object also carries candidate and debugging information that is useful in seam-level validation and trace
comparison.

## Common Validation Edges

A few checks are worth keeping in mind while you build around the API:

- total input requires `genes` plus `group1` to `group6`
- phospho input requires `uid`, `gene_names`, `gene_p_site`, `localization_prob`, `centralized_sequence`, and
  `p_group1` to `p_group6`
- file-loaded table headers are cleaned to lowercase snake case before validation, and duplicate cleaned names are
  rejected
- `gene_p_site` must be splitable into gene and site parts, such as `BTK_Y551`
- `localization_prob` and `predMat` scores must stay in `[0, 1]`
- by default, protein correction allows no silent row loss
- native workflow inputs must share phosphosite IDs across the matrix, substrate map, and sequence inputs

For the validation and parity story around those checks, see [`validation-and-parity.md`](validation-and-parity.md).
