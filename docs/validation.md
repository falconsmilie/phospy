# Validation Quick Guide

This is the short version of what PhosPy checks.

For method signatures and parameter defaults, see [`api.md`](api.md).
For parity-sensitive behaviour, see [`parity.md`](parity.md).

## File Types

- total input: TSV
- phospho input: TSV
- `predMat`: CSV, with the first column used as the phosphosite index

When total and phospho tables are loaded from files, headers are normalised to lowercase snake case first. Loading fails if two raw headers collapse to the same cleaned name.

## Required Table Shape

### Total table

Required columns:

- `genes`
- `group1` to `group6` by default, or your schema's `total_cols`

### Phospho table

Required columns:

- `uid`
- `gene_names`
- `gene_p_site`
- `localization_prob`
- `centralized_sequence`
- `p_group1` to `p_group6` by default, or your schema's `phospho_cols`

`gene_p_site` must split cleanly into gene and site parts such as `BTK_Y551`.

### `predMat`

A valid `predMat` must have:

- phosphosite IDs as the index, such as `BTK;Y551;`
- kinase names as columns
- numeric scores in `[0, 1]`
- a unique, non-null index

## Checks You Will Hit Most Often

- required columns must exist
- required identifier columns must not be null
- numeric sample columns must be numeric after coercion
- `localization_prob` must stay in `[0, 1]`
- `predMat` scores must stay in `[0, 1]` where present
- `predMat` may include `NaN` values to represent missing or unusable kinase scores, but infinite values are rejected
- downstream consumers may tighten that rule when their assignment logic requires fully finite kinase scores
- file paths must exist and point to files
- comparison pairs must use known schema groups and must not be duplicated
- downstream kinase analysis needs overlap between `predMat` and the phosphosite matrix
- that overlap must cover at least 50% of the phosphosite matrix
- native workflow runs need overlap across the matrix, substrate map, and sequence inputs
- motif-aware workflow runs need both `motif_sequences` and `site_sequences`, unless you enable `allow_profile_only_fallback=True`
- motif-aware workflow validation only requires sequence coverage for phosphosites that are actually scored and predicted

## Useful Behaviour to Know

- by default, protein correction allows no silent phosphosite row loss
- site-matrix building can drop rows with missing sequence data or incomplete corrected values
- if the same phosphosite appears more than once after correction, PhosPy keeps the row with the highest mean corrected signal
- when prediction thresholds are too strict, `PredMatWorkflow` raises `NoCandidateKinasesError` instead of returning an empty invalid `predMat`

## Good Starting Point

If you are unsure where validation happens, start here:

- `PhosphoDataset.from_files(...)` for the standard preprocessing path
- `KinaseActivityAnalyzer.run(...)` for analysis from an existing `predMat`
- `PhosRPipeline.from_files(...)` for the file-based one-shot flow
- `PredMatWorkflow.run(...)` for native `predMat` generation
- `SignalomeWorkflow.run(...)` for downstream signalome construction

## Recommended `predMat` Path

1. Load a numeric phosphosite matrix with phosphosite IDs as the index.
2. Load `site_sequences`, `substrate_map`, and `motif_sequences` keyed by the matching phosphosite or kinase identifiers.
3. Run `PredMatWorkflow.run(...)`.
4. Read the result from `result.pred_mat_result`.
5. Export with `result.pred_mat_result.to_csv("predMat.csv")`.

A runnable example lives in [`../examples/predmat_workflow_demo.py`](../examples/predmat_workflow_demo.py).
