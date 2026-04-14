# Validation Guide

This is the short version of what PhosPy checks.

For method signatures and defaults, see [`api.md`](api.md). For parity-sensitive behaviour, see [`parity.md`](parity.md).

## File Types

- total input: TSV
- phospho input: TSV
- `predMat`: CSV, with the first column used as the phosphosite index

When total and phospho tables are loaded from files, headers are cleaned to lowercase snake case first. Loading fails if two raw headers collapse to the same cleaned name.

## Required Table Shape

### Total Table

Required columns:

- `genes`
- `group1` to `group6` by default, or your schema's `total_cols`

### Phospho Table

Required columns:

- `uid`
- `gene_names`
- `gene_p_site`
- `localization_prob`
- `centralized_sequence`
- `p_group1` to `p_group6` by default, or your schema's `phospho_cols`

`gene_p_site` must split cleanly into gene and site parts such as `BTK_Y551`.
The site token must follow `<letters><digits>`, such as `S339`, `T308`, or `Y551`.

### Canonical Phosphosite Identifier Contract

PhosPy uses one canonical phosphosite ID format across preprocessing and prediction outputs:

- `ENTITY;SITE;` (for example `BTK;Y551;`)
- `ENTITY` is a non-empty token without semicolons
- `SITE` must match `<letters><digits>` (for example `S123`)

Signalome workflows use this same canonical format when deriving protein grouping from IDs.
If your IDs are not canonical, pass an explicit `site_to_protein` mapping (or use
`run_from_analysis_ready(...)` so mapping is resolved from aligned site metadata).

### `predMat`

A valid `predMat` must have:

- phosphosite IDs as the index, such as `BTK;Y551;`
- kinase names as columns
- numeric scores in `[0, 1]`
- a unique, non-null index

`NaN` values are allowed for missing or unusable kinase scores. Infinite values are rejected.

## Checks You Will Hit Most Often

- required columns must exist
- required identifier columns must not be null
- numeric sample columns must still be numeric after coercion
- `localization_prob` must be in `[0, 1]`
- `predMat` scores must be in `[0, 1]` where present
- file paths must exist and point to files
- comparison pairs must use known schema groups and must not be duplicated
- downstream kinase analysis needs overlap between `predMat` and the phosphosite matrix
- that overlap must cover at least 50% of the phosphosite matrix
- native workflow runs need overlap across the matrix, substrate map, and sequence inputs
- motif-aware runs need both `motif_sequences` and `site_sequences`, unless `allow_profile_only_fallback=True`
- motif-aware validation only requires sequence coverage for phosphosites that are actually scored and predicted
- signalome assignment needs fully finite aligned `predMat` values because each row needs a concrete top kinase

## Useful Behaviour to Know

- by default, protein correction allows no silent phosphosite row loss
- site-matrix building always drops rows with missing sequence data
- site-matrix missing-value handling is explicit through `SiteMatrixPolicy.missing_data_policy`:
  - `drop_any_missing` keeps only complete corrected rows (legacy/default behavior)
  - `retain_missing` keeps partial rows and preserves `NaN` in the matrix
  - `require_min_observed_values` keeps rows meeting a minimum observed-value count
- if the same phosphosite appears more than once after correction, PhosPy keeps the row with the highest mean corrected signal
- when prediction thresholds are too strict, prediction raises `NoCandidateKinasesError` instead of returning an empty invalid `predMat`

## Good Starting Points

Start at the boundary closest to your input:

- `PhosphoDataset.from_files(...)` for the standard preprocessing path
- `SimpleKinaseWorkflow.run(...)` for end-to-end preprocessing, prediction, and kinase activity
- `KinaseActivityAnalyzer.run(...)` for analysis from an existing `predMat`
- `SignalomeWorkflow.run(...)` for downstream signalome construction

## Quick Troubleshooting

| Problem | Usually means | Good next step |
| --- | --- | --- |
| Missing required columns | Your input headers do not match the expected schema | Check cleaned column names and, if needed, pass a custom `DatasetSchema` |
| `predMat` overlap error | The phosphosite IDs do not line up between the matrix and `predMat` | Confirm both sides use canonical `ENTITY;SITE;` IDs |
| `NoCandidateKinasesError` | Thresholds or inclusion settings filtered out every kinase | Relax the prediction config and rerun |
| Sequence coverage error | `site_sequences` or `motif_sequences` do not cover the scored sites | Check keys and confirm you passed the right reference inputs |
| Signalome top-kinase failure | Your aligned prediction values contain non-finite rows | Clean or regenerate the prediction output before signalome construction |

## Recommended Kinase-Scoring Path

1. Load phospho (and optional total) data with supported schema columns.
2. Run `SimpleKinaseWorkflow.run(...)`.
3. Read `result.pred_mat_result` and `result.kinase_activity_result`.
4. Export with `result.pred_mat_result.to_csv("predMat.csv")` when needed.

Runnable example: [`../examples/simple_workflow_demo.py`](../examples/simple_workflow_demo.py)
