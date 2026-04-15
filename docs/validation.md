# Validation Guide

This is the practical version of what PhosPy checks.

For method signatures and defaults, see [`api.md`](api.md). For parity-sensitive behaviour, see [`parity.md`](parity.md).

## Accepted File Types

- total input: TSV
- phospho input: TSV
- `predMat`: CSV, with the first column used as the phosphosite index

When total and phospho tables are loaded from files, headers are cleaned to lowercase snake case first. Loading fails if two original headers collapse to the same cleaned name.

## Minimum Required Columns

### Total table

| Column group | Required by default |
| --- | --- |
| Identifier | `genes` |
| Sample columns | `group1` to `group6` |

If you pass a custom `DatasetSchema`, the expected sample columns come from `schema.total_cols` instead.

### Phospho table

| Column group | Required by default |
| --- | --- |
| Identifiers | `uid`, `gene_names`, `gene_p_site` |
| Metadata | `localization_prob`, `centralized_sequence` |
| Sample columns | `p_group1` to `p_group6` |

If you pass a custom `DatasetSchema`, the expected sample columns come from `schema.phospho_cols` instead.

## Supported Phosphosite Formats

### Raw phospho input

`gene_p_site` must split cleanly into gene and site parts, for example:

- `BTK_Y551`
- `AKT1_T308`
- `PRKACA_S339`

The site token must follow `<letters><digits>`, such as `S339`, `T308`, or `Y551`.

### Internal and `predMat` phosphosite IDs

PhosPy uses one supported phosphosite ID format across preprocessing and prediction outputs:

- `ENTITY;SITE;`
- example: `BTK;Y551;`

Rules:

- `ENTITY` must be present
- `ENTITY` cannot contain semicolons
- `SITE` must follow `<letters><digits>`

Signalome workflows use this same format when they derive protein grouping from phosphosite IDs. If your IDs do not follow it, pass a `site_to_protein` mapping. `run_from_analysis_ready(...)` defaults to strict metadata resolution and requires `protein_id`; metadata fallback modes are opt-in.

## `predMat` Rules

A valid `predMat` must have:

- phosphosite IDs as the index, such as `BTK;Y551;`
- kinase names as columns
- numeric scores in `[0, 1]`
- a unique, non-null index

`NaN` values are allowed for missing or unusable kinase scores. Infinite values are rejected.

## Checks You Will Hit Most Often

- required columns must exist
- required identifier columns must not be null
- numeric sample columns must remain numeric after coercion
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

## Important Behaviour

- by default, protein correction allows no silent phosphosite row loss
- site-matrix missing-value handling is explicit through `SiteMatrixPolicy.missing_data_policy`
- if the same phosphosite appears more than once after correction, PhosPy keeps the row with the highest mean corrected signal
- when prediction thresholds are too strict, prediction raises `NoCandidateKinasesError` instead of returning an empty invalid `predMat`

### `SiteMatrixPolicy.missing_data_policy`

- `drop_any_missing`: keep only complete corrected rows
- `retain_missing`: keep partial rows and preserve `NaN` values
- `require_min_observed_values`: keep rows meeting a minimum observed-value count

## Good Starting Points

Start at the boundary closest to your input:

- `PhosphoDataset.from_files(...)` for preprocessing from raw tables
- `SimpleKinaseWorkflow.run(...)` for preprocessing, prediction, and kinase activity in one step
- `KinaseActivityAnalyzer.run(...)` for analysis from an existing `predMat`
- `SignalomeWorkflow.run(...)` for downstream signalome construction

## Quick Troubleshooting

| Problem | Usually means | Good next step |
| --- | --- | --- |
| Missing required columns | Your input headers do not match the expected schema | Check cleaned column names and, if needed, pass a custom `DatasetSchema` |
| `predMat` overlap error | The phosphosite IDs do not line up between the matrix and `predMat` | Confirm both sides use `ENTITY;SITE;` IDs |
| `NoCandidateKinasesError` | Thresholds or inclusion settings filtered out every kinase | Relax the prediction config and rerun |
| Sequence coverage error | `site_sequences` or `motif_sequences` do not cover the scored sites | Check keys and confirm you passed the correct reference inputs |
| Signalome top-kinase failure | Your aligned prediction values contain non-finite rows | Clean or regenerate the prediction output before signalome construction |

## Recommended Beginner Path

1. Make sure your total and phospho tables use the required columns.
2. Run `SimpleKinaseWorkflow.run(...)`.
3. Inspect `result.pred_mat_result` and `result.kinase_activity_result`.
4. Export `predMat` with `result.pred_mat_result.to_csv("predMat.csv")` when needed.

Runnable example: [`../examples/simple_workflow_demo.py`](../examples/simple_workflow_demo.py)
