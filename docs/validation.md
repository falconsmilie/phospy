# Validation Quick Guide

This is the short version of what PhosPy checks.

For method signatures and parameter-by-parameter details, see [`api.md`](api.md).
For parity scope, see [`parity.md`](parity.md).

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
- `predMat` scores must stay in `[0, 1]`
- file paths must exist and point to files
- comparison pairs must use known schema groups and must not be duplicated
- downstream kinase analysis needs overlap between `predMat` and the phosphosite matrix
- native workflow runs need overlap across the matrix, substrate map, and sequence inputs
- motif-aware native workflow runs need both `motif_sequences` and `site_sequences`

## Useful Behaviour to Know

- by default, protein correction allows no silent phosphosite row loss
- site-matrix building can drop rows with missing sequence data or incomplete corrected values
- if the same phosphosite appears more than once after correction, PhosPy keeps the row with the highest mean corrected signal

## Good Starting Point

If you are unsure where validation happens, start here:

- `PhosphoDataset.from_files(...)` for the standard preprocessing path
- `KinaseActivityAnalyzer.run(...)` for analysis from an existing `predMat`
- `PhosRPipeline.from_files(...)` for the file-based one-shot flow
- `KinaseWorkflow.run(...)` for native prediction
