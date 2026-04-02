# Validation Quick Guide

This page is a short validation guide for common PhosPy inputs.

For parity to the R `PhosR` package, see [`docs/parity.md`](parity.md).
For method signatures and parameter-by-parameter validation, see [`docs/api.md`](api.md).

## File Types

When you load from files:

- total input is read as TSV
- phospho input is read as TSV
- `predMat` is read as CSV, using the first column as the phosphosite index

File-loaded total and phospho headers are normalised to lowercase snake case before validation.
Loading fails if two raw headers collapse to the same cleaned name.

## Required Input Shape

### Total table

Required columns:

- `genes`
- `group1` to `group6` by default, or the schema's configured total columns

### Phospho table

Required columns:

- `uid`
- `gene_names`
- `gene_p_site`
- `localization_prob`
- `centralized_sequence`
- `p_group1` to `p_group6` by default, or the schema's configured phospho columns

`gene_p_site` must split cleanly into a gene and site, such as `BTK_Y551`.

### `predMat`

`predMat` must be numeric and must have:

- phosphosite IDs as the index, such as `BTK;Y551;`
- kinase names as columns
- scores in the range `[0, 1]`

## Common Validation Rules

These are the checks most users hit first:

- total and phospho sample columns must be numeric
- `localization_prob` must stay in `[0, 1]`
- `predMat` values must stay in `[0, 1]`
- file paths must exist and point to files
- comparison pairs must use known schema groups and must not be duplicated
- downstream kinase analysis requires overlap between `predMat` and the phosphosite matrix
- native workflow runs require shared phosphosite IDs across the matrix, substrate map, and sequence inputs
- motif-aware native workflow runs require both `motif_sequences` and matching `site_sequences`

## Useful Behaviour to Know

- By default, protein correction does not allow silent phosphosite row loss.
- Site-matrix building drops rows with missing sequence information or incomplete corrected values.
- If the same phosphosite appears more than once after correction, PhosPy keeps the row with the highest mean corrected signal.

## Next Step

Use [`docs/api.md`](api.md) when you need the exact validation rules for a specific class or method.
