# R Reference Fixtures

This directory contains CSV fixtures generated from the synthetic example dataset by real R/PhosR code.

These files are used as stable reference outputs for the synthetic parity layer. They are most useful for deterministic
preprocessing checks, site-matrix construction checks, and regression protection around the currently implemented
downstream wrapper flow.

Generate the fixtures from the repository root with:

```bash
Rscript scripts/generate_r_fixtures.R
```

You can also choose explicit input and output paths:

```bash
Rscript scripts/generate_r_fixtures.R \
  --total examples/data/total.tsv \
  --phospho examples/data/phospho.tsv \
  --outdir tests/fixtures/r_reference
```

Expected outputs:

- `df_total_unique.csv`
- `df_total_filtered.csv`
- `df_phospho_filtered.csv`
- `df_phospho_corrected.csv`
- `phosr_input.csv`
- `mat_phospho_corrected.csv`
- `site_sequences.csv`
- `predMat.csv`
- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `sessionInfo.txt`

These fixtures should be treated as committed reference data. If behaviour changes intentionally, regenerate the
affected files and update the related parity tests in the same line of work.