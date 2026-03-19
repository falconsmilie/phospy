# R reference fixtures

This directory is reserved for CSV fixtures generated from the synthetic example dataset by real R/PhosR code.

Generate them from the repository root with:

```bash
Rscript scripts/generate_r_fixtures.R
```

Or choose explicit paths:

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

These files are intended to become stable parity fixtures for Python-vs-R tests.
