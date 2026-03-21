# R reference fixtures for bundled PhosR L6 data

This directory is reserved for richer CSV fixtures generated from the bundled PhosR L6 myotube dataset.

Generate them from the repository root with:

```bash
Rscript scripts/generate_r_l6_fixtures.R
```

Or choose an explicit output directory:

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6
```

Expected outputs:

- `l6_phospho_matrix.csv`
- `l6_site_sequences.csv`
- `predMat.csv`
- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `sessionInfo.txt`

These fixtures are intended for richer downstream parity checks that exercise the real PhosR kinase scoring and prediction workflow on the bundled L6 example dataset.
