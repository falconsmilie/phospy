# R Reference Fixtures for Bundled PhosR L6 Data

This directory contains reference files based on PhosR’s bundled rat L6 myotube example dataset.

It serves two related purposes:

- stable R-side fixtures for downstream kinase-analysis parity checks on a more realistic dataset
- committed reference tables used by the newer native Python workflow parity layer

Generate the core L6 fixture set from the repository root with:

```bash
Rscript scripts/generate_r_l6_fixtures.R
```

You can also choose an explicit output directory:

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6
```

The generated downstream-analysis outputs include:

- `l6_phospho_matrix.csv`
- `l6_site_sequences.csv`
- `predMat.csv`
- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `sessionInfo.txt`

This directory also includes committed native-workflow reference tables used in parity tests for the Python-native
scoring and prediction seams:

- `native_substrate_map.csv`
- `native_profile_matrix.csv`
- `native_profile_scores.csv`
- `native_motif_scores.csv`
- `native_motif_sizes.csv`
- `native_combined_scores.csv`
- `native_combined_weights.csv`
- `native_candidate_substrates.csv`
- `native_prediction_top30.csv`

The `prediction_trace/` subdirectory contains R-side prediction debug traces for selected kinases. Those files are
useful for investigating where prediction-stage differences arise, but they should be treated as debugging aids rather
than blanket parity claims on their own.