# docs/fixtures.md

# Fixture and Trace Directory Guide

This document explains the fixture and trace directories used in PhosPy, what each one is for, how the files are
generated, and which outputs are treated as committed reference data.

## Overview

The repository contains three kinds of reference or debugging artefacts:

- committed R reference fixtures
- committed Python reference traces
- temporary debugging output

These do not all serve the same purpose, and they should not all be treated as part of the parity contract.

## Committed Fixture and Trace Directories

### `tests/fixtures/r_reference/`

This directory contains committed CSV fixtures generated from the small synthetic example dataset by real R/PhosR code.

Use this path for:

- deterministic preprocessing parity checks
- site-matrix construction parity checks
- regression protection around the currently implemented downstream wrapper flow

Generate these fixtures from the repository root with:

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

Expected outputs include:

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

These files are committed reference data.

### `tests/fixtures/r_reference_l6/`

This directory contains committed reference files based on PhosR’s bundled rat L6 myotube example dataset.

Use this path for:

- downstream kinase-analysis parity checks on a more realistic dataset
- committed R-side reference tables for the native workflow parity layer

Generate the core L6 fixture set from the repository root with:

```bash
Rscript scripts/generate_r_l6_fixtures.R
```

You can also choose an explicit output directory:

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6
```

Committed downstream-analysis outputs include:

- `l6_phospho_matrix.csv`
- `l6_site_sequences.csv`
- `predMat.csv`
- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `sessionInfo.txt`

This directory also includes committed native-workflow reference tables used in parity tests for Python-native scoring
and prediction seams:

- `native_substrate_map.csv`
- `native_profile_matrix.csv`
- `native_profile_scores.csv`
- `native_motif_scores.csv`
- `native_motif_sizes.csv`
- `native_combined_scores.csv`
- `native_combined_weights.csv`
- `native_candidate_substrates.csv`
- `native_prediction_top30.csv`

These files are committed reference data.

### `tests/fixtures/r_reference_l6/prediction_trace/`

This directory contains committed R-side prediction traces generated from the bundled L6 path.

Use this path for seam-level debugging of the prediction stage, including:

- candidate ranking
- initial negative sampling
- iteration-level class probabilities
- iteration-level sample identities
- final ensemble predictions
- top-ranked final outputs

Generate these traces from the repository root with:

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6 \
  --trace_kinases PRKAA1,MAPK1 \
  --trace_top_n 10
```

Typical outputs include:

- `trace_candidates.csv`
- `trace_initial_negatives.csv`
- `trace_iteration_probabilities.csv`
- `trace_iteration_samples.csv`
- `trace_final_ensemble_predictions.csv`
- `trace_final_ensemble_top.csv`

These files are committed reference traces. They are useful seam-level evidence, but they should not be treated as a
standalone claim of full workflow parity.

### `tests/fixtures/python_reference_l6/prediction_trace/`

This directory contains committed Python-side prediction traces exported from the native predictor on the L6 reference
path.

Use this path for direct comparison against the committed R-side prediction traces.

Generate these traces from the repository root with:

```bash
python scripts/export_python_prediction_traces.py \
  --trace-kinases PRKAA1 \
  --svm-mode r_parity \
  --debug-top-n 10 \
  --outdir tests/fixtures/python_reference_l6/prediction_trace
```

Typical outputs include:

- `trace_candidates.csv`
- `trace_initial_negatives.csv`
- `trace_iteration_probabilities.csv`
- `trace_iteration_samples.csv`
- `trace_final_ensemble_predictions.csv`
- `trace_final_ensemble_top.csv`

These files are committed Python reference traces for seam-level comparison.

## Temporary Debugging Output

### `tmp_trace_out/`

This directory exists as a documented scratch area for ad hoc Python prediction-trace exports during debugging.

Use this path when you want to investigate a specific kinase or learner-stage difference without immediately promoting
the results into a committed fixture directory.

Its purpose is different from the committed fixture paths:

- it is for short-lived or investigation-specific trace runs
- its contents may be regenerated or replaced freely
- it should not be treated as stable reference data
- files only become part of the parity contract when they are intentionally promoted into a committed fixture path under
  `tests/fixtures/`

Keep this directory if it is useful as a working area. Do not treat its contents as authoritative evidence unless they
are deliberately promoted.

## Which Files Are Part of the Parity Contract?

Treat the following as part of the committed parity story:

- files under `tests/fixtures/r_reference/`
- files under `tests/fixtures/r_reference_l6/`
- files under `tests/fixtures/r_reference_l6/prediction_trace/`
- files under `tests/fixtures/python_reference_l6/prediction_trace/`

Do **not** treat files under `tmp_trace_out/` as part of the parity contract by default.

## Promotion Rule

When temporary debugging output reveals an important new seam or a stable comparison case:

1. promote the relevant outputs into an appropriate committed fixture directory
2. update the related tests
3. update `docs/parity.md` or this file if the documented scope has changed
