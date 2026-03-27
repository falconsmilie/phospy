# Fixtures and Traces

This page is the quick map for the fixture and trace directories used in PhosPy.

Unless noted otherwise, commands assume:

- the repo root
- Linux or macOS
- a shell that understands standard `bash` syntax

## `tests/fixtures/r_reference`

Small synthetic R-backed fixtures for deterministic preprocessing, site-matrix construction, and the current
downstream wrapper flow.

Generate them with:

```bash
Rscript scripts/generate_r_fixtures.R
```

Key files include:

```text
df_total_unique.csv
df_total_filtered.csv
df_phospho_filtered.csv
df_phospho_corrected.csv
phosr_input.csv
mat_phospho_corrected.csv
site_sequences.csv
predMat.csv
kinase_activity_matrix.csv
ksea_scores.csv
ksea_counts.csv
kinase_target_counts.csv
```

## `tests/fixtures/r_reference_l6`

Committed R-backed references for the bundled L6 dataset. This directory supports downstream kinase-analysis parity and
selected native workflow seams.

Generate them with:

```bash
Rscript scripts/generate_r_l6_fixtures.R
```

Key files include:

```text
l6_phospho_matrix.csv
l6_site_sequences.csv
predMat.csv
kinase_activity_matrix.csv
ksea_scores.csv
ksea_counts.csv
kinase_target_counts.csv
native_substrate_map.csv
native_profile_matrix.csv
native_profile_scores.csv
native_motif_scores.csv
native_motif_sizes.csv
native_combined_scores.csv
native_combined_weights.csv
native_candidate_substrates.csv
native_prediction_top30.csv
```

## `tests/fixtures/fragile_support_reference`

A curated reference dataset used to widen evidence beyond the main L6 path and to stress native-workflow decision
boundaries.

Generate it with:

```bash
python scripts/generate_fragile_support_reference.py
```

Key files include:

```text
phospho_matrix.csv
site_sequences.csv
substrate_map.csv
motif_sequences.csv
profile_matrix.csv
profile_sizes.csv
profile_scores.csv
motif_scores.csv
motif_sizes.csv
combined_scores.csv
combined_weights.csv
candidate_substrates.csv
screening_summary.csv
README.md
```

## `tests/fixtures/r_reference_l6/prediction_trace`

Committed R prediction traces for seam-level debugging of the prediction stage.

Generate them with:

```bash
Rscript scripts/generate_r_l6_fixtures.R --outdir tests/fixtures/r_reference_l6 --trace_kinases PRKAA1,MAPK1 --trace_top_n 10
```

## `tests/fixtures/python_reference_l6/prediction_trace`

Committed Python prediction traces used alongside the R traces.

Generate them with:

```bash
python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --debug-top-n 10 --outdir tests/fixtures/python_reference_l6/prediction_trace
```

Replay the R sampling rows in Python with:

```bash
python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace --outdir tests/fixtures/python_reference_l6/prediction_trace
```

## `tmp_trace_out`

Scratch output for temporary debugging runs. This directory is not part of the committed parity contract unless you
deliberately promote something into a fixture path.

Typical use:

```bash
python scripts/export_python_prediction_traces.py --trace-kinases MAPK1 --svm-mode r_parity --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace --outdir tmp_trace_out
```

# Promotion Rule

If temporary debugging output reveals a stable and useful new seam:

1. promote the relevant files into a committed fixture directory
2. update the related tests
3. update this page or [`docs/parity.md`](parity.md) if the documented scope changes
