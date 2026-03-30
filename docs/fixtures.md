# Fixtures and Traces

This page is the quick map for the fixture and trace directories used in PhosPy.

Unless noted otherwise, commands assume:

- the repo root
- Linux or macOS
- a shell that understands standard `bash` syntax

If you just want the usual shortcuts, the `Makefile` wraps the most common commands.

## `tests/fixtures/r_reference`

Small synthetic R-backed fixtures for deterministic preprocessing, site-matrix construction, and the current
downstream wrapper flow.

Generate them with either of these:

```bash
Rscript scripts/generate_r_fixtures.R
make fixtures-r-small
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

Generate them with either of these:

```bash
Rscript scripts/generate_r_l6_fixtures.R
make fixtures-r-l6
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

By default, it is derived from `tests/fixtures/r_reference_l6`.

Generate it with either of these:

```bash
python scripts/generate_fragile_support_reference.py
make fixtures-fragile
```

If you want to point at a different L6 source directory:

```bash
python scripts/generate_fragile_support_reference.py \
  --source-dir tests/fixtures/r_reference_l6 \
  --outdir tests/fixtures/fragile_support_reference
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

Generate them with either of these:

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6 \
  --trace_kinases PRKAA1,MAPK1 \
  --trace_top_n 10
make traces-r
```

These traces are especially useful when candidate selection agrees but the later learner path starts to drift.

## `tests/fixtures/python_reference_l6/prediction_trace`

Committed Python prediction traces used alongside the R traces.

Generate them with either of these:

```bash
python scripts/export_python_prediction_traces.py \
  --trace-kinases PRKAA1,MAPK1 \
  --svm-mode r_parity \
  --debug-top-n 10 \
  --outdir tests/fixtures/python_reference_l6/prediction_trace
make traces-python
```

Replay the R sampling rows in Python with either of these:

```bash
python scripts/export_python_prediction_traces.py \
  --trace-kinases PRKAA1,MAPK1 \
  --svm-mode r_parity \
  --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace \
  --outdir tests/fixtures/python_reference_l6/prediction_trace
make traces-python-replay
```

The replay mode is handy when you want the remaining delta to be model-side rather than sampling-side.

## `tests/fixtures/r_reference_l6_seam_stress`

A smaller L6-derived fixture family used to stress narrower prediction seams without carrying the whole L6 reference
surface.

Generate it with:

```bash
python scripts/generate_l6_seam_stress_reference.py --outdir tests/fixtures/r_reference_l6_seam_stress
make fixtures-r-l6-seam-stress
```

This directory includes a reduced prediction trace subtree together with the smaller score and prediction tables used by
`tests/test_reference_workflow_seams.py`. Helpful files to inspect first are `profile_scores.csv`, `combined_scores.csv`,
`prediction_top30.csv`, and the `prediction_trace/` directory.

## `tests/fixtures/synthetic_adaptive_sampling_edge`

A fully synthetic fixture family used to exercise deterministic adaptive-sampling edge cases.

Generate it with either of these:

```bash
python scripts/generate_synthetic_adaptive_sampling_edge_fixtures.py --outdir tests/fixtures/synthetic_adaptive_sampling_edge
make fixtures-synthetic-edge
```

Key files include:

```text
combined_scores.csv
trace_candidates.csv
trace_initial_negatives.csv
trace_iteration_samples.csv
trace_iteration_probabilities.csv
trace_final_ensemble_predictions.csv
trace_final_ensemble_top.csv
README.md
```

## `tmp_trace_out`

Scratch output for temporary debugging runs. This directory is not part of the committed parity contract unless you
deliberately promote something into a fixture path.

Typical use:

```bash
python scripts/export_python_prediction_traces.py \
  --trace-kinases MAPK1 \
  --svm-mode r_parity \
  --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace \
  --outdir tmp_trace_out
```

## Rebuild Everything From Scratch

If you want the full local fixture bootstrap path, use:

```bash
make fixtures-all
```

That wraps the currently committed fixture families and trace exports in the expected order.

# Promotion Rule

If temporary debugging output reveals a stable and useful new seam:

1. promote the relevant files into a committed fixture directory
2. update the related tests
3. update this page or [`docs/parity.md`](parity.md) if the documented scope changes
