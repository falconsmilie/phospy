# Fixture and Trace Directory Guide

This document is the single source of truth for the fixture and trace directories used in PhosPy.

Unless otherwise noted, commands below assume:

- **Linux**
- **repo root**
- a shell that understands standard `bash` syntax

macOS uses the same commands unless a section says otherwise. Windows is only shown where the syntax changes.

## Directory Map

### Small Synthetic R Fixtures

**Directory**

```text
tests/fixtures/r_reference
```

**Purpose**

- deterministic preprocessing parity checks
- site-matrix construction parity checks
- regression protection around the current downstream wrapper flow

**Generate**

```bash
Rscript scripts/generate_r_fixtures.R
```

**Optional Explicit Paths**

```bash
Rscript scripts/generate_r_fixtures.R \
  --total examples/data/total.tsv \
  --phospho examples/data/phospho.tsv \
  --outdir tests/fixtures/r_reference
```

**Key Outputs**

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
sessionInfo.txt
```

**Contract Status**

- committed reference data
- part of the parity contract

### Bundled PhosR L6 R Fixtures

**Directory**

```text
tests/fixtures/r_reference_l6
```

**Purpose**

- downstream kinase-analysis parity checks on a more realistic dataset
- committed R-side reference tables for native workflow seams

**Generate**

```bash
Rscript scripts/generate_r_l6_fixtures.R
```

**Optional Explicit Output Path**

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6
```

**Key Outputs**

```text
l6_phospho_matrix.csv
l6_site_sequences.csv
predMat.csv
kinase_activity_matrix.csv
ksea_scores.csv
ksea_counts.csv
kinase_target_counts.csv
sessionInfo.txt
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

**Contract Status**

- committed reference data
- part of the parity contract

### Committed R Prediction Traces

**Directory**

```text
tests/fixtures/r_reference_l6/prediction_trace
```

**Purpose**

- seam-level debugging of the adaptive-sampling decision stage
- direct comparison against Python trace exports

**Generate**

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6 \
  --trace_kinases PRKAA1,MAPK1 \
  --trace_top_n 10
```

**Key Outputs**

```text
trace_candidates.csv
trace_selected_candidates.csv
trace_negative_pool.csv
trace_initial_negatives.csv
trace_iteration_labels.csv
trace_iteration_probabilities.csv
trace_iteration_decision_values.csv
trace_iteration_samples.csv
trace_final_ensemble_decision_values.csv
trace_final_ensemble_predictions.csv
trace_final_ensemble_top.csv
```

**Contract Status**

- committed reference traces
- part of the adaptive-sampling decision seam parity story
- not a standalone claim of full workflow parity

### Committed Python Prediction Traces

**Directory**

```text
tests/fixtures/python_reference_l6/prediction_trace
```

**Purpose**

- direct comparison against the committed R prediction traces
- seam-level inspection of adaptive-sampling learner-stage differences

**Generate**

```bash
python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --debug-top-n 10 --outdir tests/fixtures/python_reference_l6/prediction_trace
```

**Replay the R Sampling Rows**

```bash
python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace --outdir tests/fixtures/python_reference_l6/prediction_trace
```

**Key Outputs**

```text
trace_candidates.csv
trace_selected_candidates.csv
trace_negative_pool.csv
trace_initial_negatives.csv
trace_iteration_labels.csv
trace_iteration_probabilities.csv
trace_iteration_decision_values.csv
trace_iteration_samples.csv
trace_final_ensemble_decision_values.csv
trace_final_ensemble_predictions.csv
trace_final_ensemble_top.csv
```

**Contract Status**

- committed Python reference traces
- part of the adaptive-sampling decision seam parity story

### Synthetic Adaptive-Sampling Edge Fixtures

**Directory**

```text
tests/fixtures/synthetic_adaptive_sampling_edge
```

**Purpose**

- small deterministic replay fixture family for adaptive-sampling edge cases
- exercises tied candidate scores, tiny negative pools, and explicit replayed sampling rows
- complements the R-backed L6 seam without making a broader parity claim

**Generate**

```bash
PYTHONPATH=src python scripts/generate_synthetic_adaptive_sampling_edge_fixtures.py
```

**Key Outputs**

```text
combined_scores.csv
trace_candidates.csv
trace_selected_candidates.csv
trace_negative_pool.csv
trace_initial_negatives.csv
trace_iteration_labels.csv
trace_iteration_samples.csv
trace_final_ensemble_top.csv
```

**Contract Status**

- committed synthetic regression fixtures
- not part of the R parity contract
- intended to pin replay behaviour around adaptive-sampling edge cases

### Temporary Python Trace Output

**Directory**

```text
tmp_trace_out
```

**Purpose**

- short-lived debugging output
- investigation-specific trace runs
- scratch space before promotion into committed fixture paths

**Typical Use**

```bash
python scripts/export_python_prediction_traces.py --trace-kinases MAPK1 --svm-mode r_parity --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace --outdir tmp_trace_out
```

**Contract Status**

- not committed reference data by default
- not part of the parity contract unless deliberately promoted

## Promotion Rule

When temporary debugging output reveals an important new seam or a stable comparison case:

1. promote the relevant outputs into an appropriate committed fixture directory
2. update the related tests
3. update this file or [`docs/parity.md`](parity.md) if the documented scope has changed