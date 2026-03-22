# Python Prediction Trace Fixtures

This directory contains Python-side prediction debug traces exported from the native predictor on the L6 reference path.

These files are intended to be compared directly against the R trace fixtures under
`tests/fixtures/r_reference_l6/prediction_trace/` so that prediction-stage differences can be inspected step by step.

Current trace settings:

- trace kinase: `PRKAA1`
- SVM mode: `r_parity`
- per-ensemble top-N export: `10`

The generated files are:

- `trace_candidates.csv`: ranked combined-score candidates for the traced kinase
- `trace_initial_negatives.csv`: initial negative draw for each ensemble member
- `trace_iteration_probabilities.csv`: per-iteration class probabilities on the base training set
- `trace_iteration_samples.csv`: resampled site identities for each iteration and class
- `trace_final_ensemble_predictions.csv`: final per-ensemble prediction probabilities for all sites
- `trace_final_ensemble_top.csv`: final per-ensemble top-ranked sites

You can regenerate these traces from the repository root with:

```bash
python scripts/export_python_prediction_traces.py \
  --trace-kinases PRKAA1 \
  --svm-mode r_parity \
  --debug-top-n 10 \
  --outdir tests/fixtures/python_reference_l6/prediction_trace
```

These files are committed Python reference traces for seam-level comparison. They are most useful when read alongside
the R trace fixtures and the parity notes in [`docs/parity.md`](../../../../docs/parity.md).