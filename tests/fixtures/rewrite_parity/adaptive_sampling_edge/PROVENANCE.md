# Rewrite Parity Fixture Provenance (`adaptive_sampling_edge`)

These fixtures are rewrite-owned donor copies used to keep adaptive
sampling/SVM legacy science scenarios visible while that lane is still deferred
from the supported rewrite contract.

## Source

Promoted from:

`tests_legacy/fixtures/synthetic_adaptive_sampling_edge/`

on 2026-04-19.

## Why This Directory Exists

- rewrite-side donor tests and inventory checks must not resolve their normal
  fixture paths from `tests_legacy/fixtures/`
- the deferred adaptive-sampling lane still needs explicit provenance and
  discoverable donor traces for follow-on science-gap work (`SCI-GAP-05`)

## Included Files

- `combined_scores.csv`
- `trace_candidates.csv`
- `trace_selected_candidates.csv`
- `trace_negative_pool.csv`
- `trace_initial_negatives.csv`
- `trace_iteration_labels.csv`
- `trace_iteration_decision_values.csv`
- `trace_iteration_probabilities.csv`
- `trace_iteration_probability_parameters.csv`
- `trace_iteration_resampling_weights.csv`
- `trace_iteration_samples.csv`
- `trace_final_ensemble_predictions.csv`
- `trace_final_ensemble_decision_values.csv`
- `trace_final_ensemble_top.csv`
- `README.md`

## Assertion Contract Layout

`tests/parity/test_adaptive_prediction_parity.py` keeps donor parity and
cross-policy divergence as separate contracts:

- donor-vs-rewrite candidate-selection parity (trace-aligned counts and policy
  invariance checks)
- donor-vs-rewrite prediction-matrix parity for `adaptive_policy="stable"`
- donor-vs-rewrite ranked-output parity for `adaptive_policy="stable"`
- donor-vs-rewrite prediction-matrix parity for `adaptive_policy="r_parity"`
- donor-vs-rewrite ranked-output parity for `adaptive_policy="r_parity"`
- cross-policy prediction-matrix divergence (`stable` vs `r_parity`)
- cross-policy ranked-output divergence (`stable` vs `r_parity`)
- cross-policy divergence is reported explicitly as policy comparison and is
  not interpreted as donor-port failure
