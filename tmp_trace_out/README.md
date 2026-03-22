# Temporary Python Prediction Trace Output

This directory exists as a working output area for ad hoc Python prediction-trace exports during debugging.

Unlike the committed fixture directories under `tests/fixtures/`, this folder is meant for targeted, short-lived trace
runs while investigating a specific kinase or learner-stage difference. Its contents may be regenerated, replaced, or
narrowed to a single debugging scenario.

The current contents reflect a focused Python trace run for:

- trace kinase: `MAPK1`
- per-ensemble top-N export: `5`
- sampling trace override: `tests/fixtures/r_reference_l6/prediction_trace`

The files are:

- `trace_candidates.csv`: ranked combined-score candidates for the traced kinase
- `trace_initial_negatives.csv`: initial negative draw for each ensemble member
- `trace_iteration_probabilities.csv`: per-iteration class probabilities on the base training set
- `trace_iteration_samples.csv`: resampled site identities for each iteration and class
- `trace_final_ensemble_predictions.csv`: final per-ensemble prediction probabilities for all sites
- `trace_final_ensemble_top.csv`: final per-ensemble top-ranked sites

Keep this directory when it is useful as a documented scratch area for debugging. Do not treat its contents as stable
parity fixtures unless they are intentionally promoted into a committed reference path under `tests/fixtures/`.