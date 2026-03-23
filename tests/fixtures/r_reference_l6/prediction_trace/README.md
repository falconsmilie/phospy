# R prediction trace fixtures

These files are generated from the bundled PhosR L6 example path.
They are intended for direct comparison with Python-side prediction debug traces.

Trace kinases: PRKAA1, MAPK1
Per-ensemble top-N export: 10

Files:
- trace_candidates.csv: ranked combined-score candidates for the traced kinases
- trace_initial_negatives.csv: initial negative draw for each ensemble member
- trace_iteration_probabilities.csv: per-iteration class probabilities on the base train set
- trace_iteration_resampling_weights.csv: per-iteration class-specific resampling weights
- trace_iteration_samples.csv: resampled site identities for each iteration and class
- trace_final_ensemble_predictions.csv: final per-ensemble prediction probabilities for all sites
- trace_final_ensemble_top.csv: final per-ensemble top-ranked sites
