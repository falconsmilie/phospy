# R prediction trace fixtures

These files are generated from the bundled PhosR L6 example path.
They are intended for direct comparison with Python-side prediction debug traces.

Trace kinases: PRKAA1, MAPK1, MAPK9, IRAK1, TBK1, LCK
Per-ensemble top-N export: 10

Files:
- trace_candidates.csv: ranked combined-score candidates for the traced kinases
- trace_selected_candidates.csv: selected candidate substrate site IDs per traced kinase
- trace_negative_pool.csv: full negative-pool site IDs available to each traced kinase
- trace_initial_negatives.csv: initial negative draw for each ensemble member
- trace_iteration_labels.csv: per-iteration base-train labels before resampling
- trace_iteration_probabilities.csv: per-iteration class probabilities on the base train set
- trace_iteration_probability_parameters.csv: per-iteration libsvm probability-calibration parameters
- trace_iteration_decision_values.csv: per-iteration binary decision values aligned to class 1
- trace_iteration_resampling_weights.csv: per-iteration class-specific resampling weights
- trace_iteration_samples.csv: resampled site identities for each iteration and class
- trace_final_ensemble_predictions.csv: final per-ensemble prediction probabilities for all sites
- trace_final_ensemble_decision_values.csv: final per-ensemble binary decision values aligned to class 1
- trace_final_ensemble_top.csv: final per-ensemble top-ranked sites
