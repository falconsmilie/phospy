# Python prediction trace fixtures

These files are generated from the native Python predictor with debug tracing enabled.
They are intended for direct comparison with the R trace fixtures generated from the PhosR L6 example path.

Trace kinases: MAPK9, IRAK1, TBK1, LCK
SVM mode: r_parity
Sampling trace override: none
Skipped trace kinases: none

Files:
- trace_candidates.csv: ranked combined-score candidates for the traced kinases
- trace_initial_negatives.csv: initial negative draw for each ensemble member
- trace_iteration_probabilities.csv: per-iteration class probabilities on the base train set
- trace_iteration_resampling_weights.csv: per-iteration class-specific resampling weights
- trace_iteration_samples.csv: resampled site identities for each iteration and class
- trace_final_ensemble_predictions.csv: final per-ensemble prediction probabilities for all sites
- trace_final_ensemble_top.csv: final per-ensemble top-ranked sites
