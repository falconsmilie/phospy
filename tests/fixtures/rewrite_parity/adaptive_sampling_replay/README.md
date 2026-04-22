# R prediction trace fixtures (active replay slice)

These files are generated from the bundled PhosR L6 example path.
They are intended for direct comparison with Python-side prediction debug traces.

Trace kinases: PRKAA1, MAPK1, MAPK9, IRAK1, TBK1, LCK
Per-ensemble top-N export: 10

Active files:
- trace_candidates.csv: ranked combined-score candidates for the traced kinases
- trace_initial_negatives.csv: initial negative draw for each ensemble member
- trace_iteration_samples.csv: resampled site identities for each iteration and class
- trace_final_ensemble_predictions.csv: final per-ensemble prediction probabilities for all sites
- trace_final_ensemble_top.csv: final per-ensemble top-ranked sites

Archived non-gated debug tables were moved on 2026-04-22 to:
`tests/fixtures/archive/adaptive_sampling_replay_trace_debug/`.
