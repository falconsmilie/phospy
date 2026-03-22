# tests/fixtures/r_reference_l6/prediction_trace/README.md

# R Prediction Trace Fixtures

This directory contains R-side prediction debug traces generated from the bundled PhosR L6 example path.

These files are intended for direct comparison with Python-side trace exports so that learner-stage differences can be
inspected more concretely, rather than inferred from final ranking gaps alone.

Current trace settings:

- trace kinases: `PRKAA1`, `MAPK1`
- per-ensemble top-N export: `10`

The generated files are:

- `trace_candidates.csv`: ranked combined-score candidates for the traced kinases
- `trace_initial_negatives.csv`: initial negative draw for each ensemble member
- `trace_iteration_probabilities.csv`: per-iteration class probabilities on the base training set
- `trace_iteration_samples.csv`: resampled site identities for each iteration and class
- `trace_final_ensemble_predictions.csv`: final per-ensemble prediction probabilities for all sites
- `trace_final_ensemble_top.csv`: final per-ensemble top-ranked sites

Regenerate these traces from the repository root with:

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6 \
  --trace_kinases PRKAA1,MAPK1 \
  --trace_top_n 10
```

These files are committed reference traces for debugging and seam-level comparison. They are useful evidence for
specific learner and sampling steps, but they are not a standalone claim of full workflow parity.