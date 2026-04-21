# Rewrite Parity Fixture Provenance (`r_reference_l6_prediction`)

This directory contains rewrite-owned parity references for the promoted L6
downstream scoring/prediction lane.

## Source

Initial promotion source (2026-04-20):

- `tests_legacy/fixtures/r_reference_l6/native_profile_scores.csv`
- `tests_legacy/fixtures/r_reference_l6/native_combined_scores.csv`
- `tests_legacy/fixtures/r_reference_l6/native_combined_weights.csv`
- `tests_legacy/fixtures/r_reference_l6/native_candidate_substrates.csv`
- `tests_legacy/fixtures/r_reference_l6/native_prediction_top30.csv`
- `tests_legacy/fixtures/r_reference_l6/predMat.csv`
- `tests_legacy/fixtures/r_reference_l6/l6_phospho_matrix.csv`

Rewrite refresh source (2026-04-21):

- generated from supported rewrite workflow execution:
  - `dataset`: rewrite-owned rat L6 dataset fixture (`build_rat_l6_dataset(n_sites=None)`)
  - `references`: `ReferencePreset.AUTO` (resolved rat bundled references)
  - `scoring_config`: `min_substrates=2`, `include_diagnostic_scoring_tables=True`
  - `prediction_config`:
    `mode="adaptive_ensemble"`, `adaptive_policy="stable"`, `top_k=30`,
    `ensemble_size=10`, `n_iterations=5`, `random_state=1`
- refreshed tables:
  - `native_combined_scores.csv`
  - `native_combined_weights.csv`
  - `native_candidate_substrates.csv`
  - `native_prediction_top30.csv`
  - `predMat.csv`
- `native_profile_scores.csv` remains promoted from the donor lane and is still
  parity-locked separately on the shared surface.

## Generation Notes

- 2026-04-20 promotion copied donor files without modification.
- 2026-04-21 refresh intentionally re-promoted stable-lane rewrite outputs to
  remove donor-fit scoring assumptions from the supported scoring contract.
- Active parity tests consume this directory via
  `tests/support/rewrite_fixture_data.py`.

## Comparison Policy Notes

The active rewrite parity test family uses explicit per-surface policies:

- profile scores: strict shared-surface numeric parity (very tight tolerance).
- combined scores / weights: strict parity against promoted rewrite references.
- candidate substrates: strict overlap against promoted rewrite references.
- prediction rankings:
  - stable lane: strict ranked-reference parity floors against promoted rewrite
    reference tables.
  - `r_parity` lane: bounded divergence floors versus the promoted stable
    reference lane.

Release-gate thresholds for this lane are enforced in
`tests/parity/test_l6_prediction_parity.py`.
