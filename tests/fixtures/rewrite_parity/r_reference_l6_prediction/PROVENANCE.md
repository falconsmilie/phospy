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
  - `predMat.csv` (full prediction matrix surface)
  - `native_prediction_top30.csv` (candidate-restricted ranked export surface)
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
- candidate-set parity:
  - compare rewrite `substrate_list` kinase/site pairs with
    `native_candidate_substrates.csv`.
- prediction matrix ranking parity:
  - derive ranked sites directly from rewrite `pred_mat` and fixture
    `predMat.csv` (like-for-like matrix surface).
- top-k ranked export parity:
  - compare rewrite `substrate_list` ranked output with
    `native_prediction_top30.csv` (candidate-restricted export surface).
- `r_parity` policy checks are reported on both ranking surfaces against the
  promoted stable fixture lane and are interpreted as bounded divergence gates.

## Assertion Contract Layout

`tests/parity/test_l6_prediction_parity.py` is organized as separate contract
assertion groups so failure semantics are unambiguous:

- donor-vs-rewrite scoring table parity (`profile`, `combined`, `weights`)
- donor-vs-rewrite prediction-matrix numeric parity (`pred_mat` vs `predMat.csv`)
- donor-vs-rewrite prediction-matrix ranking parity (rank-order derived from
  `pred_mat` vs rank-order derived from `predMat.csv`)
- donor-vs-rewrite candidate-set parity
- donor-vs-rewrite ranked top-k export parity
- cross-policy divergence on prediction-matrix surface (`stable` vs `r_parity`)
- cross-policy divergence on ranked top-k export surface (`stable` vs
  `r_parity`)
- all cross-policy divergence checks are reported separately from donor parity

Release-gate thresholds for this lane are enforced in
`tests/parity/test_l6_prediction_parity.py`.

As of 2026-04-22, donor-vs-rewrite ranking hard gates in that file are:

- prediction-matrix and top-k export mean Spearman rank correlation `>= 0.96`
- prediction-matrix and top-k export mean top-20 overlap `>= 0.85`
- prediction-matrix and top-k export mean top-30 overlap `>= 0.88`
- prediction-matrix and top-k export good-top10-count `>= 20`

These values intentionally restore the legacy release-quality ranking bar after
the parity comparison surfaces were corrected to be like-for-like. Reporting
tests that emit parity summaries are informational diagnostics and are not used
as threshold-bearing release blockers.
