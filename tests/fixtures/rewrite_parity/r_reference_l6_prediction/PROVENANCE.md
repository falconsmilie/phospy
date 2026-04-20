# Rewrite Parity Fixture Provenance (`r_reference_l6_prediction`)

This directory contains rewrite-owned parity references for the promoted L6
downstream scoring/prediction lane.

## Source

Promoted from archival donor fixtures:

- `tests_legacy/fixtures/r_reference_l6/native_profile_scores.csv`
- `tests_legacy/fixtures/r_reference_l6/native_combined_scores.csv`
- `tests_legacy/fixtures/r_reference_l6/native_combined_weights.csv`
- `tests_legacy/fixtures/r_reference_l6/native_candidate_substrates.csv`
- `tests_legacy/fixtures/r_reference_l6/native_prediction_top30.csv`
- `tests_legacy/fixtures/r_reference_l6/predMat.csv`
- `tests_legacy/fixtures/r_reference_l6/l6_phospho_matrix.csv`

Promotion date: 2026-04-20.

## Generation Notes

- Files were copied without modification into rewrite-owned fixture paths.
- Active parity tests consume this directory via
  `tests/support/rewrite_fixture_data.py`.

## Comparison Policy Notes

The active rewrite parity test family uses explicit per-surface policies:

- profile scores: strict shared-surface numeric parity (very tight tolerance).
- combined scores / weights: explicit numeric tolerance + correlation floors.
- candidate substrates: overlap metrics (precision/recall/F1) and kinase-identity
  checks.
- prediction rankings: rank-correlation + top-10/20/30 overlap summaries for
  `adaptive_policy="stable"` and `adaptive_policy="r_parity"`.
