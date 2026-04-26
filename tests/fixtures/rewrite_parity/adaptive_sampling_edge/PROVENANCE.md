# Rewrite Parity Fixture Provenance (`adaptive_sampling_edge`)

These fixtures are rewrite-owned historical-baseline copies used by active
adaptive prediction parity tests.

## Source

- Promoted from historical project snapshots on 2026-04-19.

## Included Files

- `rank_weighted_fusion_scores.csv`
- `trace_candidates.csv`
- `trace_final_ensemble_predictions.csv`
- `trace_final_ensemble_top.csv`
- `README.md`

## Assertion Contract Layout

`tests/parity/test_adaptive_prediction_parity.py` keeps baseline parity and
cross-policy divergence as separate contracts.
