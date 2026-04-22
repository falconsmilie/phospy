# Rewrite Parity Fixture Provenance (`r_reference_l6_prediction`)

This directory contains rewrite-owned parity references for the promoted L6
scoring/prediction lane.

## Source

- Initial promotion from historical project snapshots (2026-04-20).
- Rewrite refresh from supported workflow execution (2026-04-21).

## Comparison Policy Notes

The active parity family uses explicit per-surface policies:

- profile-scores parity
- combined-scores/weights parity
- candidate-set parity
- prediction-matrix ranking parity
- top-k ranked-export parity
- cross-policy divergence checks (`stable` vs `r_parity`) reported separately

Release-gate thresholds for this lane are enforced in
`tests/parity/test_l6_prediction_parity.py`.
