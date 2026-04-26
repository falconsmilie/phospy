# Synthetic adaptive-sampling decision seam fixtures

These fixtures are intentionally small and deterministic.
They complement the R-backed L6 parity traces by pinning edge-case replay behaviour:

- tied candidate scores with stable mergesort ordering
- tiny candidate and negative-pool sizes
- explicit per-iteration sampling overrides
- exact final top-site decisions on a deterministic replay path

Active parity-gate files in this directory are:

- `rank_weighted_fusion_scores.csv`
- `trace_candidates.csv`
- `trace_final_ensemble_predictions.csv`
- `trace_final_ensemble_top.csv`

Non-gated seam-debug trace tables were moved on 2026-04-22 to:
`tests/fixtures/archive/adaptive_sampling_edge_trace_debug/`.

They are not a standalone claim of PhosR parity.
