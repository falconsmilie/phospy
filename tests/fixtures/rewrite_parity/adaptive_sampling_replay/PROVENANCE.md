# Rewrite Parity Fixture Provenance (`adaptive_sampling_replay`)

This directory contains rewrite-owned replay-trace baseline material used by
active parity tests.

## Source

- Promoted from historical project snapshots on 2026-04-20.

## Included Files

- `combined_scores.csv`
- `trace_candidates.csv`
- `trace_initial_negatives.csv`
- `trace_iteration_samples.csv`
- `trace_final_ensemble_predictions.csv`
- `trace_final_ensemble_top.csv`
- `README.md`

## Comparison Policy Notes

Active replay parity assertions include:

- initial negative-pool overlap
- per-iteration sample-membership overlap
- final ensemble probability correlation/MAE
- top-rank and top-set overlap metrics
- deterministic replay checks under fixed seed for both adaptive policies
