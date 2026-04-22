# Rewrite Parity Fixture Provenance (`adaptive_sampling_replay`)

This directory contains rewrite-owned adaptive replay-trace donor material used
by active parity tests.

## Source

Promoted from archival donor fixtures:

- `tests_legacy/fixtures/r_reference_l6/native_combined_scores.csv`
- `tests_legacy/fixtures/r_reference_l6/prediction_trace/*`

Promotion date: 2026-04-20.

## Included Files

- `combined_scores.csv`
- `trace_candidates.csv`
- `trace_initial_negatives.csv`
- `trace_iteration_samples.csv`
- `trace_final_ensemble_predictions.csv`
- `trace_final_ensemble_top.csv`
- `README.md`

## Archived Debug Tables

On 2026-04-22, non-gated seam-debug trace tables were moved to:

- `tests/fixtures/archive/adaptive_sampling_replay_trace_debug/`

These files are retained as historical provenance only and are not part of the
active replay parity gate lane.

## Generation Notes

- Files were promoted as donor replay references; active tests execute rewrite
  adaptive sampling and compare replay surfaces against these references.

## Comparison Policy Notes

Active replay parity assertions include:

- initial negative-pool overlap (set-based replay surface)
- per-iteration sample-membership overlap
- final ensemble probability correlation/MAE
- top-rank and top-set overlap metrics
- deterministic replay checks under fixed seed for both
  `adaptive_policy="stable"` and `adaptive_policy="r_parity"`

This fixture lane is used as a hard rewrite parity gate in
`tests/parity/test_adaptive_replay_parity.py` with donor-level overlap and
correlation thresholds.
