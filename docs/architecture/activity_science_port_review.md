# Activity Science Port Review (Legacy -> Rewrite)

- Date: 2026-04-18
- Scope:
  - Legacy donor science:
    - `legacy_archive/phospy_legacy/activities/scoring.py`
    - `legacy_archive/phospy_legacy/activities/analysis.py`
    - `legacy_archive/phospy_legacy/activities/results.py`
  - Rewrite implementation:
    - `src/phospy/activities/scoring.py`
    - `src/phospy/activities/models.py`
    - `src/phospy/validation/workflows/activity.py`
    - `src/phospy/workflows/kinase/executor.py`

This review is science-focused. Architectural differences (legacy analyzer/writer
classes vs rewrite validator/interpreter/executor boundaries) are intentionally
excluded from gap classification.

## Scientific Verification Matrix

| Scientific behavior | Legacy reference | Rewrite status | Gap |
| --- | --- | --- | --- |
| Weighted kinase activity kernel | weighted mean over top-N predicted substrates, overlap-aware, sample-wise NaN handling | Matched in `_compute_weighted_kinase_activity` + `_nan_aware_weighted_average` | None |
| KSEA-style score kernel | thresholded substrate selection + per-sample mean with NaN skip | Matched in `_compute_ksea_scores` + `_nan_aware_mean_array` | None |
| Threshold semantics | strict `>` thresholding for KSEA/target outputs | Matched in `_prediction_mask`, `_prediction_mask_array`, `_build_kinase_target_table` | None |
| Minimum substrate rules | candidate gating by `min_substrates` | Matched for both weighted and KSEA paths | None |
| pred_mat/phospho overlap policy | require overlap count/fraction floors | Matched in `KinaseActivityInputValidator` (same default floors) | None |
| Missing-value handling per sample | weighted and mean kernels ignore NaN per sample | Matched; covered in unit and parity tests | None |
| Top-N deterministic tie handling | deterministic ordering for tied prediction scores | Matched via stable argsort (`kind="stable"`) | None |
| Target count behavior | thresholded count per kinase, sorted descending | Matched in `_count_predicted_targets` | None |
| Target edge-table derivation | thresholded long-form edges sorted by kinase/score | Matched in `_build_kinase_target_table` | None |

## Rewrite-Specific Boundary Policy (Intentional)

Rewrite keeps legacy science kernels but adds workflow-boundary diagnostics:

- `seam=kinase.activity.input_overlap` for invalid overlap between prediction and
  phospho matrices.
- `seam=kinase.activity.valid_candidates` when activity filters remove all
  weighted/KSEA candidates.

This is not a scientific divergence; it is explicit product-boundary behavior for
the supported workflow contract.

## Completion Decision

Activity science port is complete for the supported rewrite scope:

- Scientific kernels are legacy-aligned.
- Behavior is parity-backed on rewrite-owned fixtures.
- Activity remains part of `KinaseWorkflow` (no standalone public activity
  workflow).
- `activity_result` is a supported stage output when activity is enabled.
