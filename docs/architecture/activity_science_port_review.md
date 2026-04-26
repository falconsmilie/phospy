# Activity Science Port Review (Legacy -> Rewrite)

- Date: 2026-04-20
- Scope:
  - Historical baseline science:
    - promoted from archived project snapshots in git history
  - Rewrite implementation:
    - `src/phospy/activities/scoring.py`
    - `src/phospy/activities/models.py`
    - `src/phospy/validation/workflows/activity.py`
    - `src/phospy/workflows/kinase/executor.py`

> Audience: maintainers auditing activity-science parity evidence.
> For first-time usage, use [Quickstart](../getting-started/quickstart-first-workflow.md).

This review is science-focused. Architectural differences (historical analyzer/writer
classes vs rewrite validator/interpreter/executor boundaries) are intentionally
excluded from gap classification.

## Scientific Verification Matrix

| Scientific behavior | Historical baseline | Rewrite status | Gap |
| --- | --- | --- | --- |
| Weighted kinase activity kernel | weighted mean over top-N predicted substrates, overlap-aware, sample-wise NaN handling | Matched in `_compute_weighted_kinase_activity` + `_nan_aware_weighted_average` | None |
| Thresholded substrate-mean activity kernel | thresholded substrate selection + per-sample mean with NaN skip | Matched in `_compute_thresholded_substrate_mean_activity` + `_nan_aware_mean_array` | None |
| Threshold semantics | strict `>` thresholding for thresholded-mean/target outputs | Matched in `_prediction_mask`, `_prediction_mask_array`, `_build_kinase_target_table` | None |
| Minimum substrate rules | candidate gating by `min_substrates` | Matched for both weighted and thresholded-mean paths | None |
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
  weighted/thresholded-mean candidates.

Historical donor terminology called this KSEA-style. In the rewrite public API
this is named `thresholded_substrate_mean_activity` because it computes a
thresholded substrate mean and does not implement full KSEA enrichment
statistics.

This is not a scientific divergence; it is explicit product-boundary behavior for
the supported workflow contract.

## Completion Decision

Activity science port is complete for the supported rewrite scope:

- Scientific kernels are baseline-aligned.
- Behavior is parity-backed on rewrite-owned fixtures.
- Active parity execution is rewrite-only (`tests/parity/test_activity_stage_parity.py`
  compares rewrite runtime outputs to committed fixtures and does not import or
  execute archived-runtime modules).
- Activity parity has an explicit hard CI regression gate
  (`activity-parity-gate`; marker selection `parity and activity_parity`).
- Fixture provenance and promotion policy are maintained in
  `tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md`.
- Activity remains part of `KinaseWorkflow` (no standalone public activity
  workflow).
- `activity_result` is a supported stage output when activity is enabled.

## Parity Independence and Cutover Status

- Rewrite-owned fixture expectations exist and are the active authority under
  `tests/fixtures/rewrite_parity/r_reference_l6/`.
- The active parity gate does not execute historical donor code;
  it compares rewrite runtime outputs to committed rewrite fixtures.
- Legacy donor references remain only as archival provenance and lock-style
  checks (for example, targeted donor-inventory/unit tests), not as the active
  parity execution path.
- Remaining cutover for supported activity behavior: none in this scoped pass.
  Future scientific additions should be tracked as new scoped gaps, not as
  unfinished activity-port debt.

## Where Next

- Broader parity tier policy: [Parity to PhosR](../parity.md)
- Full inventory context: [Legacy science gap audit](legacy_science_gap_audit.md)
- Maintainer entry: [Contributor and maintainer docs](../contributor/index.md)
