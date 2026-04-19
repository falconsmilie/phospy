# Legacy Science Gap Audit: Rewrite vs Legacy Archive

- Date: 2026-04-19
- Purpose: governance truth source for scientific parity status.
- This revision explicitly separates historical gap labels from current
  implementation reality.

## Status Vocabulary (Normative)

- **Landed**: implemented in rewrite-native modules and covered by tests in the
  supported lane.
- **Open**: still missing in the supported rewrite lane.
- **Historical gap label**: previous ticket language retained for traceability,
  not proof that the gap is still open.

## Current Truth Snapshot (2026-04-19)

- The prior `SCI-GAP-01/05/06/08/09/10/11` items in this file are now landed in
  rewrite-native code paths.
- Prediction and signalome consume the same downstream score lane
  (`combined_scores` first, `profile_scores` fallback).
- Legacy modules remain scientific donors only; rewrite contract authority is
  `src/phospy/*` plus active ADR boundaries.

## Former Gap -> Current State Matrix

| Scientific area | Historical wording (former audit state) | Current rewrite truth | Rewrite-native evidence (code/tests) | Status now |
| --- | --- | --- | --- | --- |
| Adaptive ensemble prediction (`SCI-GAP-05`) | Deferred / not implemented | Adaptive ensemble mode is implemented and executed in kinase workflow (`mode="adaptive_ensemble"`). | `src/phospy/prediction/execution.py`, `src/phospy/prediction/sampling_core.py`, `src/phospy/workflows/kinase/executor.py`, `tests/unit/test_prediction_adaptive_sampling.py`, `tests/parity/test_adaptive_prediction_parity.py`, `tests/unit/test_prediction_mode_regression.py` | Landed |
| Clustering + module-count selection (`SCI-GAP-06`) | Deferred / dominant grouping only | Signalome clustering, module-count selection diagnostics, and protein module derivation are implemented and wired in executor flow. | `src/phospy/signalomes/clustering.py`, `src/phospy/workflows/signalome/executor.py`, `tests/unit/test_signalome_module_selection.py`, `tests/parity/test_signalome_workflow_parity.py` | Landed |
| Weighted-top assignment (`SCI-GAP-08`) | Metadata exists but policy not implemented | `assignment_policy="weighted_top"` is implemented for module shares and expanded signalome support attribution. | `src/phospy/api/configs.py`, `src/phospy/signalomes/assignments.py`, `src/phospy/signalomes/modules.py`, `tests/unit/test_signalome_science.py`, `tests/unit/test_signalome_bundle_compatibility.py` | Landed |
| Expanded signalome outputs (`SCI-GAP-10`) | Contract only, workflow returned `None` | Expanded signalome table is materialized in supported workflow output and carried through publishing/bundle paths. | `src/phospy/signalomes/expanded.py`, `src/phospy/workflows/signalome/executor.py`, `src/phospy/api/results.py`, `tests/parity/test_signalome_workflow_parity.py`, `tests/integration/test_signalome_workflow_integration.py`, `tests/integration/test_signalome_bundle_integration.py` | Landed |
| Signalome network policies (`SCI-GAP-09`) | Follow-on gap list still treated variants as open | `positive_only`, `absolute_threshold`, and `signed` policies are implemented, validated, and parity-tested. | `src/phospy/api/configs.py`, `src/phospy/signalomes/network.py`, `tests/unit/test_signalome_science.py`, `tests/parity/test_signalome_workflow_parity.py`, `tests/unit/test_validator_boundaries.py` | Landed |
| Profile missing-value strategy (`SCI-GAP-01`) | Partially ported (strict only) | Both `strict` and `median_skipna` are supported through public config and scoring behavior. | `src/phospy/api/configs.py`, `src/phospy/workflows/kinase/science.py`, `tests/unit/test_kinase_science.py`, `tests/parity/test_kinase_workflow_parity.py`, `tests/integration/test_kinase_workflow_integration.py` | Landed |
| Activity parity (`SCI-GAP-11`) | Kept as ongoing follow-on lock item | Activity weighted/KSEA kernels are implemented in supported lane and parity-locked against legacy fixtures. | `src/phospy/activities/scoring.py`, `src/phospy/workflows/kinase/executor.py`, `tests/parity/test_activity_stage_parity.py`, `docs/architecture/activity_science_port_review.md` | Landed (regression lock remains active) |

## Open Scientific Gaps (This Audit Scope)

- None confirmed in this scoped pass.
- If new scientific deltas appear, open new concrete gaps against rewrite-native
  seams rather than reviving superseded "deferred" wording for landed items.

## Historical Notes (Preserved)

- Earlier wording in this file that marked adaptive ensemble, clustering,
  weighted-top, expanded signalome, network policy variants, and
  `median_skipna` profile handling as deferred/partial is now superseded.
- Historical `SCI-GAP-*` labels are retained for traceability, but the items
  above are resolved in the current rewrite implementation.
