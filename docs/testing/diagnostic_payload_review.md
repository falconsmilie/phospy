# Diagnostic Payload Review (TST-DIAG-002)

This document records human review decisions for diagnostic payload test candidates from
`docs/testing/diagnostic_payload_test_report.md` before any consolidation or rewrites.

## Decision Vocabulary

- `keep-public-contract`
- `keep-scientific-reproducibility`
- `rewrite-to-public-behaviour`
- `narrow-assertion-scope`
- `keep-internal-with-rationale`
- `defer`

## Recommended Action Vocabulary

- `keep`
- `rewrite`
- `narrow`
- `defer`

## Candidate Decisions

| Test Path | Test Name / Candidate Group | Current Classification (report) | Review Decision | Protected Contract (if any) | Recommended Action | Follow-up Ticket ID | Notes | Implementation Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `tests/unit/test_signalome_workflow_diagnostics.py` | `test_executor_uses_preconditioned_scores_when_missing_rows_are_present` | mixed/needs review | `narrow-assertion-scope` | Workflow diagnostic contract for preconditioning status and missing-row accounting exposed in boundary diagnostics | `narrow` | `TST-DIAG-003` | Keep stable workflow-facing fields; trim incidental table-shape and implementation-coupled payload details. | Completed in `TST-DIAG-003` on 2026-05-09. |
| `tests/unit/test_kinase_workflow_diagnostics.py` | `test_activity_stage_returns_weighted_thresholded_mean_and_target_outputs` | mixed/needs review | `keep-scientific-reproducibility` | Scientific/reproducibility contract for activity-stage diagnostic outputs and threshold summary behavior | `keep` | `TST-DIAG-002` | Dense assertions are acceptable here where they guard reproducible scientific diagnostics; avoid expanding to new incidental fields. | Out of scope in `TST-DIAG-003` (kept by prior review decision). |
| `tests/unit/test_kinase_workflow_diagnostics.py` | `test_boundary_error_reports_empty_eligible_kinase_set_counts` | mixed/needs review | `keep-public-contract` | Public boundary-error diagnostic contract: empty eligible kinase-set counts and failure explanation are part of workflow error surface | `keep` | `TST-DIAG-002` | Preserve key count fields and message contract; if fields expand later, treat as additive-only. | Out of scope in `TST-DIAG-003` (kept by prior review decision). |
| `tests/unit/test_kinase_workflow_diagnostics.py` | `test_boundary_error_reports_unusable_reference_coverage_counts` | mixed/needs review | `keep-public-contract` | Public boundary-error diagnostic contract for unusable reference coverage counts | `keep` | `TST-DIAG-002` | Keep count/coverage semantics stable; avoid locking unrelated formatting details. | Out of scope in `TST-DIAG-003` (kept by prior review decision). |
| `tests/unit/test_signalome_workflow_diagnostics.py` | `test_boundary_error_reports_network_failure_modes` | mixed/needs review | `narrow-assertion-scope` | Workflow failure-mode contract for network diagnostics and reason signaling | `narrow` | `TST-DIAG-003` | Retain stable failure-mode keys and reason classes; narrow brittle key-presence breadth to contract-relevant fields. | Completed in `TST-DIAG-003` on 2026-05-09. |
| `tests/unit/test_signalome_workflow_diagnostics.py` | `test_boundary_error_reports_context_table_site_membership_failure_seam` | internal/unstable | `rewrite-to-public-behaviour` | None (current assertions are seam-coupled) | `rewrite` | `TST-DIAG-003` | Rewrite around public workflow error behavior and externally visible diagnostics, not private context-table seam structure. | Completed in `TST-DIAG-003` on 2026-05-09. |
| `tests/unit/test_signalome_workflow_diagnostics.py` | `test_boundary_error_reports_context_table_protein_context_failure_seam` | internal/unstable | `rewrite-to-public-behaviour` | None (current assertions are seam-coupled) | `rewrite` | `TST-DIAG-003` | Replace seam-level payload shape checks with public failure contract assertions at workflow boundary. | Completed in `TST-DIAG-003` on 2026-05-09. |
| `tests/unit/test_kinase_workflow_diagnostics.py` | `test_boundary_error_reports_activity_overlap_edge_case` | internal/unstable | `keep-internal-with-rationale` | Internal scientific failure-mode guard: overlap attrition counts are needed to detect silent collapse paths not otherwise observable via current public outputs | `keep` | `TST-DIAG-002` | Keep as targeted internal guard until equivalent public contract surface exists; do not broaden internal payload locking. | Out of scope in `TST-DIAG-003` (keep decision retained). |
| `tests/unit/test_kinase_workflow_diagnostics.py` | `test_boundary_error_reports_prediction_ensemble_collapse_counts` | internal/unstable | `defer` | None (public ownership of ensemble-collapse diagnostics is not yet explicitly defined) | `defer` | `TST-FOLLOWUP-UNASSIGNED` | Defer rewrite/keep decision until diagnostic ownership policy explicitly states whether ensemble-collapse counts are public workflow contract or internal instrumentation. | Out of scope in `TST-DIAG-003` (deferred by review decision). |

## Regeneration Alignment Note (TST-AUDIT-007)

After regenerating `docs/testing/diagnostic_payload_test_report.md` on 2026-05-09, the remaining
candidate clusters are kinase-diagnostic cases already covered by the decisions in this document:

- `keep-scientific-reproducibility`
- `keep-public-contract`
- `keep-internal-with-rationale`
- `defer`

No additional signalome rewrite/narrow candidates remained after `TST-DIAG-003`.
