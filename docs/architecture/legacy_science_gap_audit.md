# Legacy Science Gap Audit: Rewrite vs Legacy Archive

- Date: 2026-04-20
- Purpose: governance truth source for legacy-science coverage status.
- Scope: whole legacy science surface that is currently known, not only promoted
  donor lanes.

## Status Vocabulary (Normative)

- `PORTED`: legacy scientific area is implemented in supported rewrite lanes and
  guarded by rewrite-owned tests.
- `INTENTIONALLY_RETIRED`: legacy area is intentionally unsupported and should
  not be represented as an active parity target.
- `OPEN_GAP`: legacy scientific area is not yet ported in the supported rewrite
  lane.
- `CONTRACT_CHANGED`: rewrite intentionally narrows or reshapes behavior versus
  legacy science.
- The current inventory has no `INTENTIONALLY_RETIRED` rows.

## Scope Clarification

- Closed historical `SCI-GAP-*` tickets are only one subset of this audit.
- Closed tracked tickets do not imply full legacy-science parity.
- Full parity claims require this full inventory to have no `OPEN_GAP` entries.

## Legacy Science Coverage Inventory

| Legacy science area | Status | Science-gap ticket | Current rewrite truth |
| --- | --- | --- | --- |
| profile policy behavior | PORTED | `SCI-GAP-01` | `strict` + `median_skipna` profile behavior is supported and parity-tested. |
| core kinase scoring/prediction lane | PORTED | `SCI-GAP-12` | Candidate/ranking/replay behavior is parity-gated in rewrite-owned tests. |
| adaptive sampling / svm_mode | CONTRACT_CHANGED | `SCI-GAP-05` | Adaptive science is ported, but public contract uses `adaptive_policy` rather than legacy `svm_mode` naming. |
| signalome clustering/module selection | PORTED | `SCI-GAP-06` | Clustering and module-count diagnostics are implemented and parity-backed. |
| weighted-top assignment behavior | PORTED | `SCI-GAP-08` | Weighted-top assignment and fractional support propagation are implemented. |
| network policy variants | PORTED | `SCI-GAP-09` | `positive_only`, `absolute_threshold`, and `signed` are implemented and tested. |
| expanded signalome outputs | PORTED | `SCI-GAP-10` | `expanded_signalome` is materialized in the supported workflow path. |
| activity parity lock | PORTED | `SCI-GAP-11` | Activity/KSEA science is rewrite-ported and guarded by parity CI gates. |
| preprocessing transformation establishment | CONTRACT_CHANGED | - | Builder preprocessing is intentionally narrow: pass-through linear transformation establishment plus limited missing-data policy. |
| total/protein correction | PORTED | - | `total_protein_correction.policy="ratio_to_total"` is supported in builder preprocessing with strict phospho/total matching checks. |
| site-matrix construction | OPEN_GAP | - | Legacy site-matrix construction policy surface is not fully ported. |
| comparison-building | OPEN_GAP | - | Legacy pairwise comparison-building workflow is not in supported rewrite lane. |
| site-to-protein resolution fallback behavior | CONTRACT_CHANGED | - | Signalome requires explicit `site_metadata.protein_id` and does not apply legacy fallback to site-id prefix. |
| signalome input route contraction | CONTRACT_CHANGED | - | Supported signalome entrypoint is contracted to `SignalomeWorkflowRequest(kinase_result=...)`. |
| dataset-vs-reference sequence authority decisions | CONTRACT_CHANGED | - | Motif sequence authority in supported kinase lane is `references.site_sequences`, not dataset-sequence fallback. |

## Open Legacy-Science Areas

- `site-matrix construction`
- `comparison-building`

These are open legacy-science gaps even though currently tracked `SCI-GAP-*`
tickets are closed.

## Historical Notes

- `SCI-GAP-01/05/06/08/09/10/11/12` remain preserved for traceability.
- Historical labels are not treated as complete inventory coverage.
