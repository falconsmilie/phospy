# Legacy Science Gap Audit: Rewrite vs Legacy Archive

- Date: 2026-04-19
- Scope reviewed:
  - Legacy science donors under `legacy_archive/phospy_legacy/{prediction,signalomes,activities}`
  - Rewrite implementation under `src/phospy/`
  - Parity/unit/integration tests under `tests/`
- Contract baselines:
  - ADR-012 (fresh-start rewrite): architecture is rewrite-native
  - ADR-013 (scientific parity policy): parity is scientific, not structural
  - ADR-001/002/003/005/011: active public/workflow/dataset/result contracts

This document is governance-facing. It classifies legacy science components by
their *current rewrite truth* and separates scientific reference value from
active downstream contract authority.

## Status Vocabulary (Normative)

- **Supported now**: implemented and part of the supported public lane.
- **Implemented but not authoritative downstream**: implemented output exists,
  but it is not the decision-driving downstream lane by itself.
- **Partially ported**: meaningful subset is implemented; parity-scoped behavior
  remains missing.
- **Deferred**: intentionally not in the supported lane yet.
- **Historical / superseded**: reference context only; not active guidance.

## Governance Snapshot (2026-04-19)

- The supported kinase downstream score lane is `combined_scores` first, with
  `profile_scores` fallback only when combined scores are unavailable.
- Prediction and signalome consume the same authoritative downstream score lane.
- Legacy scientific modules remain science donors only; legacy architecture is
  not a target structure.
- No ADR-vs-code conflict was found in this audit pass. If a future conflict is
  found, it must be called out explicitly instead of harmonized implicitly.

## Science Matrix (Evidence-Backed)

| Scientific area | Current rewrite state | Status | Evidence (code/tests) | Follow-on |
| --- | --- | --- | --- | --- |
| Profile scoring kernel + kinase profile build | Implemented and active in kinase scoring route. | Supported now | `src/phospy/workflows/kinase/science.py`, `tests/parity/test_kinase_workflow_parity.py` | Baseline guardrails only |
| Optional profile missing-value policy variants from legacy (`median_skipna` lane) | Rewrite supports strict current default (`skipna=False`) only. | Partially ported | `build_kinase_profiles` in `src/phospy/workflows/kinase/science.py` | `SCI-GAP-01` |
| Motif scoring output (`motif_scores`) | Implemented and published in scoring result. | Implemented but not authoritative downstream | `src/phospy/prediction/motif_scoring.py`, `tests/parity/test_prediction_science_parity.py` | Closed for supported lane |
| Profile+motif weighted combination (`combined_scores`, `weights`) | Implemented and selected as authoritative downstream lane when present. | Supported now | `src/phospy/prediction/scoring.py`, `src/phospy/workflows/kinase/executor.py`, `tests/parity/test_kinase_workflow_parity.py` | Closed for supported lane |
| Candidate selection/ranking/prediction assembly | Implemented against resolved downstream matrix. | Supported now | `src/phospy/prediction/candidates.py`, `src/phospy/workflows/kinase/science.py` | Closed for supported lane |
| Adaptive ensemble sampling prediction core | Not implemented; rewrite uses deterministic score ranking lane. | Deferred | No rewrite equivalent of legacy `prediction/sampling_core.py` + `prediction/execution.py` | `SCI-GAP-05` |
| Activity weighted/KSEA kernels | Implemented and parity-backed in supported kinase workflow. | Supported now | `src/phospy/activities/scoring.py`, `tests/parity/test_activity_stage_parity.py`, `docs/architecture/activity_science_port_review.md` | Keep parity fixture lock (`SCI-GAP-11`) |
| Signalome consumption of upstream downstream score lane | Implemented; uses same resolved lane as prediction. | Supported now | `src/phospy/workflows/signalome/interpreter.py`, `tests/unit/test_signalome_workflow_diagnostics.py` | Closed for supported lane |
| Site/protein assignment ambiguity metadata (`top_kinase_candidates`, tie counts) | Implemented with deterministic lexicographic tie policy and ambiguity diagnostics. | Supported now | `src/phospy/workflows/signalome/science.py`, `tests/unit/test_signalome_science_ambiguity.py`, parity contract fixture checks | `SCI-GAP-07` closed for current lane |
| `top_kinase_weights` propagation into weighted-top assignment policy | Weights metadata exists, but weighted-top assignment policy is not implemented. | Implemented but not authoritative downstream | `top_kinase_weights` in `src/phospy/workflows/signalome/science.py`; no weighted-top policy path in executor | `SCI-GAP-08` |
| Signalome module-count auto-selection/clustering lane | Not implemented; rewrite uses dominant-kinase grouping route. | Deferred | No rewrite equivalent of `legacy_archive/phospy_legacy/signalomes/clustering.py` | `SCI-GAP-06` |
| Signalome network policy variants (`positive_only`, etc.) | Signed absolute-threshold behavior is implemented; broader policy surface is absent. | Partially ported | `build_kinase_network` in `src/phospy/workflows/signalome/science.py` | `SCI-GAP-09` |
| Expanded signalome output | Output contract exists but workflow sets `expanded_signalome=None`. | Deferred | `SignalomeWorkflowResult.expanded_signalome` in `src/phospy/api/results.py`; executor returns `None` | `SCI-GAP-10` |

## Authoritative Downstream Lanes

- **Kinase prediction lane**: authoritative matrix is resolved by
  `select_downstream_score_matrix(...)` to `combined_scores` first, then
  `profile_scores` fallback.
- **Signalome score-driven lane**: interpreter resolves and passes the same
  downstream matrix source used by kinase prediction.
- **Not authoritative by itself**:
  - `motif_scores` as a standalone lane
  - `top_kinase_weights` metadata without weighted-top policy execution

## Legacy Scientific Reference vs Active Rewrite Contract

**Legacy scientific reference (allowed):**

- scoring kernels and motifs from legacy prediction modules
- adaptive sampling concepts
- clustering/module-selection heuristics
- weighted-top assignment semantics
- activity/KSEA kernels

**Active rewrite contract (binding):**

- Public API/request/result boundaries in ADR-001/005/011 and implemented
  `src/phospy/api/*`
- Validator/interpreter/executor workflow staging in ADR-002 and
  `src/phospy/workflows/*`
- Fresh-start architecture boundary in ADR-012

Legacy structure under `legacy_archive/phospy_legacy/` is **historical /
superseded architecture context**, not an active contract source.

## Follow-on Queue (Still Valid)

- `SCI-GAP-01`: optional profile missing-value strategy lane.
- `SCI-GAP-05`: adaptive ensemble sampling lane.
- `SCI-GAP-06`: signalome clustering/module-count selection diagnostics.
- `SCI-GAP-08`: weighted-top assignment policy.
- `SCI-GAP-09`: explicit network policy variants.
- `SCI-GAP-10`: `expanded_signalome` population and contract fixtures.
- `SCI-GAP-11`: keep activity parity fixtures as blocking regressions.

## Historical / Superseded Notes

- Earlier wording that implied “motif/combined science is missing” is
  superseded.
- Earlier wording that implied “profile-only downstream is intended” is
  superseded.
- Earlier wording that treated legacy package structure as migration target is
  superseded by ADR-012 and current `src/phospy/` implementation boundaries.
