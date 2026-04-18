# Legacy Science Gap Audit: Rewrite vs Legacy Archive

- Date: 2026-04-18
- Priority: P1
- Scope reviewed:
  - `legacy_archive/phospy_legacy/prediction/*`
  - `legacy_archive/phospy_legacy/signalomes/*`
  - `legacy_archive/phospy_legacy/activities/*`
  - `src/phospy/workflows/*`
  - rewrite parity fixtures and tests
- ADR alignment: ADR-012 (fresh-start rewrite), ADR-013 (scientific parity policy)

This note defines which legacy scientific components are still valid donors for the rewrite and which legacy structures must remain archived.

## P1 Decision Record (2026-04-18)

This ticket explicitly resolves the three open scientific-route questions:

- **Does motif scoring belong in the supported route now?**  
  **Yes.** Rewrite-native motif scoring is now part of kinase scoring outputs.
- **Should combined profile/motif scoring return?**  
  **Yes.** Rewrite now publishes `combined_scores` and per-kinase `weights`.
- **Is adaptive or richer candidate selection scientifically required now?**  
  **Richer candidate filtering: Yes (ported).**  
  **Adaptive sampling: Not yet required for the supported route.** It remains
  deferred until dependency/runtime policy and reproducibility contracts are
  formalized for the rewrite package.

## Classification Legend

- `port as-is with adaptation`
- `port conceptually but reimplement cleanly`
- `do not port`

## Gap Matrix

| Legacy science component | Current rewrite status | Classification | Parity risk and fixture impact | Follow-on ticket |
| --- | --- | --- | --- | --- |
| Profile scoring kernel (`prediction/scoring.py::score_phosphosite_profiles`) | Implemented in `src/phospy/workflows/kinase/science.py::score_profile_correlations` | `port as-is with adaptation` (already landed) | Low risk. Keep `tests/fixtures/rewrite_parity/r_reference_l6/native_profile_scores.csv` stable. | `SCI-GAP-00` (baseline guardrails only) |
| Kinase profile construction + substrate support counting (`prediction/profiles.py`) | Implemented with strict missing propagation; no policy surface | `port as-is with adaptation` for current default; extend policy in clean API later | Low/medium. Any future `median_skipna` option will change score baselines and `predMat` expectations. | `SCI-GAP-01` |
| Motif scoring (frequency matrices, sequence window scoring, min-max scaling) (`prediction/motif_scoring.py`) | Not implemented in rewrite workflow path | `port conceptually but reimplement cleanly` | High for scoring outputs when enabled. Promote motif fixture set from `tests_legacy/fixtures/r_reference_l6/native_motif_scores.csv` and `native_motif_sizes.csv` into rewrite parity fixtures. | `SCI-GAP-02` |
| Profile + motif weighted combination (`prediction/scoring.py::combine_profile_and_motif_scores`) | Not implemented | `port conceptually but reimplement cleanly` | High for downstream prediction. Add rewrite fixtures for `native_combined_scores.csv` and `native_combined_weights.csv`. | `SCI-GAP-03` |
| Candidate substrate selection (`prediction/candidates.py`: `top`, `score_threshold`, `inclusion`) | Rewrite prediction uses direct top-k per selected kinase only | `port as-is with adaptation` | Medium/high. Will alter sparsity and rank in `predMat` and substrate list outputs. Impacts `tests/fixtures/rewrite_parity/r_reference_l6/predMat.csv`. | `SCI-GAP-04` |
| Adaptive sampling ensemble prediction core (`prediction/sampling_core.py`, `prediction/execution.py`) | Not implemented (rewrite uses deterministic score ranking) | `port conceptually but reimplement cleanly` | High. Changes full prediction distribution and all downstream activity/signalome fixtures (`predMat`, `kinase_activity_matrix.csv`, `ksea_scores.csv`, signalome contract tables). | `SCI-GAP-05` |
| Signalome module-count selection and clustering (`signalomes/clustering.py`) | Not implemented (rewrite derives modules from dominant kinase grouping) | `port conceptually but reimplement cleanly` | High. Module IDs, assignment distribution, and network composition will shift. Impacts `signalome_rewrite_l6_contract.json`, modules, assignments, and edge fixtures. | `SCI-GAP-06` |
| Site/protein assignment tie metadata (`signalomes/assignments.py` top candidates + ambiguity counts) | Partially implemented (candidates + ambiguity metadata present) | `port as-is with adaptation` (already partial; finish parity semantics) | Medium. Tie-distribution counts in `signalome_rewrite_l6_contract.json` must remain deterministic. | `SCI-GAP-07` |
| Weighted-top tie handling (`top_kinase_weights` + `assignment_policy="weighted_top"`) | Not implemented in rewrite | `port conceptually but reimplement cleanly` | High for ambiguous assignments and module shares. Add explicit tie fixtures from `tests_legacy` weighted-top scenarios. | `SCI-GAP-08` |
| Signalome network threshold policies (`positive_only`, `absolute_threshold`, `signed`) | Rewrite supports one signed absolute-threshold behavior | `port conceptually but reimplement cleanly` for optional policies; keep current default | Medium. Only affected when policy surface is expanded; current fixtures remain stable if default unchanged. | `SCI-GAP-09` |
| Expanded signalome generation (`signalomes/assignments.py::build_expanded_signalomes`) | Explicitly not implemented in rewrite (`expanded_signalome=None`) | `port conceptually but reimplement cleanly` | Medium/high for new output contract; add dedicated workflow fixtures instead of mutating existing selected-point contract files. | `SCI-GAP-10` |
| Activity scoring and KSEA-style calculations (`activities/scoring.py`) | Implemented in rewrite with equivalent kernels and parity tests | `port as-is with adaptation` (already landed) | Low. Keep `kinase_activity_matrix.csv`, `ksea_scores.csv`, `ksea_counts.csv`, and target count fixtures as regression anchors. | `SCI-GAP-11` |

## Science Worth Reusing

- Profile correlation scoring and profile construction logic.
- Motif frequency scoring and profile/motif weighting concepts.
- Candidate filtering semantics (`top`, threshold, inclusion) before prediction.
- Adaptive sampling ensemble concept and deterministic RNG policy seams.
- Signalome clustering/module selection heuristics and diagnostics.
- Weighted-top tie semantics for ambiguous site assignment.
- Activity weighted-average and KSEA kernels (already adopted).

## Legacy Structure Not To Reuse

- Legacy package orchestration/wrapper layers in `legacy_archive/phospy_legacy/api/*` and `internal/*`.
- Compatibility-heavy result facades and mutable alias surfaces in legacy `signalomes/results.py`.
- Trace sink ownership plumbing and replay/export packaging as architecture templates.
- Legacy workflow composition classes as direct module templates.

The archive remains a science donor only; implementation must stay in rewrite validator/interpreter/executor boundaries.

## Follow-on Ticket Queue

- `SCI-GAP-01`: Add optional profile missing-value strategy lane (`propagate_any_missing` default, optional `median_skipna`).
- `SCI-GAP-02`: Introduce motif scoring stage inputs/outputs in rewrite scoring domain.
- `SCI-GAP-03`: Add clean profile+motif weighted combination and surface `combined_scores`/`weights`.
- `SCI-GAP-04`: Replace direct prediction top-k selection with candidate filter (`top`, threshold, inclusion).
- `SCI-GAP-05`: Implement adaptive ensemble prediction lane behind rewrite-native contract.
- `SCI-GAP-06`: Introduce signalome clustering and module-count selection with diagnostics.
- `SCI-GAP-07`: Tighten tie metadata parity semantics for deterministic assignment diagnostics.
- `SCI-GAP-08`: Add weighted-top assignment policy and support matrix propagation.
- `SCI-GAP-09`: Add explicit signalome network policy enum and execution paths.
- `SCI-GAP-10`: Implement `expanded_signalome` output model and fixture-backed regression tests.
- `SCI-GAP-11`: Keep activity/KSEA parity fixtures as blocking regressions while prediction lane evolves.

## Parity and Fixture Strategy

- Preserve current rewrite fixture baselines for already-landed science (`native_profile_scores`, activity outputs) as non-regression anchors.
- Promote legacy scientific fixtures into rewrite-owned fixture paths before enabling new lanes:
  - motif/combined scoring fixtures from `tests_legacy/fixtures/r_reference_l6/`
  - weighted-top and clustering edge cases from `tests_legacy` synthetic scenarios
- Expect staged fixture churn in this order:
  1. scoring fixtures (`motif`, `combined`)
  2. prediction fixtures (`predMat`, candidate and trace-derived subsets if retained)
  3. downstream activity fixtures (`kinase_activity_matrix`, `ksea_scores`, counts)
  4. signalome contract fixtures (`signalome_rewrite_l6_*`)

No parity update should require copying legacy architecture; only scientific outputs and stable contracts are parity targets.
