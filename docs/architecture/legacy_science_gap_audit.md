# Legacy Science Gap Audit: Rewrite vs Legacy Archive

- Date: 2026-04-22
- Purpose: governance truth source for legacy-science coverage status in
  rewrite-native code.
- This audit is intentionally conservative and scoped. It is not a blanket
  claim of whole-project legacy parity.

> Audience: maintainers and advanced users tracking scientific evidence posture.
> If you are new to PhosPy usage, start at [Getting started](../getting-started/index.md)
> and use this page later for governance-level detail.

## Coverage Tier Vocabulary (Normative)

This audit uses the same tier model as `docs/parity.md`:

- `PARITY_GATED_ACTIVE_SCIENCE`: active parity-focused regression protection.
  This tier requires explicit active rewrite-owned `tests/parity/...` evidence
  in this audit table for the claimed area.
- `DONOR_BACKED_REWRITE_COVERAGE`: rewrite-owned unit/integration support plus
  donor-backed evidence, but not promoted to parity-gated tier.
- `CONTRACT_CHANGED_SUPPORTED_LANE`: intentionally narrowed or reshaped rewrite
  contract relative to legacy behavior.
- `OPEN_SCIENTIFIC_GAP`: unresolved scientific/parity coverage decision.

## Status Vocabulary (Governance Axis)

Legacy status labels are retained for inventory governance:

- `PORTED`: implemented in supported rewrite lanes and backed by rewrite-owned
  tests.
- `INTENTIONALLY_RETIRED`: intentionally unsupported legacy area.
- `OPEN_GAP`: not yet fully ported in the supported rewrite lane.
- `CONTRACT_CHANGED`: rewrite intentionally narrows/reshapes behavior versus
  legacy science.

Status labels and coverage tiers are intentionally separate so closure claims do
not exceed evidence strength.

## Audit Boundaries

Audited in this pass:

- kinase downstream scoring/prediction lane (including adaptive policy behavior)
- signalome lane (clustering, assignment, network, expanded output)
- activity lane
- dataset preprocessing lane (`missing_data`, `total_protein_correction`,
  `site_matrix`, `comparisons`)
- rewrite contract seams that intentionally changed legacy behavior

Out of scope for this pass:

- any legacy science surface not listed in the inventory table below
- runtime/performance parity and implementation details of `legacy_archive/`
- claims of full package equivalence beyond the audited seams in this document

## Current Snapshot (2026-04-22)

Current audited inventory split by tier:

- `PARITY_GATED_ACTIVE_SCIENCE`:
  `profile policy behavior`, `signalome clustering/module selection`,
  `weighted-top assignment behavior`, `network policy variants`,
  `expanded signalome outputs`, `activity parity lock`,
  `total/protein correction`, `site-matrix construction`,
  `comparison-building`
- `DONOR_BACKED_REWRITE_COVERAGE`:
  `core kinase scoring/prediction lane`
- `CONTRACT_CHANGED_SUPPORTED_LANE`:
  `adaptive sampling / svm_mode`,
  `preprocessing transformation establishment`,
  `site-to-protein resolution fallback behavior`,
  `signalome input route contraction`,
  `dataset-vs-reference sequence authority decisions`
- `OPEN_SCIENTIFIC_GAP`:
  no rows in this audited inventory snapshot

Historical ticket labels are retained for traceability only:
`SCI-GAP-01/05/06/08/09/10/11/12` are historical labels, not a substitute for
evidence-tier classification.

## Evidence Posture (2026-04-22)

- Strong evidence of preserved science: profile policy behavior, signalome
  clustering/assignment/network/expanded outputs, activity parity lock, and
  preprocessing lanes (`total_protein_correction`, `site_matrix`,
  `comparisons`) remain parity-gated with active rewrite-owned tests.
- Unresolved parity-infrastructure issues: kinase prediction ranking confidence
  is currently limited by mixed donor-vs-rewrite and policy-vs-policy
  comparison surfaces plus interim loosened ranking gates; this is tracked
  under Ticket 1 and Ticket 2.
- Intentional contract changes: adaptive naming (`adaptive_policy`),
  candidate-filtering defaults, profile-only fallback behavior, motif-missing
  handling, transformation-establishment narrowing, signalome
  protein-identity/entrypoint narrowing, and reference-sequence authority
  choices are intentional supported-lane differences.

## Kinase Scoring/Prediction Rewrite-vs-Legacy Classification (2026-04-22)

| Difference | Rewrite contract | Legacy baseline | Classification |
| --- | --- | --- | --- |
| Profile-only fallback in score combine | Enabled in supported workflow scoring path | Legacy default disabled | Intentional and supported |
| Missing motif-value handling | Profile score is preserved when motif value is missing for that kinase/site cell | No explicit profile-rescue path in legacy combine | Intentional and supported |
| Workflow candidate filtering defaults | Candidate selection uses `score_threshold=0.0`, `inclusion=1` with caller-owned `top_k` | Legacy defaults used `score_threshold=0.8`, `inclusion=20`, `top=50` | Intentional and supported |
| Request/config knobs | Legacy knobs (`allow_profile_only_fallback`, `score_threshold`, `inclusion`, `min_motif_size`, `svm_mode`, `profile_policy`) are out of public contract | Legacy prediction config exposed those knobs | Intentional and supported |
| Ranking parity measurement surface | Ranking checks still include mixed donor-vs-rewrite and policy-vs-policy comparison surfaces in governance reporting | Closure-grade ranking parity requires like-for-like source and policy surfaces | Provisional evidence only until Ticket 1 and Ticket 2 are completed |
| Parity wording drift in this repo | Historical labels still use `donor-vs-rewrite` wording on rewrite-owned fixture surfaces | N/A | Temporary and should be removed for parity clarity |

Unresolved design decisions in this lane at this audit snapshot:
- Ticket 1: normalize ranking comparison surfaces so governance bars are
  source-consistent and policy-consistent.
- Ticket 2: restore closure-grade ranking gates only after Ticket 1 is complete
  and remove interim loosened-threshold logic from closure decisions.

## Legacy Science Coverage Inventory (Code + Test Evidence)

| Legacy science area | Status | Coverage tier | Contract relation | Science-gap ticket | Rewrite-native code evidence | Active rewrite test evidence | Current rewrite truth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| profile policy behavior | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `SCI-GAP-01` | `src/phospy/workflows/kinase/science.py` (`build_kinase_profiles`) | `tests/unit/test_legacy_donor_science.py::test_profile_policy_donor_locks_strict_median_behavior_and_contract_surface`; `tests/parity/test_kinase_workflow_parity.py::test_profile_missing_value_policy_changes_downstream_lane_for_mixed_missing_input` | `strict` + `median_skipna` profile behavior is supported and parity-tested. |
| core kinase scoring/prediction lane | CONTRACT_CHANGED | DONOR_BACKED_REWRITE_COVERAGE | Contract changed (supported defaults and fallback policy differ from legacy) | `SCI-GAP-12` | `src/phospy/workflows/kinase/executor.py`; `src/phospy/prediction/` | `tests/parity/test_l6_prediction_parity.py::test_l6_prediction_matrix_numeric_parity_matches_promoted_reference_surfaces`; `tests/parity/test_l6_prediction_parity.py::test_l6_prediction_matrix_ranking_parity_matches_promoted_reference_surfaces`; `tests/parity/test_l6_prediction_parity.py::test_l6_candidate_selection_parity_matches_promoted_reference_surfaces`; `tests/parity/test_l6_prediction_parity.py::test_l6_ranked_topk_export_parity_matches_promoted_reference_surfaces`; `tests/parity/test_adaptive_replay_parity.py::test_adaptive_replay_trace_parity_matches_promoted_trace_surfaces` | Numeric scoring/candidate/replay evidence is strong on promoted fixtures, but ranking parity remains provisional while comparison surfaces are inconsistent. Do not treat this lane as parity-closed for ranking until Ticket 1 and Ticket 2 are completed. |
| adaptive sampling / svm_mode | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (`adaptive_policy` replaces legacy `svm_mode` naming) | `SCI-GAP-05` | `src/phospy/prediction/execution.py`; `src/phospy/api/configs.py` (`adaptive_policy`) | `tests/unit/test_legacy_donor_science.py::test_adaptive_sampling_donor_is_archival_and_svm_mode_is_not_rewrite_contract`; `tests/parity/test_adaptive_prediction_parity.py::test_adaptive_prediction_matrix_divergence_cross_policy_stable_vs_r_parity`; `tests/parity/test_adaptive_prediction_parity.py::test_adaptive_ranked_output_divergence_cross_policy_stable_vs_r_parity` | Adaptive science is implemented and tested, but rewrite public contract intentionally differs from legacy naming. |
| signalome clustering/module selection | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `SCI-GAP-06` | `src/phospy/workflows/signalome/executor.py`; `src/phospy/signalomes/clustering.py` | `tests/unit/test_legacy_donor_science.py::test_signalome_clustering_donor_locks_rewrite_dominant_module_assignment_behavior`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_module_assignments_match_l6_full_fixture_table` | Clustering and module-count diagnostics are implemented in the supported signalome workflow lane and parity-tested. |
| weighted-top assignment behavior | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `SCI-GAP-08` | `src/phospy/signalomes/modules.py`; `src/phospy/signalomes/assignments.py` | `tests/unit/test_legacy_donor_science.py::test_weighted_top_assignment_donor_locks_fractional_metadata_and_non_fractional_module_selection`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_module_assignments_match_l6_full_fixture_table` | Weighted-top assignment and fractional support propagation are implemented and parity-backed. |
| network policy variants | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `SCI-GAP-09` | `src/phospy/signalomes/network.py`; `src/phospy/api/configs.py` | `tests/unit/test_legacy_donor_science.py::test_network_policy_variant_donor_locks_signed_edges_and_narrow_config_surface`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_network_policy_variants_match_fixed_matrix_expectations` | `positive_only`, `absolute_threshold`, and `signed` policies are implemented and parity-tested. |
| expanded signalome outputs | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `SCI-GAP-10` | `src/phospy/signalomes/expanded.py`; `src/phospy/workflows/signalome/executor.py` | `tests/unit/test_legacy_donor_science.py::test_expanded_signalome_donor_locks_supported_lane_to_materialized_output`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_expanded_signalome_matches_l6_full_fixture_table_with_tolerance` | `expanded_signalome` is materialized in the supported workflow executor path and parity-tested. |
| activity parity lock | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | `SCI-GAP-11` | `src/phospy/activities/scoring.py`; `src/phospy/workflows/kinase/executor.py` | `tests/unit/test_legacy_donor_science.py::test_activity_parity_lock_donor_uses_rewrite_owned_fixture_path`; `tests/parity/test_activity_stage_parity.py::test_weighted_activity_matches_rewrite_reference_fixture` | Activity/KSEA science is rewrite-ported and protected by rewrite parity gates. |
| preprocessing transformation establishment | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (narrow builder state-establishment policy) | - | `src/phospy/datasets/builders/transformation_resolver.py`; `src/phospy/transformations/transformers/identity.py` | `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_establishes_transformation_state_via_supported_path`; `tests/unit/test_dataset_transformation_state_establishment.py::test_identity_transformer_is_strict_passthrough_establisher` | Transformation establishment in the supported builder lane is intentionally narrow (`linear` identity establishment). |
| total/protein correction | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | - | `src/phospy/datasets/preprocessing/stages/total_protein_correction.py`; `src/phospy/datasets/preprocessing/pipeline.py` | `tests/parity/test_preprocessing_science_parity.py::test_ratio_to_total_total_protein_correction_matches_rewrite_reference_fixture`; `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_total_protein_correction_matches_legacy_donor_fixture`; `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_applies_total_protein_correction_when_requested` | `total_protein_correction.policy="ratio_to_total"` is parity-gated with rewrite-owned fixture expectations and strict phospho/total matching checks. |
| site-matrix construction | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | - | `src/phospy/datasets/preprocessing/stages/site_matrix.py`; `src/phospy/validation/datasets/preprocessing.py`; `src/phospy/api/configs.py` | `tests/parity/test_preprocessing_science_parity.py::test_site_matrix_build_from_metadata_matches_rewrite_reference_fixture`; `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_site_matrix_build_matches_legacy_donor_fixture`; `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_site_matrix_supports_min_observed_and_duplicate_aggregate_mean`; `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_rejects_site_matrix_duplicate_rows_in_error_mode` | Supported site-matrix policy controls are parity-gated on rewrite-owned fixture outputs (matrix shape/values, row-retention diagnostics, and site identity). |
| comparison-building | PORTED | PARITY_GATED_ACTIVE_SCIENCE | Legacy-equivalent in supported lane | - | `src/phospy/datasets/preprocessing/stages/comparisons.py`; `src/phospy/datasets/preprocessing/pipeline.py` | `tests/parity/test_preprocessing_science_parity.py::test_comparison_building_explicit_pair_matches_rewrite_reference_fixture`; `tests/parity/test_preprocessing_science_parity.py::test_comparison_building_inferred_pairs_match_rewrite_reference_fixture`; `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_comparison_building_matches_legacy_pairwise_expectation`; `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_builds_inferred_comparisons_from_sample_metadata` | Sample-metadata comparison construction is parity-gated in supported explicit and inferred pair lanes, including pair identity/order and expected numeric output. |
| site-to-protein resolution fallback behavior | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (no legacy site-id-prefix fallback) | - | `src/phospy/workflows/signalome/interpreter.py` (`_resolve_site_to_protein`) | `tests/unit/test_signalome_workflow_diagnostics.py::test_interpreter_does_not_fallback_to_site_id_prefix_when_protein_id_column_missing`; `tests/integration/test_signalome_workflow_integration.py::test_signalome_workflow_requires_explicit_dataset_site_metadata_protein_id` | Signalome requires explicit `site_metadata.protein_id`; no legacy site-id-prefix fallback. |
| signalome input route contraction | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (entrypoint narrowed) | - | `src/phospy/workflows/signalome/public.py`; `src/phospy/api/requests.py` | `tests/unit/test_public_contract_workflows.py::test_workflow_requests_keep_ingestion_outside_workflows`; `tests/integration/test_signalome_workflow_integration.py::test_signalome_workflow_runs_dataset_to_kinase_to_signalome_path` | Supported signalome entrypoint remains contracted to `SignalomeWorkflowRequest(kinase_result=...)`. |
| dataset-vs-reference sequence authority decisions | CONTRACT_CHANGED | CONTRACT_CHANGED_SUPPORTED_LANE | Contract changed (reference bundle is sequence authority) | - | `src/phospy/workflows/kinase/interpreter.py` (`_resolve_scoring_site_index`) | `tests/integration/test_kinase_workflow_integration.py::test_kinase_workflow_runs_without_dataset_site_sequence_column`; `tests/integration/test_kinase_workflow_integration.py::test_kinase_workflow_runs_dataset_to_kinase_path` | Motif sequence authority in the supported kinase lane is `references.site_sequences`, not dataset-sequence fallback. |

## Open Legacy-Science Areas

- No `OPEN_GAP` rows exist in this audited inventory snapshot as of
  `2026-04-22`.
- No `OPEN_GAP` row here should be read as ranking-parity closure for the core
  kinase lane; Ticket 1 and Ticket 2 remain active governance caveats.
- `OPEN_SCIENTIFIC_GAP` remains a required tier for future unresolved surfaces.
- This statement is scoped to this inventory and its audit boundaries only.

## Maintenance Rule (Governance)

- New science areas must be added to this inventory with both a `Status` value
  and an explicit `Coverage tier` value.
- An area should not be described as parity-gated unless active parity-focused
  tests protect it in supported rewrite lanes.
- For `PARITY_GATED_ACTIVE_SCIENCE`, the `Active rewrite test evidence` column
  must include at least one active `tests/parity/...` reference for that area.
- If parity-focused gate evidence is not yet active, classify the area as
  `DONOR_BACKED_REWRITE_COVERAGE`,
  `CONTRACT_CHANGED_SUPPORTED_LANE`, or `OPEN_SCIENTIFIC_GAP` until promoted.
- Do not promote ranking-focused kinase claims to
  `PARITY_GATED_ACTIVE_SCIENCE` while comparison surfaces remain mixed or
  threshold policy remains provisional.
- `PORTED` alone is not sufficient to imply parity-gated closure; use the
  coverage-tier column.
- Closed ticket state alone is insufficient to mark an area `PORTED`.
- Areas omitted from this audit must be listed as explicit exclusions here or
  tracked as open elsewhere.
- `tests/unit/test_legacy_donor_inventory.py` keeps area/status inventory rows
  synchronized across parity governance docs and guards conservative
  parity-gate evidence rules.

## Historical Notes

- Historical ticket labels are retained for traceability and changelog
  continuity.
- Historical ticket closure is not interpreted as whole-surface scientific
  parity closure.

## Where Next

- Tier definitions and active parity gates: [Parity to PhosR](../parity.md)
- Directional project status: [Roadmap](../roadmap.md)
- Maintainer entry hub: [Contributor and maintainer docs](../contributor/index.md)
