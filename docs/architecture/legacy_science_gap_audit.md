# Legacy Science Gap Audit: Rewrite vs Legacy Archive

- Date: 2026-04-21
- Purpose: governance truth source for legacy-science coverage status in rewrite-native code.
- This audit is intentionally conservative: it is not a blanket claim of whole-project legacy parity.

## Status Vocabulary (Normative)

- `PORTED`: implemented in supported rewrite lanes and backed by rewrite-owned tests.
- `INTENTIONALLY_RETIRED`: intentionally unsupported legacy area.
- `OPEN_GAP`: not yet fully ported in the supported rewrite lane.
- `CONTRACT_CHANGED`: rewrite intentionally narrows/reshapes behavior versus legacy science.
- The current inventory has no `INTENTIONALLY_RETIRED` rows.

## Audit Boundaries

Audited in this pass:

- kinase downstream scoring/prediction lane (including adaptive policy behavior)
- signalome lane (clustering, assignment, network, expanded output)
- activity lane
- dataset preprocessing lane (`missing_data`, `total_protein_correction`, `site_matrix`, `comparisons`)
- rewrite contract seams that intentionally changed legacy behavior

Out of scope for this pass:

- any legacy science surface not listed in the inventory table below
- runtime/performance parity and implementation details of `legacy_archive/`
- claims of full package equivalence beyond the audited seams in this document

## Current Truth Snapshot (2026-04-21)

Landed in current rewrite code (`PORTED`):

- `profile policy behavior`
- `core kinase scoring/prediction lane`
- `signalome clustering/module selection`
- `weighted-top assignment behavior`
- `network policy variants`
- `expanded signalome outputs`
- `activity parity lock`
- `total/protein correction`
- `site-matrix construction`
- `comparison-building`

Out of scope for this audit pass:

- non-inventoried legacy surfaces and non-science parity dimensions (see boundaries above)

Historical gap labels retained for traceability only:

- `SCI-GAP-01/05/06/08/09/10/11/12` are historical ticket labels, not complete coverage inventory.

## Legacy Science Coverage Inventory (Code + Test Evidence)

| Legacy science area | Status | Science-gap ticket | Rewrite-native code evidence | Active rewrite test evidence | Current rewrite truth |
| --- | --- | --- | --- | --- | --- |
| profile policy behavior | PORTED | `SCI-GAP-01` | `src/phospy/workflows/kinase/science.py` (`build_kinase_profiles`) | `tests/unit/test_legacy_donor_science.py::test_profile_policy_donor_locks_strict_median_behavior_and_contract_surface`; `tests/parity/test_kinase_workflow_parity.py::test_profile_missing_value_policy_changes_downstream_lane_for_mixed_missing_input` | `strict` + `median_skipna` profile behavior is supported in rewrite-native scoring and parity-tested. |
| core kinase scoring/prediction lane | PORTED | `SCI-GAP-12` | `src/phospy/workflows/kinase/executor.py`; `src/phospy/prediction/` | `tests/parity/test_l6_prediction_parity.py::test_l6_full_prediction_and_scoring_parity_against_promoted_reference_tables`; `tests/parity/test_adaptive_replay_parity.py::test_adaptive_replay_trace_parity_matches_promoted_trace_surfaces` | Candidate selection, ranking/top-k, and replay surfaces are parity-gated in active rewrite tests. |
| adaptive sampling / svm_mode | CONTRACT_CHANGED | `SCI-GAP-05` | `src/phospy/prediction/execution.py`; `src/phospy/api/configs.py` (`adaptive_policy`) | `tests/unit/test_legacy_donor_science.py::test_adaptive_sampling_donor_is_archival_and_svm_mode_is_not_rewrite_contract`; `tests/parity/test_adaptive_prediction_parity.py::test_adaptive_ensemble_outputs_match_promoted_fixture_tolerances` | Adaptive science is ported, but rewrite public contract uses `adaptive_policy` rather than legacy `svm_mode` naming. |
| signalome clustering/module selection | PORTED | `SCI-GAP-06` | `src/phospy/workflows/signalome/executor.py`; `src/phospy/signalomes/clustering.py` | `tests/unit/test_legacy_donor_science.py::test_signalome_clustering_donor_locks_rewrite_dominant_module_assignment_behavior`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_module_assignments_match_l6_full_fixture_table` | Clustering and module-count diagnostics are implemented in the supported signalome workflow lane. |
| weighted-top assignment behavior | PORTED | `SCI-GAP-08` | `src/phospy/signalomes/modules.py`; `src/phospy/signalomes/assignments.py` | `tests/unit/test_legacy_donor_science.py::test_weighted_top_assignment_donor_locks_fractional_metadata_and_non_fractional_module_selection`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_module_assignments_match_l6_full_fixture_table` | Weighted-top assignment and fractional support propagation are implemented and parity-backed. |
| network policy variants | PORTED | `SCI-GAP-09` | `src/phospy/signalomes/network.py`; `src/phospy/api/configs.py` | `tests/unit/test_legacy_donor_science.py::test_network_policy_variant_donor_locks_signed_edges_and_narrow_config_surface`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_network_policy_variants_match_fixed_matrix_expectations` | `positive_only`, `absolute_threshold`, and `signed` policies are implemented and tested. |
| expanded signalome outputs | PORTED | `SCI-GAP-10` | `src/phospy/signalomes/expanded.py`; `src/phospy/workflows/signalome/executor.py` | `tests/unit/test_legacy_donor_science.py::test_expanded_signalome_donor_locks_supported_lane_to_materialized_output`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_expanded_signalome_matches_l6_full_fixture_table_with_tolerance` | `expanded_signalome` is materialized in the supported workflow executor path. |
| activity parity lock | PORTED | `SCI-GAP-11` | `src/phospy/activities/scoring.py`; `src/phospy/workflows/kinase/executor.py` | `tests/unit/test_legacy_donor_science.py::test_activity_parity_lock_donor_uses_rewrite_owned_fixture_path`; `tests/parity/test_activity_stage_parity.py::test_weighted_activity_matches_rewrite_reference_fixture` | Activity/KSEA science is rewrite-ported and protected by rewrite parity gates. |
| preprocessing transformation establishment | CONTRACT_CHANGED | - | `src/phospy/datasets/builders/transformation_resolver.py`; `src/phospy/transformations/transformers/identity.py` | `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_establishes_transformation_state_via_supported_path`; `tests/unit/test_dataset_transformation_state_establishment.py::test_identity_transformer_is_strict_passthrough_establisher` | Transformation establishment in the supported builder lane is intentionally narrow (`linear` identity establishment). |
| total/protein correction | PORTED | - | `src/phospy/datasets/preprocessing/stages/total_protein_correction.py`; `src/phospy/datasets/preprocessing/pipeline.py` | `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_total_protein_correction_matches_legacy_donor_fixture`; `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_applies_total_protein_correction_when_requested` | `total_protein_correction.policy="ratio_to_total"` is supported with strict phospho/total matching checks. |
| site-matrix construction | PORTED | - | `src/phospy/datasets/preprocessing/stages/site_matrix.py`; `src/phospy/validation/datasets/preprocessing.py`; `src/phospy/api/configs.py` | `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_site_matrix_build_matches_legacy_donor_fixture`; `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_site_matrix_supports_min_observed_and_duplicate_aggregate_mean`; `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_rejects_site_matrix_duplicate_rows_in_error_mode`; `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_supports_site_matrix_duplicate_aggregation_policy` | Supported rewrite now exposes legacy-equivalent site-matrix policy controls (`missing_data_policy`, `minimum_observed_values`, `duplicate_site_strategy`) and enforces deterministic site identity, duplicate handling, and row-drop diagnostics in preprocessing. |
| comparison-building | PORTED | - | `src/phospy/datasets/preprocessing/stages/comparisons.py`; `src/phospy/datasets/preprocessing/pipeline.py` | `tests/unit/test_dataset_preprocessing_subsystem.py::test_dataset_preprocessor_comparison_building_matches_legacy_pairwise_expectation`; `tests/integration/test_dataset_builder_integration.py::test_dataset_builder_builds_inferred_comparisons_from_sample_metadata` | Sample-metadata-based pairwise comparison construction is supported (explicit or inferred pairs). |
| site-to-protein resolution fallback behavior | CONTRACT_CHANGED | - | `src/phospy/workflows/signalome/interpreter.py` (`_resolve_site_to_protein`) | `tests/unit/test_signalome_workflow_diagnostics.py::test_interpreter_does_not_fallback_to_site_id_prefix_when_protein_id_column_missing`; `tests/integration/test_signalome_workflow_integration.py::test_signalome_workflow_requires_explicit_dataset_site_metadata_protein_id` | Signalome requires explicit `site_metadata.protein_id`; no legacy site-id-prefix fallback. |
| signalome input route contraction | CONTRACT_CHANGED | - | `src/phospy/workflows/signalome/public.py`; `src/phospy/api/requests.py` | `tests/unit/test_public_contract_workflows.py::test_workflow_requests_keep_ingestion_outside_workflows`; `tests/integration/test_signalome_workflow_integration.py::test_signalome_workflow_runs_dataset_to_kinase_to_signalome_path` | Supported signalome entrypoint remains contracted to `SignalomeWorkflowRequest(kinase_result=...)`. |
| dataset-vs-reference sequence authority decisions | CONTRACT_CHANGED | - | `src/phospy/workflows/kinase/interpreter.py` (`_resolve_scoring_site_index`) | `tests/integration/test_kinase_workflow_integration.py::test_kinase_workflow_runs_without_dataset_site_sequence_column`; `tests/integration/test_kinase_workflow_integration.py::test_kinase_workflow_runs_dataset_to_kinase_path` | Motif sequence authority in the supported kinase lane is `references.site_sequences`, not dataset-sequence fallback. |

## Open Legacy-Science Areas

- No `OPEN_GAP` areas remain in this audited inventory as of 2026-04-21.

## Maintenance Rule (Governance)

- A science area may be marked `PORTED` only when this audit lists:
  rewrite-native code paths and active rewrite tests for the supported lane.
- Closed ticket state alone is insufficient to mark an area `PORTED`.
- Areas omitted from this audit must be listed as explicit exclusions here or tracked as open elsewhere.
- `tests/unit/test_legacy_donor_inventory.py` is the guard that keeps area/status inventory rows synchronized across parity governance docs.

## Historical Notes

- Historical ticket labels are retained for traceability and changelog continuity.
- Historical ticket closure is not interpreted as whole-surface scientific parity.
