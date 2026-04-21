from __future__ import annotations

from dataclasses import dataclass

STATUS_PORTED = "PORTED"
STATUS_INTENTIONALLY_RETIRED = "INTENTIONALLY_RETIRED"
STATUS_OPEN_GAP = "OPEN_GAP"
STATUS_CONTRACT_CHANGED = "CONTRACT_CHANGED"

LEGACY_SCIENCE_STATUS_VALUES = (
    STATUS_PORTED,
    STATUS_INTENTIONALLY_RETIRED,
    STATUS_OPEN_GAP,
    STATUS_CONTRACT_CHANGED,
)

REQUIRED_LEGACY_SCIENCE_AREAS = (
    "profile policy behavior",
    "core kinase scoring/prediction lane",
    "adaptive sampling / svm_mode",
    "signalome clustering/module selection",
    "weighted-top assignment behavior",
    "network policy variants",
    "expanded signalome outputs",
    "activity parity lock",
    "preprocessing transformation establishment",
    "total/protein correction",
    "site-matrix construction",
    "comparison-building",
    "site-to-protein resolution fallback behavior",
    "signalome input route contraction",
    "dataset-vs-reference sequence authority decisions",
)

TRACKED_SCIENCE_GAP_TICKETS = (
    "SCI-GAP-01",
    "SCI-GAP-12",
    "SCI-GAP-05",
    "SCI-GAP-06",
    "SCI-GAP-08",
    "SCI-GAP-09",
    "SCI-GAP-10",
    "SCI-GAP-11",
)

# Governance truth source:
# - tracked SCI-GAP tickets are a historical subset only
# - legacy-science coverage status lives on LEGACY_SCIENCE_AREAS below
OPEN_SCIENCE_GAP_TICKETS: tuple[str, ...] = ()
CLOSED_SCIENCE_GAP_TICKETS: tuple[str, ...] = TRACKED_SCIENCE_GAP_TICKETS


@dataclass(frozen=True, slots=True)
class LegacyScienceAreaInventory:
    area: str
    status: str
    status_summary: str
    science_gap_ticket: str | None
    rewrite_unit_tests: tuple[str, ...]
    rewrite_parity_tests: tuple[str, ...]
    rewrite_integration_tests: tuple[str, ...]
    archival_only_tests: tuple[str, ...]
    promoted_fixture_paths: tuple[str, ...]
    provenance_paths: tuple[str, ...]


LEGACY_SCIENCE_AREAS: tuple[LegacyScienceAreaInventory, ...] = (
    LegacyScienceAreaInventory(
        area="profile policy behavior",
        status=STATUS_PORTED,
        status_summary=(
            "Profile missing-value strategies from legacy donor behavior are "
            "ported in the supported rewrite lane."
        ),
        science_gap_ticket="SCI-GAP-01",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_profile_policy_donor_locks_strict_median_behavior_and_contract_surface",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_kinase_workflow_parity.py::"
            "test_profile_missing_value_policy_changes_downstream_lane_for_mixed_missing_input",
        ),
        rewrite_integration_tests=(),
        archival_only_tests=(
            "tests_legacy/test_profiles.py::"
            "test_build_kinase_substrate_profiles_can_skip_missing_values_when_requested",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/rewrite_parity/r_reference_l6/l6_phospho_matrix.csv",
        ),
        provenance_paths=(
            "tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="core kinase scoring/prediction lane",
        status=STATUS_PORTED,
        status_summary=(
            "Core downstream scoring, candidate selection, prediction-matrix "
            "ranking, top-k export ranking, and replay surfaces are parity-gated "
            "on promoted rewrite fixtures."
        ),
        science_gap_ticket="SCI-GAP-12",
        rewrite_unit_tests=(),
        rewrite_parity_tests=(
            "tests/parity/test_l6_prediction_parity.py::"
            "test_l6_full_prediction_and_scoring_parity_against_promoted_reference_tables",
            "tests/parity/test_adaptive_replay_parity.py::"
            "test_adaptive_replay_trace_parity_matches_promoted_trace_surfaces",
        ),
        rewrite_integration_tests=(),
        archival_only_tests=(
            "tests_legacy/test_parity-with_metrics.py::"
            "test_l6_native_prediction_rankings_agree_with_r_reference",
            "tests_legacy/test_parity-with_metrics.py::"
            "test_l6_replayed_prediction_trace_matches_r_sampling_path",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/rewrite_parity/r_reference_l6_prediction/native_profile_scores.csv",
            "tests/fixtures/rewrite_parity/r_reference_l6_prediction/native_combined_scores.csv",
            "tests/fixtures/rewrite_parity/r_reference_l6_prediction/native_combined_weights.csv",
            "tests/fixtures/rewrite_parity/r_reference_l6_prediction/native_candidate_substrates.csv",
            "tests/fixtures/rewrite_parity/r_reference_l6_prediction/native_prediction_top30.csv",
            "tests/fixtures/rewrite_parity/adaptive_sampling_replay/trace_initial_negatives.csv",
            "tests/fixtures/rewrite_parity/adaptive_sampling_replay/trace_iteration_samples.csv",
            "tests/fixtures/rewrite_parity/adaptive_sampling_replay/trace_final_ensemble_predictions.csv",
            "tests/fixtures/rewrite_parity/adaptive_sampling_replay/trace_final_ensemble_top.csv",
            "src/phospy/data/reference_bundles/rat/l6_native/motif_scores.csv",
            "src/phospy/data/reference_bundles/rat/l6_native/motif_sizes.csv",
        ),
        provenance_paths=(
            "tests/fixtures/rewrite_parity/r_reference_l6_prediction/PROVENANCE.md",
            "tests/fixtures/rewrite_parity/adaptive_sampling_replay/PROVENANCE.md",
            "docs/parity.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="adaptive sampling / svm_mode",
        status=STATUS_CONTRACT_CHANGED,
        status_summary=(
            "Adaptive ensemble science is ported, but legacy svm_mode naming is "
            "not the rewrite public contract."
        ),
        science_gap_ticket="SCI-GAP-05",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_adaptive_sampling_donor_is_archival_and_svm_mode_is_not_rewrite_contract",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_adaptive_prediction_parity.py::"
            "test_adaptive_ensemble_outputs_match_promoted_fixture_tolerances",
        ),
        rewrite_integration_tests=(),
        archival_only_tests=(
            "tests_legacy/test_prediction.py::test_predict_accepts_explicit_r_parity_mode",
            "tests_legacy/test_prediction.py::"
            "test_resolve_prediction_sampling_policy_maps_public_modes",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/rewrite_parity/adaptive_sampling_edge/combined_scores.csv",
            "tests/fixtures/rewrite_parity/adaptive_sampling_edge/trace_candidates.csv",
        ),
        provenance_paths=(
            "tests/fixtures/rewrite_parity/adaptive_sampling_edge/PROVENANCE.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="signalome clustering/module selection",
        status=STATUS_PORTED,
        status_summary=(
            "Signalome clustering and module-count selection diagnostics are "
            "implemented and parity-locked."
        ),
        science_gap_ticket="SCI-GAP-06",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_signalome_clustering_donor_locks_rewrite_dominant_module_assignment_behavior",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_module_assignments_match_l6_full_fixture_table",
        ),
        rewrite_integration_tests=(),
        archival_only_tests=(
            "tests_legacy/test_signalomes.py::"
            "test_select_module_count_builds_one_cluster_tree_for_candidate_scoring",
            "tests_legacy/test_signalomes.py::"
            "test_signalome_workflow_accepts_explicit_module_selection_policy",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/public_workflow_reference/signalome_rewrite_l6_modules.csv",
            "tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json",
        ),
        provenance_paths=("docs/parity.md",),
    ),
    LegacyScienceAreaInventory(
        area="weighted-top assignment behavior",
        status=STATUS_PORTED,
        status_summary=(
            "Weighted-top assignment and fractional support propagation are "
            "implemented in the supported signalome lane."
        ),
        science_gap_ticket="SCI-GAP-08",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_weighted_top_assignment_donor_locks_fractional_metadata_and_non_fractional_module_selection",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_module_assignments_match_l6_full_fixture_table",
        ),
        rewrite_integration_tests=(),
        archival_only_tests=(
            "tests_legacy/test_signalomes.py::"
            "test_weighted_top_assignment_policy_propagates_fractional_module_shares",
            "tests_legacy/test_signalomes.py::"
            "test_build_signalome_support_matrix_supports_weighted_top_policy",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/public_workflow_reference/"
            "signalome_rewrite_l6_module_assignments.csv",
            "tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json",
        ),
        provenance_paths=("docs/parity.md",),
    ),
    LegacyScienceAreaInventory(
        area="network policy variants",
        status=STATUS_PORTED,
        status_summary=(
            "Signed, positive_only, and absolute-threshold network policies are "
            "implemented and parity-tested."
        ),
        science_gap_ticket="SCI-GAP-09",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_network_policy_variant_donor_locks_signed_edges_and_narrow_config_surface",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_network_edges_match_l6_full_fixture_table_with_tolerance",
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_network_policy_variants_match_fixed_matrix_expectations",
        ),
        rewrite_integration_tests=(),
        archival_only_tests=(
            "tests_legacy/test_signalomes.py::"
            "test_build_kinase_network_policies_apply_expected_thresholding",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/public_workflow_reference/"
            "signalome_rewrite_l6_network_edges.csv",
            "tests/fixtures/public_workflow_reference/"
            "signalome_rewrite_l6_network_nodes.csv",
        ),
        provenance_paths=("docs/parity.md",),
    ),
    LegacyScienceAreaInventory(
        area="expanded signalome outputs",
        status=STATUS_PORTED,
        status_summary=(
            "Expanded signalome output population is active in the supported "
            "workflow executor path."
        ),
        science_gap_ticket="SCI-GAP-10",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_expanded_signalome_donor_locks_supported_lane_to_materialized_output",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_expanded_signalome_matches_l6_full_fixture_table_with_tolerance",
        ),
        rewrite_integration_tests=(
            "tests/integration/test_signalome_workflow_integration.py::"
            "test_signalome_workflow_runs_dataset_to_kinase_to_signalome_path",
        ),
        archival_only_tests=(
            "tests_legacy/test_signalomes.py::"
            "test_build_expanded_signalomes_uses_neighbor_map_and_preserves_site_order",
            "tests_legacy/test_signalomes.py::"
            "test_signalome_result_expanded_signalomes_materialize_with_parity",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/public_workflow_reference/"
            "signalome_rewrite_l6_expanded_signalome.csv",
            "tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json",
        ),
        provenance_paths=("docs/parity.md",),
    ),
    LegacyScienceAreaInventory(
        area="activity parity lock",
        status=STATUS_PORTED,
        status_summary=(
            "Activity/KSEA kernels are ported and protected by rewrite-owned "
            "parity regression gates."
        ),
        science_gap_ticket="SCI-GAP-11",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_activity_parity_lock_donor_uses_rewrite_owned_fixture_path",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_activity_stage_parity.py::"
            "test_activity_parity_fixture_set_is_present_readable_and_provenanced",
            "tests/parity/test_activity_stage_parity.py::"
            "test_weighted_activity_matches_rewrite_reference_fixture",
            "tests/parity/test_activity_stage_parity.py::"
            "test_ksea_outputs_match_rewrite_reference_fixture",
            "tests/parity/test_activity_stage_parity.py::"
            "test_target_outputs_match_rewrite_reference_fixture",
        ),
        rewrite_integration_tests=(),
        archival_only_tests=(
            "tests_legacy/test_activities.py::"
            "test_weighted_activity_matches_reference_on_large_sparse_input",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/rewrite_parity/r_reference_l6/kinase_activity_matrix.csv",
            "tests/fixtures/rewrite_parity/r_reference_l6/ksea_scores.csv",
            "tests/fixtures/rewrite_parity/r_reference_l6/ksea_counts.csv",
            "tests/fixtures/rewrite_parity/r_reference_l6/kinase_target_counts.csv",
            "tests/fixtures/rewrite_parity/r_reference_l6/kinase_target_table.csv",
        ),
        provenance_paths=(
            "tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="preprocessing transformation establishment",
        status=STATUS_CONTRACT_CHANGED,
        status_summary=(
            "Rewrite builder preprocessing is intentionally narrow: pass-through "
            "linear transformation establishment plus limited missing-data policy."
        ),
        science_gap_ticket=None,
        rewrite_unit_tests=(
            "tests/unit/test_validator_boundaries.py::"
            "test_dataset_build_request_rejects_unknown_missing_data_policy",
            "tests/unit/test_validator_boundaries.py::"
            "test_dataset_build_request_rejects_impute_policy_without_min_observed_values",
            "tests/unit/test_validator_boundaries.py::"
            "test_dataset_build_request_rejects_min_observed_values_for_forbid_policy",
        ),
        rewrite_parity_tests=(),
        rewrite_integration_tests=(
            "tests/integration/test_dataset_builder_integration.py::"
            "test_dataset_builder_establishes_transformation_state_via_supported_path",
            "tests/integration/test_dataset_builder_integration.py::"
            "test_dataset_builder_supports_row_median_missing_data_preprocessing_policy",
        ),
        archival_only_tests=(
            "tests_legacy/test_preprocessing.py::"
            "test_dataset_preprocessing_run_analysis_ready_uses_example_fixture_data",
        ),
        promoted_fixture_paths=(),
        provenance_paths=(
            "docs/api.md",
            "docs/validation.md",
            "docs/parity.md",
            "docs/architecture/legacy_science_gap_audit.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="total/protein correction",
        status=STATUS_PORTED,
        status_summary=(
            "Legacy-style total/protein correction is implemented in the "
            "supported builder preprocessing lane behind explicit policy."
        ),
        science_gap_ticket=None,
        rewrite_unit_tests=(
            "tests/unit/test_dataset_preprocessing_subsystem.py::"
            "test_dataset_preprocessor_total_protein_correction_matches_legacy_donor_fixture",
            "tests/unit/test_dataset_preprocessing_subsystem.py::"
            "test_dataset_preprocessor_rejects_correction_when_proteins_are_unmatched",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_preprocessing_science_parity.py::"
            "test_ratio_to_total_total_protein_correction_matches_rewrite_reference_fixture",
        ),
        rewrite_integration_tests=(
            "tests/integration/test_dataset_builder_integration.py::"
            "test_dataset_builder_applies_total_protein_correction_when_requested",
        ),
        archival_only_tests=(
            "tests_legacy/test_preprocessing.py::"
            "test_correct_phospho_to_protein_and_pairwise_comparisons",
            "tests_legacy/test_preprocessing.py::"
            "test_protein_correction_service_applies_correction_and_pairwise_augmentation",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/rewrite_parity/protein_correction/"
            "legacy_r_reference_corrected_matrix.csv",
        ),
        provenance_paths=(
            "tests/fixtures/rewrite_parity/protein_correction/PROVENANCE.md",
            "docs/parity.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="site-matrix construction",
        status=STATUS_PORTED,
        status_summary=(
            "Legacy site-matrix construction policy surface is ported in the "
            "supported rewrite builder preprocessing lane, including duplicate "
            "and missing-data policy controls."
        ),
        science_gap_ticket=None,
        rewrite_unit_tests=(
            "tests/unit/test_dataset_preprocessing_subsystem.py::"
            "test_dataset_preprocessor_site_matrix_build_matches_legacy_donor_fixture",
            "tests/unit/test_dataset_preprocessing_subsystem.py::"
            "test_dataset_preprocessor_site_matrix_supports_min_observed_and_duplicate_aggregate_mean",
            "tests/unit/test_dataset_preprocessing_subsystem.py::"
            "test_dataset_preprocessor_rejects_site_matrix_duplicate_rows_in_error_mode",
            "tests/unit/test_validator_boundaries.py::"
            "test_dataset_build_request_allows_site_matrix_policy_overrides",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_preprocessing_science_parity.py::"
            "test_site_matrix_build_from_metadata_matches_rewrite_reference_fixture",
        ),
        rewrite_integration_tests=(
            "tests/integration/test_dataset_builder_integration.py::"
            "test_dataset_builder_supports_site_matrix_build_from_metadata_policy",
            "tests/integration/test_dataset_builder_integration.py::"
            "test_dataset_builder_supports_site_matrix_duplicate_aggregation_policy",
        ),
        archival_only_tests=(
            "tests_legacy/test_preprocessing.py::"
            "test_core_preprocessing_config_normalizes_site_matrix_policy_mapping",
            "tests_legacy/test_matrices.py::"
            "test_build_site_matrix_can_require_minimum_observed_values",
            "tests_legacy/test_matrices.py::"
            "test_build_site_matrix_can_aggregate_duplicate_rows_by_mean",
            "tests_legacy/test_preprocessing.py::"
            "test_analysis_ready_builder_full_inputs_reuses_dataset_preprocessing_seam",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/rewrite_parity/site_matrix/"
            "legacy_r_reference_phospho_corrected.csv",
            "tests/fixtures/rewrite_parity/site_matrix/"
            "legacy_r_reference_expected_matrix.csv",
            "tests/fixtures/rewrite_parity/site_matrix/"
            "legacy_r_reference_expected_phosr_input.csv",
        ),
        provenance_paths=(
            "tests/fixtures/rewrite_parity/site_matrix/PROVENANCE.md",
            "docs/parity.md",
            "docs/architecture/legacy_science_gap_audit.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="comparison-building",
        status=STATUS_PORTED,
        status_summary=(
            "Legacy pairwise comparison-building is supported in the builder "
            "preprocessing lane via sample-metadata grouping policy."
        ),
        science_gap_ticket=None,
        rewrite_unit_tests=(
            "tests/unit/test_dataset_preprocessing_subsystem.py::"
            "test_dataset_preprocessor_builds_inferred_comparisons_from_sample_groups",
            "tests/unit/test_dataset_preprocessing_subsystem.py::"
            "test_dataset_preprocessor_comparison_building_matches_legacy_pairwise_expectation",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_preprocessing_science_parity.py::"
            "test_comparison_building_explicit_pair_matches_rewrite_reference_fixture",
            "tests/parity/test_preprocessing_science_parity.py::"
            "test_comparison_building_inferred_pairs_match_rewrite_reference_fixture",
        ),
        rewrite_integration_tests=(
            "tests/integration/test_dataset_builder_integration.py::"
            "test_dataset_builder_builds_inferred_comparisons_from_sample_metadata",
        ),
        archival_only_tests=(
            "tests_legacy/test_preprocessing.py::"
            "test_add_pairwise_comparisons_uses_schema_group_names",
            "tests_legacy/test_preprocessing.py::"
            "test_add_pairwise_comparisons_rejects_reverse_duplicate_pairs_with_custom_mapping",
        ),
        promoted_fixture_paths=(
            "tests/fixtures/rewrite_parity/comparison_building/legacy_pairwise_expected.csv",
        ),
        provenance_paths=(
            "tests/fixtures/rewrite_parity/comparison_building/PROVENANCE.md",
            "docs/parity.md",
            "docs/architecture/legacy_science_gap_audit.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="site-to-protein resolution fallback behavior",
        status=STATUS_CONTRACT_CHANGED,
        status_summary=(
            "Rewrite signalome requires explicit dataset protein identity and "
            "does not use legacy site-id-prefix fallback behavior."
        ),
        science_gap_ticket=None,
        rewrite_unit_tests=(
            "tests/unit/test_signalome_workflow_diagnostics.py::"
            "test_interpreter_uses_explicit_site_metadata_protein_id_when_present",
            "tests/unit/test_signalome_workflow_diagnostics.py::"
            "test_interpreter_does_not_fallback_to_site_id_prefix_when_protein_id_column_missing",
        ),
        rewrite_parity_tests=(),
        rewrite_integration_tests=(
            "tests/integration/test_signalome_workflow_integration.py::"
            "test_signalome_workflow_requires_explicit_dataset_site_metadata_protein_id",
            "tests/integration/test_signalome_workflow_integration.py::"
            "test_signalome_workflow_uses_explicit_dataset_protein_identity_when_present",
        ),
        archival_only_tests=(
            "tests_legacy/test_datasets_models_site_to_protein_mapping.py::"
            "test_resolve_site_to_protein_mapping_falls_back_to_next_complete_column",
        ),
        promoted_fixture_paths=(),
        provenance_paths=(
            "docs/api.md",
            "docs/parity.md",
            "docs/architecture/legacy_science_gap_audit.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="signalome input route contraction",
        status=STATUS_CONTRACT_CHANGED,
        status_summary=(
            "Supported signalome entrypoint is contracted to "
            "SignalomeWorkflowRequest(kinase_result=...) rather than broader "
            "legacy-style direct inputs."
        ),
        science_gap_ticket=None,
        rewrite_unit_tests=(
            "tests/unit/test_public_contract_workflows.py::"
            "test_workflow_requests_keep_ingestion_outside_workflows",
            "tests/unit/test_public_contract_workflows.py::"
            "test_workflow_run_type_contracts_are_request_to_result",
        ),
        rewrite_parity_tests=(),
        rewrite_integration_tests=(
            "tests/integration/test_signalome_workflow_integration.py::"
            "test_signalome_workflow_runs_dataset_to_kinase_to_signalome_path",
        ),
        archival_only_tests=(),
        promoted_fixture_paths=(),
        provenance_paths=(
            "docs/api.md",
            "docs/parity.md",
            "docs/architecture/legacy_science_gap_audit.md",
        ),
    ),
    LegacyScienceAreaInventory(
        area="dataset-vs-reference sequence authority decisions",
        status=STATUS_CONTRACT_CHANGED,
        status_summary=(
            "In the supported kinase lane, motif sequence authority is the "
            "resolved reference bundle (`references.site_sequences`), not dataset "
            "site metadata fallback."
        ),
        science_gap_ticket=None,
        rewrite_unit_tests=(),
        rewrite_parity_tests=(),
        rewrite_integration_tests=(
            "tests/integration/test_kinase_workflow_integration.py::"
            "test_kinase_workflow_runs_without_dataset_site_sequence_column",
            "tests/integration/test_kinase_workflow_integration.py::"
            "test_kinase_workflow_runs_dataset_to_kinase_path",
        ),
        archival_only_tests=(
            "tests_legacy/test_workflow.py::"
            "test_kinase_workflow_limits_motif_and_prediction_outputs_to_sites_with_sequences",
        ),
        promoted_fixture_paths=(),
        provenance_paths=(
            "docs/api.md",
            "docs/validation.md",
            "docs/parity.md",
            "docs/architecture/legacy_science_gap_audit.md",
        ),
    ),
)

OPEN_LEGACY_SCIENCE_AREAS: tuple[str, ...] = tuple(
    entry.area for entry in LEGACY_SCIENCE_AREAS if entry.status == STATUS_OPEN_GAP
)

PORTED_LEGACY_SCIENCE_AREAS: tuple[str, ...] = tuple(
    entry.area for entry in LEGACY_SCIENCE_AREAS if entry.status == STATUS_PORTED
)

CONTRACT_CHANGED_LEGACY_SCIENCE_AREAS: tuple[str, ...] = tuple(
    entry.area
    for entry in LEGACY_SCIENCE_AREAS
    if entry.status == STATUS_CONTRACT_CHANGED
)

INTENTIONALLY_RETIRED_LEGACY_SCIENCE_AREAS: tuple[str, ...] = tuple(
    entry.area
    for entry in LEGACY_SCIENCE_AREAS
    if entry.status == STATUS_INTENTIONALLY_RETIRED
)
