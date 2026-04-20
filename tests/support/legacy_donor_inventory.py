from __future__ import annotations

from dataclasses import dataclass

REQUIRED_DONOR_AREAS = (
    "profile policy behavior",
    "adaptive sampling / svm_mode",
    "signalome clustering/module selection",
    "weighted-top assignment behavior",
    "network policy variants",
    "expanded signalome outputs",
    "activity parity lock",
)

OPEN_SCIENCE_GAP_TICKETS = (
    "SCI-GAP-01",
    "SCI-GAP-05",
    "SCI-GAP-06",
    "SCI-GAP-08",
    "SCI-GAP-09",
    "SCI-GAP-10",
    "SCI-GAP-11",
)


@dataclass(frozen=True, slots=True)
class LegacyDonorAreaInventory:
    area: str
    science_gap_ticket: str
    rewrite_unit_tests: tuple[str, ...]
    rewrite_parity_tests: tuple[str, ...]
    rewrite_integration_tests: tuple[str, ...]
    archival_only_tests: tuple[str, ...]
    promoted_fixture_paths: tuple[str, ...]
    provenance_paths: tuple[str, ...]


LEGACY_DONOR_AREAS: tuple[LegacyDonorAreaInventory, ...] = (
    LegacyDonorAreaInventory(
        area="profile policy behavior",
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
    LegacyDonorAreaInventory(
        area="adaptive sampling / svm_mode",
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
    LegacyDonorAreaInventory(
        area="signalome clustering/module selection",
        science_gap_ticket="SCI-GAP-06",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_signalome_clustering_donor_locks_rewrite_dominant_module_assignment_behavior",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_module_assignments_match_selected_l6_regression_points",
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
    LegacyDonorAreaInventory(
        area="weighted-top assignment behavior",
        science_gap_ticket="SCI-GAP-08",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_weighted_top_assignment_donor_locks_fractional_metadata_and_non_fractional_module_selection",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_module_assignments_match_selected_l6_regression_points",
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
            "signalome_rewrite_l6_module_assignments_selected.csv",
            "tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json",
        ),
        provenance_paths=("docs/parity.md",),
    ),
    LegacyDonorAreaInventory(
        area="network policy variants",
        science_gap_ticket="SCI-GAP-09",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_network_policy_variant_donor_locks_signed_edges_and_narrow_config_surface",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_network_edges_match_l6_fixture_pairs_and_sign_counts",
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
            "signalome_rewrite_l6_network_edges_selected.csv",
            "tests/fixtures/public_workflow_reference/"
            "signalome_rewrite_l6_network_nodes.csv",
        ),
        provenance_paths=("docs/parity.md",),
    ),
    LegacyDonorAreaInventory(
        area="expanded signalome outputs",
        science_gap_ticket="SCI-GAP-10",
        rewrite_unit_tests=(
            "tests/unit/test_legacy_donor_science.py::"
            "test_expanded_signalome_donor_locks_supported_lane_to_materialized_output",
        ),
        rewrite_parity_tests=(
            "tests/parity/test_signalome_workflow_parity.py::"
            "test_signalome_expanded_slice_matches_l6_selected_akt1_fixture",
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
            "tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json",
        ),
        provenance_paths=("docs/parity.md",),
    ),
    LegacyDonorAreaInventory(
        area="activity parity lock",
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
)
