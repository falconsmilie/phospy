from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.science.signalomes.clustering as clustering_module
from phospy.errors import SignalomeModuleCountValidationError
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.science.signalomes.clustering import (
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE,
    cluster_sites_with_diagnostics,
    fit_cluster_labels,
    select_module_count_with_diagnostics,
)
from phospy.science.signalomes.clustering import exact_python as exact_clustering
from phospy.science.signalomes.clustering.candidate_scoring import (
    _CandidateClusterScoreResult,
    compute_candidate_cluster_scores,
    summarize_profile_degeneracy,
)
from phospy.science.signalomes.clustering.scientific_policies import (
    build_signalome_module_candidate_score_policy,
)
from phospy.science.signalomes.clustering.tree_building import (
    prepare_scoring_values_for_clustering,
    prepare_signalome_clustering_matrix,
    summarize_clustering_missing_value_diagnostics,
)


def test_module_count_selection_accepts_explicit_request_equal_to_site_count() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [1.1, 2.1, 3.1],
            [0.9, 1.9, 2.9],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=3,
    )

    assert diagnostics.strategy == "explicit_module_count"
    assert diagnostics.selected_module_count == 3
    assert diagnostics.requested_module_count == 3
    assert diagnostics.reason == "module_count was provided explicitly by the caller"
    assert diagnostics.candidate_scores == {}
    assert diagnostics.stability_report.status == "not_computable"
    assert (
        diagnostics.stability_report.seed_policy
        == "not_applicable_no_resampling_performed"
    )
    assert diagnostics.stability_report.not_computable_reason == (
        "module_count was provided explicitly; automatic module-count stability "
        "was not evaluated"
    )


def test_module_count_selection_accepts_explicit_request_less_than_site_count() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [1.1, 2.1, 3.1],
            [0.9, 1.9, 2.9],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=2,
    )

    assert diagnostics.strategy == "explicit_module_count"
    assert diagnostics.selected_module_count == 2
    assert diagnostics.requested_module_count == 2


def test_module_count_selection_accepts_explicit_request_of_one() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [1.1, 2.1, 3.1],
            [0.9, 1.9, 2.9],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=1,
    )

    assert diagnostics.strategy == "explicit_module_count"
    assert diagnostics.selected_module_count == 1
    assert diagnostics.requested_module_count == 1


def test_module_count_selection_rejects_explicit_request_above_site_count() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [1.1, 2.1, 3.1],
            [0.9, 1.9, 2.9],
        ],
        dtype=float,
    )

    with pytest.raises(SignalomeModuleCountValidationError) as exc_info:
        select_module_count_with_diagnostics(
            scoring_values=scoring_values,
            requested_module_count=5,
        )

    message = str(exc_info.value)
    assert "field=signalome workflow request config.clustering.module_count" in message
    assert "requested_module_count=5" in message
    assert "available_clustering_site_count=3" in message
    assert (
        "choose a module count between 1 and the number of available clustering sites"
        in message
    )


@pytest.mark.parametrize("requested_module_count", [0, -1])
def test_module_count_selection_rejects_non_positive_explicit_request(
    requested_module_count: int,
) -> None:
    scoring_values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [1.1, 2.1, 3.1],
            [0.9, 1.9, 2.9],
        ],
        dtype=float,
    )

    with pytest.raises(SignalomeModuleCountValidationError) as exc_info:
        select_module_count_with_diagnostics(
            scoring_values=scoring_values,
            requested_module_count=requested_module_count,
        )

    message = str(exc_info.value)
    assert "field=signalome workflow request config.clustering.module_count" in message
    assert f"requested_module_count={requested_module_count}" in message
    assert "available_clustering_site_count=3" in message
    assert (
        "choose a module count between 1 and the number of available clustering sites"
        in message
    )


def test_invalid_explicit_module_count_fails_before_candidate_scoring() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [1.1, 2.1, 3.1],
            [0.9, 1.9, 2.9],
        ],
        dtype=float,
    )
    with pytest.raises(SignalomeModuleCountValidationError):
        select_module_count_with_diagnostics(
            scoring_values=scoring_values,
            requested_module_count=5,
        )


def test_module_count_automatic_selection_surfaces_stable_diagnostics() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 1.1, 1.2],
            [1.1, 1.0, 1.2],
            [0.9, 1.2, 1.0],
            [-1.0, -1.1, -1.2],
            [-1.1, -1.0, -1.2],
            [-0.9, -1.2, -1.0],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(scoring_values=scoring_values)

    assert diagnostics.used_automatic_selection
    assert diagnostics.strategy == "correlation_thresholds"
    assert diagnostics.selected_module_count == 2
    assert diagnostics.threshold_used in {0.5, 0.1}
    assert diagnostics.max_clusters_evaluated >= 2
    assert diagnostics.reason
    assert 2 in diagnostics.candidate_scores


def test_module_selection_stability_report_marks_stable_synthetic_clusters() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 1.1, 1.2],
            [1.1, 1.0, 1.2],
            [0.9, 1.2, 1.0],
            [-1.0, -1.1, -1.2],
            [-1.1, -1.0, -1.2],
            [-0.9, -1.2, -1.0],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        module_selection_stability_seed=17,
        module_selection_stability_perturbations=4,
    )

    report = diagnostics.stability_report
    assert diagnostics.selected_module_count == 2
    assert report.evaluation_method == "seeded_score_perturbation_and_threshold_grid"
    assert report.evaluation_version == "1"
    assert report.seed_policy == "caller_supplied_fixed_seed"
    assert report.random_seed == 17
    assert report.perturbation_count == 4
    assert report.selected_count_frequency == {2: 4}
    assert report.assignment_similarity_metric == "pairwise_coassignment_agreement"
    assert report.assignment_similarity.minimum == pytest.approx(1.0)
    assert report.threshold_sensitivity.selected_count_frequency == {2: 9}
    assert report.threshold_sensitivity.disagrees_with_selected_count is False
    assert report.status == "stable"
    assert report.not_computable_reason is None
    assert report.limitations


def test_module_selection_stability_report_marks_unstable_boundary_case() -> None:
    scoring_values = np.asarray(
        [
            [-2.104189335288666, 1.4263441455628871, -1.5360373292295497],
            [-2.1075141394612467, 0.6644065769633847, 0.506942562289195],
            [1.194306099359353, 0.7167795557673138, -0.08088287860039405],
            [0.6727266535834475, -1.100749224284884, 0.15098996598149328],
            [-0.46293956688440885, 0.23189491925591632, 0.16299788350861275],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        primary_threshold=0.45,
        fallback_threshold=0.1,
        max_clusters=5,
        module_selection_stability_seed=23,
        module_selection_stability_perturbations=4,
    )

    report = diagnostics.stability_report
    assert diagnostics.selected_module_count == 2
    assert report.status == "unstable"
    assert report.selected_count_frequency == {1: 2, 2: 2}
    assert report.assignment_similarity.minimum == pytest.approx(0.4)
    assert report.threshold_sensitivity.disagrees_with_selected_count is True
    assert report.threshold_sensitivity.selected_count_frequency == {2: 6, 1: 3}


def test_module_selection_stability_report_not_computable_for_insufficient_samples() -> (
    None
):
    scoring_values = np.asarray(
        [
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        module_selection_stability_seed=17,
        module_selection_stability_perturbations=4,
    )

    report = diagnostics.stability_report
    assert report.status == "not_computable"
    assert report.seed_policy == "not_applicable_no_resampling_performed"
    assert report.selected_count_frequency == {}
    assert report.assignment_similarity.evaluated_perturbations == 0
    assert report.assignment_similarity.minimum is None
    assert report.threshold_sensitivity.records == ()
    assert report.not_computable_reason == (
        "fewer than three phosphosite profiles are available"
    )


def test_module_selection_stability_report_is_deterministic_under_fixed_seed() -> None:
    scoring_values = np.asarray(
        [
            [-2.104189335288666, 1.4263441455628871, -1.5360373292295497],
            [-2.1075141394612467, 0.6644065769633847, 0.506942562289195],
            [1.194306099359353, 0.7167795557673138, -0.08088287860039405],
            [0.6727266535834475, -1.100749224284884, 0.15098996598149328],
            [-0.46293956688440885, 0.23189491925591632, 0.16299788350861275],
        ],
        dtype=float,
    )

    first = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        primary_threshold=0.45,
        fallback_threshold=0.1,
        max_clusters=5,
        module_selection_stability_seed=23,
        module_selection_stability_perturbations=4,
    )
    second = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        primary_threshold=0.45,
        fallback_threshold=0.1,
        max_clusters=5,
        module_selection_stability_seed=23,
        module_selection_stability_perturbations=4,
    )

    assert first.stability_report == second.stability_report


def test_module_selection_stability_report_tracks_threshold_sensitivity() -> None:
    scoring_values = np.asarray(
        [
            [-2.104189335288666, 1.4263441455628871, -1.5360373292295497],
            [-2.1075141394612467, 0.6644065769633847, 0.506942562289195],
            [1.194306099359353, 0.7167795557673138, -0.08088287860039405],
            [0.6727266535834475, -1.100749224284884, 0.15098996598149328],
            [-0.46293956688440885, 0.23189491925591632, 0.16299788350861275],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        primary_threshold=0.45,
        fallback_threshold=0.1,
        max_clusters=5,
        module_selection_stability_seed=23,
        module_selection_stability_perturbations=4,
    )

    sensitivity = diagnostics.stability_report.threshold_sensitivity
    assert sensitivity.method == "primary_fallback_threshold_grid"
    assert len(sensitivity.records) == 9
    assert {record.primary_threshold for record in sensitivity.records} == {
        0.4,
        0.45,
        0.5,
    }
    assert {record.fallback_threshold for record in sensitivity.records} == {
        0.05,
        0.1,
        0.15,
    }
    assert sensitivity.disagrees_with_selected_count is True


def test_cluster_sites_matches_two_pass_partition_for_selected_module_count() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0, -1.0, -2.0, -3.0],
            "K2": [2.0, 4.0, 6.0, -2.0, -4.0, -6.0],
            "K3": [3.0, 6.0, 9.0, -3.0, -6.0, -9.0],
        },
        index=[f"P{index};S{index};" for index in range(1, 7)],
    )
    scoring_values = scoring_matrix.to_numpy(dtype=float)
    diagnostics = select_module_count_with_diagnostics(scoring_values=scoring_values)
    module_count = int(diagnostics.selected_module_count)
    if module_count == 1:
        baseline_labels = np.ones(scoring_values.shape[0], dtype=int)
    else:
        baseline_labels = fit_cluster_labels(scoring_values, module_count) + 1

    clustered = cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=None,
    )
    observed_labels = clustered.site_clusters.to_numpy(dtype=int, copy=False)
    observed_partition = observed_labels[:, None] == observed_labels[None, :]
    baseline_partition = baseline_labels[:, None] == baseline_labels[None, :]

    assert clustered.module_selection_diagnostics == diagnostics
    assert np.array_equal(observed_partition, baseline_partition)


def _candidate_scoring_test_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "K1": [1.0, 0.9, -1.0, -0.8],
            "K2": [0.0, 0.1, 1.0, 0.8],
            "K3": [0.5, 0.4, -0.5, -0.4],
        },
        index=["P1;S1;", "P2;S2;", "P3;S3;", "P4;S4;"],
    )


def test_sampled_candidate_scoring_auto_module_count_marks_evaluated() -> None:
    clustered = cluster_sites_with_diagnostics(
        scoring_matrix=_candidate_scoring_test_matrix(),
        requested_module_count=None,
        candidate_scoring_policy=(
            clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
        ),
    )

    assert (
        clustered.candidate_scoring_mode
        == clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert clustered.candidate_scoring_evaluated is True
    assert clustered.candidate_scoring_skip_reason is None
    assert isinstance(clustered.candidate_scoring_sampling, dict)


def test_sampled_candidate_scoring_explicit_module_count_marks_skip_reason() -> None:
    clustered = cluster_sites_with_diagnostics(
        scoring_matrix=_candidate_scoring_test_matrix(),
        requested_module_count=2,
        candidate_scoring_policy=(
            clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
        ),
    )

    assert (
        clustered.candidate_scoring_mode
        == clustering_module.SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    )
    assert clustered.candidate_scoring_evaluated is False
    assert (
        clustered.candidate_scoring_skip_reason
        == clustering_module.SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT
    )
    assert clustered.candidate_scoring_sampling is None


def test_full_candidate_scoring_auto_module_count_marks_evaluated() -> None:
    clustered = cluster_sites_with_diagnostics(
        scoring_matrix=_candidate_scoring_test_matrix(),
        requested_module_count=None,
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    )

    assert (
        clustered.candidate_scoring_mode
        == clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    )
    assert clustered.candidate_scoring_evaluated is True
    assert clustered.candidate_scoring_skip_reason is None
    assert clustered.candidate_scoring_sampling is None


def test_explicit_module_count_skips_candidate_scoring_for_full_backend() -> None:
    clustered = cluster_sites_with_diagnostics(
        scoring_matrix=_candidate_scoring_test_matrix(),
        requested_module_count=2,
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    )

    assert (
        clustered.candidate_scoring_mode
        == clustering_module.SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    )
    assert clustered.candidate_scoring_evaluated is False
    assert (
        clustered.candidate_scoring_skip_reason
        == clustering_module.SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT
    )
    assert clustered.candidate_scoring_sampling is None


def test_candidate_scoring_helper_returns_stable_shape_across_all_branches() -> None:
    values = np.asarray(
        [
            [1.0, 0.0, 0.5],
            [0.9, 0.1, 0.4],
            [0.0, 1.0, 0.5],
            [0.1, 0.9, 0.6],
        ],
        dtype=float,
    )
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)

    common_kwargs = {
        "clustering_values": clustering_values,
        "correlation_values": values,
        "profile_degeneracy": profile_degeneracy,
        "n_sites": int(values.shape[0]),
        "scoring_mode": clustering_module.SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        "tree_engine": clustering_module.SIGNALOME_TREE_ENGINE_EXACT,
        "max_exact_tree_sites": None,
        "max_full_candidate_scoring_sites": 10,
    }

    empty_result = compute_candidate_cluster_scores(
        candidate_range=range(2, 2),
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        **common_kwargs,
    )
    full_result = compute_candidate_cluster_scores(
        candidate_range=range(2, 3),
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        **common_kwargs,
    )
    sampled_result = compute_candidate_cluster_scores(
        candidate_range=range(2, 3),
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        **common_kwargs,
    )

    assert isinstance(empty_result, _CandidateClusterScoreResult)
    assert isinstance(full_result, _CandidateClusterScoreResult)
    assert isinstance(sampled_result, _CandidateClusterScoreResult)

    assert empty_result.candidate_scores == {}
    assert empty_result.candidate_labels == {}
    assert (
        empty_result.candidate_scoring_mode
        == clustering_module.SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    )
    assert empty_result.candidate_scoring_evaluated is False
    assert empty_result.candidate_scoring_skip_reason is None
    assert full_result.candidate_scoring_mode == "full"
    assert full_result.candidate_scoring_evaluated is True
    assert full_result.candidate_scoring_skip_reason is None
    assert sampled_result.candidate_scoring_mode == "sampled"
    assert sampled_result.candidate_scoring_evaluated is True
    assert sampled_result.candidate_scoring_skip_reason is None
    assert 2 in full_result.candidate_scores
    assert 2 in sampled_result.candidate_scores


def test_prepare_scoring_values_imputes_partial_missing_with_column_median() -> None:
    values = np.asarray(
        [
            [1.0, np.nan, 3.0],
            [2.0, 2.0, np.nan],
            [3.0, 4.0, 9.0],
        ],
        dtype=float,
    )

    prepared = prepare_scoring_values_for_clustering(values)

    assert np.array_equal(np.isfinite(prepared), np.ones_like(prepared, dtype=bool))
    np.testing.assert_allclose(
        prepared,
        np.asarray(
            [
                [1.0, 3.0, 3.0],
                [2.0, 2.0, 6.0],
                [3.0, 4.0, 9.0],
            ],
            dtype=float,
        ),
    )


def test_prepare_scoring_values_drops_fully_missing_columns() -> None:
    values = np.asarray(
        [
            [1.0, np.nan, np.inf],
            [2.0, np.nan, -np.inf],
            [3.0, np.nan, np.nan],
        ],
        dtype=float,
    )

    prepared = prepare_scoring_values_for_clustering(values)
    diagnostics = summarize_clustering_missing_value_diagnostics(values)

    assert np.array_equal(np.isfinite(prepared), np.ones_like(prepared, dtype=bool))
    np.testing.assert_allclose(
        prepared,
        np.asarray(
            [
                [1.0],
                [2.0],
                [3.0],
            ],
            dtype=float,
        ),
    )
    assert diagnostics.preparation_policy_id == (
        SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE
    )
    assert diagnostics.retained_dimension_labels == ("dimension_0",)
    assert diagnostics.dropped_fully_missing_dimension_count == 2
    assert diagnostics.dropped_fully_missing_dimension_labels == (
        "dimension_1",
        "dimension_2",
    )
    assert diagnostics.dropped_fully_missing_dimension_preview == (
        "dimension_1",
        "dimension_2",
    )
    assert diagnostics.dropped_fully_missing_value_count == 6
    assert diagnostics.non_finite_input_value_count == 6
    assert diagnostics.missing_after_non_finite_normalization_count == 6
    assert diagnostics.imputed_value_count == 0
    assert diagnostics.imputed_value_counts_by_dimension == {"dimension_0": 0}
    assert diagnostics.prepared_matrix_fingerprint is not None


def test_prepare_signalome_clustering_matrix_preserves_retained_column_order() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "K2": [np.nan, np.nan, np.nan],
            "K1": [1.0, np.nan, 3.0],
            "K3": [7.0, 8.0, 9.0],
        },
        index=["P1;S1;", "P2;S2;", "P3;S3;"],
    )

    prepared = prepare_signalome_clustering_matrix(scoring_matrix)

    assert prepared.retained_column_labels == ("K1", "K3")
    assert prepared.prepared_matrix.columns.tolist() == ["K1", "K3"]
    assert prepared.dropped_fully_missing_column_labels == ("K2",)
    np.testing.assert_allclose(
        prepared.values,
        np.asarray(
            [
                [1.0, 7.0],
                [2.0, 8.0],
                [3.0, 9.0],
            ],
            dtype=float,
        ),
    )


def test_prepare_signalome_clustering_matrix_tracks_dropped_and_imputed_cells() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0],
            "K2": [np.nan, np.nan, np.nan],
            "K3": [4.0, np.nan, 8.0],
        },
        index=["P1;S1;", "P2;S2;", "P3;S3;"],
    )

    prepared = prepare_signalome_clustering_matrix(scoring_matrix)
    diagnostics = prepared.to_diagnostics()

    assert diagnostics.dropped_fully_missing_dimension_labels == ("K2",)
    assert diagnostics.dropped_fully_missing_value_count == 3
    assert diagnostics.imputed_value_count == 1
    assert diagnostics.imputed_value_counts_by_dimension == {"K1": 0, "K3": 1}
    np.testing.assert_allclose(
        prepared.prepared_matrix.loc[:, "K3"].to_numpy(dtype=float),
        np.asarray([4.0, 6.0, 8.0], dtype=float),
    )


def test_prepare_signalome_clustering_matrix_fingerprint_is_value_deterministic() -> (
    None
):
    scoring_matrix = pd.DataFrame(
        {
            "K1": [1.0, np.nan, 3.0],
            "K2": [4.0, 5.0, 6.0],
        },
        index=["P1;S1;", "P2;S2;", "P3;S3;"],
    )
    changed_matrix = scoring_matrix.copy(deep=True)
    changed_matrix.loc["P3;S3;", "K2"] = 7.0

    first = prepare_signalome_clustering_matrix(scoring_matrix)
    second = prepare_signalome_clustering_matrix(scoring_matrix)
    changed = prepare_signalome_clustering_matrix(changed_matrix)

    assert first.prepared_matrix_fingerprint == second.prepared_matrix_fingerprint
    assert first.prepared_matrix_fingerprint != changed.prepared_matrix_fingerprint


def test_cluster_sites_rejects_all_missing_dimensions_before_backend() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "K1": [np.nan, np.nan, np.nan],
            "K2": [np.nan, np.nan, np.nan],
        },
        index=["P1;S1;", "P2;S2;", "P3;S3;"],
    )
    backend_invoked = False

    class _FailIfInvokedTreeOperations:
        def build_cluster_tree(self, scoring_values: np.ndarray) -> object:
            nonlocal backend_invoked
            backend_invoked = True
            return object()

        def build_cluster_labels_from_tree(
            self,
            *,
            cluster_tree: object,
            cluster_counts: object,
        ) -> dict[int, np.ndarray]:
            del cluster_tree, cluster_counts
            return {}

    with pytest.raises(ValueError, match="retained no kinase/dimension columns"):
        exact_clustering.cluster_sites_with_diagnostics(
            scoring_matrix=scoring_matrix,
            requested_module_count=2,
            cluster_tree_operations=_FailIfInvokedTreeOperations(),  # type: ignore[arg-type]
        )

    assert backend_invoked is False


def test_cluster_sites_uses_imputed_values_for_final_tree_input() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0],
            "K2": [np.nan, 2.0, 4.0],
            "K3": [np.nan, np.nan, np.nan],
        },
        index=["P1;S1;", "P2;S2;", "P3;S3;"],
    )
    captured: dict[str, np.ndarray] = {}

    class _CaptureTreeOperations:
        def build_cluster_tree(self, scoring_values: np.ndarray) -> object:
            captured["scoring_values"] = np.asarray(scoring_values, dtype=float).copy()
            return object()

        def build_cluster_labels_from_tree(
            self,
            *,
            cluster_tree: object,
            cluster_counts: object,
        ) -> dict[int, np.ndarray]:
            del cluster_tree
            labels: dict[int, np.ndarray] = {}
            for count in [int(value) for value in cluster_counts]:
                if count == 2:
                    labels[count] = np.asarray([0, 0, 1], dtype=int)
                else:
                    labels[count] = np.zeros(3, dtype=int)
            return labels

    exact_clustering.cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=2,
        cluster_tree_operations=_CaptureTreeOperations(),  # type: ignore[arg-type]
    )

    prepared_values = captured.get("scoring_values")
    assert prepared_values is not None
    assert prepared_values.shape == (3, 2)
    assert np.array_equal(
        np.isfinite(prepared_values), np.ones_like(prepared_values, dtype=bool)
    )
    np.testing.assert_allclose(
        prepared_values,
        np.asarray(
            [
                [1.0, 3.0],
                [2.0, 2.0],
                [3.0, 4.0],
            ],
            dtype=float,
        ),
    )


def test_cluster_sites_is_invariant_to_adding_all_missing_dimension() -> None:
    base_matrix = pd.DataFrame(
        {
            "K1": [1.0, 1.2, 0.9, -1.0, -1.2, -0.9],
            "K2": [2.0, 2.2, 1.8, -2.0, -2.2, -1.8],
        },
        index=[f"P{index};S{index};" for index in range(1, 7)],
    )
    with_all_missing = base_matrix.assign(K3=np.nan)

    baseline = cluster_sites_with_diagnostics(
        scoring_matrix=base_matrix,
        requested_module_count=2,
    )
    observed = cluster_sites_with_diagnostics(
        scoring_matrix=with_all_missing,
        requested_module_count=2,
    )

    pd.testing.assert_series_equal(observed.site_clusters, baseline.site_clusters)
    assert (
        observed.clustering_preparation_diagnostics.dropped_fully_missing_dimension_labels
        == ("K3",)
    )
    assert (
        baseline.clustering_preparation_diagnostics.dropped_fully_missing_dimension_labels
        == ()
    )


def test_module_selection_survives_empty_candidate_range_branch() -> None:
    values = np.asarray(
        [
            [1.0, 0.0, 0.5],
            [0.9, 0.1, 0.4],
            [0.0, 1.0, 0.5],
            [0.1, 0.9, 0.6],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=values,
        max_clusters=1,
    )

    assert diagnostics.selected_module_count == 1
    assert diagnostics.candidate_scores == {}
    assert "fewer than two cluster counts are available" in diagnostics.reason


def test_candidate_scoring_policy_record_matches_clustering_execution_metadata() -> (
    None
):
    clustered = cluster_sites_with_diagnostics(
        scoring_matrix=_candidate_scoring_test_matrix(),
        requested_module_count=None,
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    )
    record = build_signalome_module_candidate_score_policy(
        requested_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        candidate_scoring_mode=str(clustered.candidate_scoring_mode),
        max_exact_tree_sites=clustering_module.MAX_FULL_CORRELATION_SITE_COUNT,
        max_full_candidate_scoring_sites=clustering_module.MAX_FULL_CORRELATION_SITE_COUNT,
        candidate_scoring_evaluated=bool(clustered.candidate_scoring_evaluated),
        candidate_scoring_skip_reason=clustered.candidate_scoring_skip_reason,
    )

    assert record.id == ScientificPolicyId.SIGNALOME_MODULE_CANDIDATE_SCORE
    assert record.parameters["candidate_scoring_mode"] == "full"
    assert record.parameters["candidate_scoring_evaluated"] is True
