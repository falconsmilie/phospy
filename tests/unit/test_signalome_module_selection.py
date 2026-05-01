from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.signalomes.clustering as clustering_module
from phospy.errors import SignalomeModuleCountValidationError
from phospy.scientific_policies import (
    ScientificPolicyId,
    build_signalome_module_candidate_score_policy,
)
from phospy.signalomes.clustering import (
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS,
    cluster_sites_with_diagnostics,
    fit_cluster_labels,
    select_module_count_with_diagnostics,
)
from phospy.signalomes.clustering import exact_python as exact_clustering
from phospy.signalomes.clustering.tree_building import (
    prepare_scoring_values_for_clustering,
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
    profile_degeneracy = clustering_module.summarize_profile_degeneracy(values)
    clustering_values = clustering_module._prepare_scoring_values_for_clustering(values)

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

    empty_result = clustering_module._compute_candidate_cluster_scores(
        candidate_range=range(2, 2),
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        **common_kwargs,
    )
    full_result = clustering_module._compute_candidate_cluster_scores(
        candidate_range=range(2, 3),
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        **common_kwargs,
    )
    sampled_result = clustering_module._compute_candidate_cluster_scores(
        candidate_range=range(2, 3),
        candidate_scoring_policy=clustering_module.SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        **common_kwargs,
    )

    assert isinstance(empty_result, clustering_module._CandidateClusterScoreResult)
    assert isinstance(full_result, clustering_module._CandidateClusterScoreResult)
    assert isinstance(sampled_result, clustering_module._CandidateClusterScoreResult)

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


def test_prepare_scoring_values_imputes_non_finite_and_fully_missing_columns() -> None:
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
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
            ],
            dtype=float,
        ),
    )
    assert diagnostics.policy == (
        SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS
    )
    assert diagnostics.non_finite_input_value_count == 6
    assert diagnostics.missing_after_non_finite_normalization_count == 6
    assert diagnostics.imputed_value_count == 6
    assert diagnostics.fully_missing_column_count == 2


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
    assert np.array_equal(
        np.isfinite(prepared_values), np.ones_like(prepared_values, dtype=bool)
    )
    np.testing.assert_allclose(
        prepared_values,
        np.asarray(
            [
                [1.0, 3.0, 0.0],
                [2.0, 2.0, 0.0],
                [3.0, 4.0, 0.0],
            ],
            dtype=float,
        ),
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
