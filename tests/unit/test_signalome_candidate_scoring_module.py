from __future__ import annotations

import numpy as np
import pytest

from phospy.signalomes.clustering.candidate_scoring import (
    _CandidateClusterScoreResult,
    compute_candidate_cluster_scores,
    resolve_candidate_scoring_policy,
    summarize_profile_degeneracy,
)
from phospy.signalomes.clustering.policies import (
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_TREE_ENGINE_EXACT,
)
from phospy.signalomes.clustering.tree_building import (
    prepare_scoring_values_for_clustering,
)


def _values() -> np.ndarray:
    return np.asarray(
        [
            [1.0, 0.0, 0.5],
            [0.9, 0.1, 0.4],
            [0.0, 1.0, 0.5],
            [0.1, 0.9, 0.6],
        ],
        dtype=float,
    )


def test_resolve_candidate_scoring_policy_auto_switches_by_site_limit() -> None:
    assert (
        resolve_candidate_scoring_policy(
            scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
            candidate_scoring_policy=None,
            n_sites=10,
            max_full_candidate_scoring_sites=10,
        )
        == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    )
    assert (
        resolve_candidate_scoring_policy(
            scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
            candidate_scoring_policy=None,
            n_sites=11,
            max_full_candidate_scoring_sites=10,
        )
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )


def test_compute_candidate_cluster_scores_has_independent_full_and_sampled_paths() -> (
    None
):
    values = _values()
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)

    full_result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 3),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        max_exact_tree_sites=10,
        max_full_candidate_scoring_sites=10,
    )
    sampled_result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 3),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=10,
        max_full_candidate_scoring_sites=10,
    )

    assert isinstance(full_result, _CandidateClusterScoreResult)
    assert isinstance(sampled_result, _CandidateClusterScoreResult)
    assert full_result.candidate_scoring_mode == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    assert (
        sampled_result.candidate_scoring_mode
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert full_result.candidate_scoring_evaluated is True
    assert sampled_result.candidate_scoring_evaluated is True
    assert full_result.candidate_scoring_sampling is None
    assert isinstance(sampled_result.candidate_scoring_sampling, dict)


@pytest.mark.parametrize(
    ("values", "expected_pair_count", "minimum_score"),
    [
        (np.asarray([[1.0, 2.0, 3.0]], dtype=float), 0, 0.0),
        (
            np.asarray(
                [
                    [1.0, 2.0, 3.0],
                    [1.0, 2.0, 3.1],
                ],
                dtype=float,
            ),
            1,
            0.99,
        ),
    ],
)
def test_candidate_scoring_handles_one_and_two_site_inputs(
    values: np.ndarray,
    expected_pair_count: int,
    minimum_score: float,
) -> None:
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)
    result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(1, 2),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )

    assert result.candidate_scoring_mode == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    assert result.candidate_scoring_sampling is not None
    assert result.candidate_scoring_sampling["actual_sampled_pair_count"] == (
        expected_pair_count
    )
    assert result.candidate_scores[1].mean_median_correlation >= minimum_score


def test_candidate_scoring_excludes_degenerate_rows_and_scores_near_identical_profiles() -> (
    None
):
    values = np.asarray(
        [
            [1.0, 1.0, 1.0],  # constant
            [np.nan, np.nan, np.nan],  # all-missing
            [1.0, 1.0000001, 1.0000002],  # near-constant
            [0.0, 1.0, 2.0],  # informative
            [0.1, 1.1, 2.1],  # near-identical informative profile
        ],
        dtype=float,
    )
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)
    result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(1, 2),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )

    assert profile_degeneracy.zero_variance_count == 1
    assert profile_degeneracy.near_constant_count == 1
    assert profile_degeneracy.excluded_count == 3
    assert result.candidate_scores[1].min_median_correlation == pytest.approx(1.0)
    assert result.candidate_scores[1].mean_median_correlation == pytest.approx(1.0)


def test_candidate_scoring_sampled_and_full_paths_match_on_small_input() -> None:
    values = _values()
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)

    full_result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 4),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )
    sampled_result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 4),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )

    assert set(full_result.candidate_scores) == set(sampled_result.candidate_scores)
    for cluster_count, full_score in full_result.candidate_scores.items():
        sampled_score = sampled_result.candidate_scores[cluster_count]
        assert sampled_score.min_median_correlation == pytest.approx(
            full_score.min_median_correlation
        )
        assert sampled_score.mean_median_correlation == pytest.approx(
            full_score.mean_median_correlation
        )


def test_candidate_scoring_sampled_path_is_deterministic_for_fixed_input() -> None:
    x = np.linspace(0.0, 1.0, 300, dtype=float)
    values = np.column_stack((x, x**2, np.sin(3.0 * x)))
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)

    first = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 3),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=2000,
        max_full_candidate_scoring_sites=10,
    )
    second = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 3),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=2000,
        max_full_candidate_scoring_sites=10,
    )

    assert first.candidate_scores[2] == second.candidate_scores[2]
    assert first.candidate_scoring_sampling == second.candidate_scoring_sampling
