from __future__ import annotations

import numpy as np

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
