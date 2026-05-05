"""Score preparation and candidate-scoring collaborators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phospy.signalomes.clustering.candidate_scoring import (
    _CandidateClusterScoreResult,
    _ProfileDegeneracySummary,
    build_correlation_exclusion_note,
    compute_candidate_cluster_scores,
    resolve_candidate_scoring_policy,
    summarize_profile_degeneracy,
)
from phospy.signalomes.clustering.policies import (
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
)
from phospy.signalomes.clustering.tree_building import (
    ClusterTreeOperations,
    prepare_scoring_values_for_clustering,
)


@dataclass(frozen=True, slots=True)
class ScorePreconditioner:
    """Prepare matrix values and profile diagnostics for clustering/scoring."""

    def for_clustering(self, scoring_values: np.ndarray) -> np.ndarray:
        return prepare_scoring_values_for_clustering(scoring_values)

    def profile_degeneracy(
        self, scoring_values: np.ndarray
    ) -> _ProfileDegeneracySummary:
        return summarize_profile_degeneracy(scoring_values)

    def exclusion_note(self, summary: _ProfileDegeneracySummary) -> str:
        return build_correlation_exclusion_note(summary)


@dataclass(frozen=True, slots=True)
class ModuleScorer:
    """Resolve candidate-scoring policy and score candidate module counts."""

    def resolve_policy(
        self,
        *,
        scoring_mode: SignalomeClusteringScoringMode,
        candidate_scoring_policy: SignalomeCandidateScoringPolicy | None,
        n_sites: int,
        max_full_candidate_scoring_sites: int,
    ) -> SignalomeCandidateScoringPolicy:
        return resolve_candidate_scoring_policy(
            scoring_mode=scoring_mode,
            candidate_scoring_policy=candidate_scoring_policy,
            n_sites=n_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        )

    def score_candidates(
        self,
        *,
        clustering_values: np.ndarray,
        correlation_values: np.ndarray,
        candidate_range: range,
        profile_degeneracy: _ProfileDegeneracySummary,
        n_sites: int,
        scoring_mode: SignalomeClusteringScoringMode,
        tree_engine: SignalomeTreeEngine,
        candidate_scoring_policy: SignalomeCandidateScoringPolicy,
        max_exact_tree_sites: int | None,
        max_full_candidate_scoring_sites: int,
        cluster_tree_operations: ClusterTreeOperations | None,
    ) -> _CandidateClusterScoreResult:
        return compute_candidate_cluster_scores(
            clustering_values=clustering_values,
            correlation_values=correlation_values,
            candidate_range=candidate_range,
            profile_degeneracy=profile_degeneracy,
            n_sites=n_sites,
            scoring_mode=scoring_mode,
            tree_engine=tree_engine,
            candidate_scoring_policy=candidate_scoring_policy,
            max_exact_tree_sites=max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
            cluster_tree_operations=cluster_tree_operations,
        )


__all__ = ["ModuleScorer", "ScorePreconditioner"]
