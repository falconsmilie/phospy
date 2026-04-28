"""Module-count selection facade and helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.signalomes.clustering.orchestration import (
    ClusterTreeOperations,
    SignalomeCandidateScoringBackend,
    SignalomeClusteringScoringMode,
    SignalomeClusterTreeBackend,
    _CandidateScoringMode,
    _ProfileDegeneracySummary,
)
from phospy.signalomes.models import (
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionDiagnostics,
)


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_cluster_tree_sites: int | None = 2000,
) -> int:
    """Select a module count from a scoring matrix."""

    return select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    ).selected_module_count


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_cluster_tree_sites: int | None = 2000,
) -> SignalomeModuleSelectionDiagnostics:
    """Select a module count and return diagnostics."""

    from phospy.signalomes.clustering import orchestration

    return orchestration.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    )


def _compute_module_selection(
    *,
    scoring_values: np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    cluster_tree_backend: SignalomeClusterTreeBackend = "exact",
    candidate_scoring_backend: SignalomeCandidateScoringBackend | None = None,
    max_exact_cluster_tree_sites: int | None = 2000,
    max_full_correlation_sites: int = 2000,
    cluster_tree_operations: ClusterTreeOperations | None = None,
):
    """Compatibility forwarding hook for internal selection computation."""

    from phospy.signalomes.clustering import orchestration

    return orchestration._compute_module_selection(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
        max_full_correlation_sites=max_full_correlation_sites,
        cluster_tree_operations=cluster_tree_operations,
    )


def _resolve_pre_scoring_module_selection(
    *,
    requested_module_count: int | None,
    n_sites: int,
    max_clusters: int,
    profile_degeneracy: _ProfileDegeneracySummary,
    correlation_exclusion_note: str,
):
    """Compatibility forwarding hook for internal pre-scoring resolution."""

    from phospy.signalomes.clustering import orchestration

    return orchestration._resolve_pre_scoring_module_selection(
        requested_module_count=requested_module_count,
        n_sites=n_sites,
        max_clusters=max_clusters,
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
    )


def _compute_candidate_cluster_scores(
    *,
    clustering_values: np.ndarray,
    correlation_values: np.ndarray,
    candidate_range: range,
    profile_degeneracy: _ProfileDegeneracySummary,
    n_sites: int,
    scoring_mode: SignalomeClusteringScoringMode,
    cluster_tree_backend: SignalomeClusterTreeBackend,
    candidate_scoring_backend: SignalomeCandidateScoringBackend,
    max_exact_cluster_tree_sites: int | None,
    max_full_correlation_sites: int,
    cluster_tree_operations: ClusterTreeOperations | None = None,
):
    """Compatibility forwarding hook for candidate score computation."""

    from phospy.signalomes.clustering import orchestration

    return orchestration._compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=correlation_values,
        candidate_range=candidate_range,
        profile_degeneracy=profile_degeneracy,
        n_sites=n_sites,
        scoring_mode=scoring_mode,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
        max_full_correlation_sites=max_full_correlation_sites,
        cluster_tree_operations=cluster_tree_operations,
    )


def _resolve_candidate_scoring_backend(
    *,
    scoring_mode: SignalomeClusteringScoringMode,
    candidate_scoring_backend: SignalomeCandidateScoringBackend | None,
    n_sites: int,
    max_full_correlation_sites: int,
):
    """Compatibility forwarding hook for candidate scoring backend resolution."""

    from phospy.signalomes.clustering import orchestration

    return orchestration._resolve_candidate_scoring_backend(
        scoring_mode=scoring_mode,
        candidate_scoring_backend=candidate_scoring_backend,
        n_sites=n_sites,
        max_full_correlation_sites=max_full_correlation_sites,
    )


def _select_best_candidate_count(candidate_scores: dict[int, float]) -> int:
    from phospy.signalomes.clustering import orchestration

    return orchestration._select_best_candidate_count(candidate_scores)


def _select_threshold_candidate(
    *,
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    candidate_labels: dict[int, np.ndarray],
    max_clusters: int,
    threshold: float,
    requested_module_count: int | None,
    reason: str,
    profile_degeneracy: _ProfileDegeneracySummary,
    correlation_exclusion_note: str,
    approximation_note: str,
    cluster_tree_backend: str,
    candidate_scoring_mode: _CandidateScoringMode,
    exact_cluster_tree_built: bool,
    candidate_scoring_evaluated: bool,
    candidate_scoring_skip_reason: str | None,
    candidate_scoring_sampling: dict[str, object] | None,
):
    from phospy.signalomes.clustering import orchestration

    return orchestration._select_threshold_candidate(
        candidate_scores=candidate_scores,
        candidate_labels=candidate_labels,
        max_clusters=max_clusters,
        threshold=threshold,
        requested_module_count=requested_module_count,
        reason=reason,
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
        approximation_note=approximation_note,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        candidate_scoring_evaluated=candidate_scoring_evaluated,
        candidate_scoring_skip_reason=candidate_scoring_skip_reason,
        candidate_scoring_sampling=candidate_scoring_sampling,
    )


def filter_cluster_candidates(
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    *,
    threshold: float,
) -> dict[int, float]:
    from phospy.signalomes.clustering import orchestration

    return orchestration.filter_cluster_candidates(
        candidate_scores=candidate_scores,
        threshold=threshold,
    )


__all__ = [
    "filter_cluster_candidates",
    "select_module_count",
    "select_module_count_with_diagnostics",
]
