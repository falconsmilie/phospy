"""Module-count selection logic for signalome clustering."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phospy.signalomes.clustering.candidate_scoring import _ProfileDegeneracySummary
from phospy.signalomes.clustering.diagnostic_schemas import (
    SignalomeCandidateScoringSamplingDiagnostics,
)
from phospy.signalomes.clustering.diagnostics import build_module_selection_diagnostics
from phospy.signalomes.clustering.policies import (
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT,
    SIGNALOME_TREE_ENGINE_EXACT,
    _CandidateScoringMode,
)
from phospy.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionDiagnostics,
    SignalomeModuleSelectionStrategy,
)


@dataclass(frozen=True, slots=True)
class _ModuleSelectionComputation:
    diagnostics: SignalomeModuleSelectionDiagnostics
    candidate_labels: dict[int, np.ndarray]
    tree_engine: str
    candidate_scoring_mode: _CandidateScoringMode
    exact_cluster_tree_built: bool
    candidate_scoring_evaluated: bool
    candidate_scoring_skip_reason: str | None
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None


def build_module_selection_result(
    *,
    strategy: SignalomeModuleSelectionStrategy,
    selected_module_count: int,
    requested_module_count: int | None,
    threshold_used: float | None,
    max_clusters_evaluated: int,
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    reason: str,
    profile_degeneracy: _ProfileDegeneracySummary,
    excluded_from_correlation_count: int,
    candidate_labels: dict[int, np.ndarray],
    tree_engine: str = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_mode: _CandidateScoringMode = (
        SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    ),
    exact_cluster_tree_built: bool = False,
    candidate_scoring_evaluated: bool = False,
    candidate_scoring_skip_reason: str | None = None,
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None = (
        None
    ),
) -> _ModuleSelectionComputation:
    return _ModuleSelectionComputation(
        diagnostics=build_module_selection_diagnostics(
            strategy=strategy,
            selected_module_count=selected_module_count,
            requested_module_count=requested_module_count,
            threshold_used=threshold_used,
            max_clusters_evaluated=max_clusters_evaluated,
            candidate_scores=candidate_scores,
            reason=reason,
            zero_variance_profile_count=profile_degeneracy.zero_variance_count,
            near_constant_profile_count=profile_degeneracy.near_constant_count,
            excluded_from_correlation_count=excluded_from_correlation_count,
        ),
        candidate_labels=candidate_labels,
        tree_engine=str(tree_engine),
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=bool(exact_cluster_tree_built),
        candidate_scoring_evaluated=bool(candidate_scoring_evaluated),
        candidate_scoring_skip_reason=(
            None
            if candidate_scoring_skip_reason is None
            else str(candidate_scoring_skip_reason)
        ),
        candidate_scoring_sampling=candidate_scoring_sampling,
    )


def resolve_pre_scoring_module_selection(
    *,
    requested_module_count: int | None,
    n_sites: int,
    max_clusters: int,
    profile_degeneracy: _ProfileDegeneracySummary,
    correlation_exclusion_note: str,
) -> tuple[_ModuleSelectionComputation | None, int]:
    if n_sites <= 1:
        return (
            build_module_selection_result(
                strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
                selected_module_count=1,
                requested_module_count=requested_module_count,
                threshold_used=None,
                max_clusters_evaluated=1,
                candidate_scores={},
                reason="single phosphosite input only supports one signalome module",
                profile_degeneracy=profile_degeneracy,
                excluded_from_correlation_count=0,
                candidate_labels={},
            ),
            1,
        )

    if requested_module_count is not None:
        resolved_count = int(requested_module_count)
        return (
            build_module_selection_result(
                strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
                selected_module_count=resolved_count,
                requested_module_count=int(requested_module_count),
                threshold_used=None,
                max_clusters_evaluated=min(int(max_clusters), n_sites),
                candidate_scores={},
                reason="module_count was provided explicitly by the caller",
                profile_degeneracy=profile_degeneracy,
                excluded_from_correlation_count=0,
                candidate_labels={},
                candidate_scoring_evaluated=False,
                candidate_scoring_skip_reason=(
                    SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT
                ),
            ),
            1,
        )

    resolved_max_clusters = min(int(max_clusters), n_sites)
    if resolved_max_clusters < 2:
        return (
            build_module_selection_result(
                strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
                selected_module_count=1,
                requested_module_count=None,
                threshold_used=None,
                max_clusters_evaluated=resolved_max_clusters,
                candidate_scores={},
                reason="fewer than two cluster counts are available for evaluation",
                profile_degeneracy=profile_degeneracy,
                excluded_from_correlation_count=0,
                candidate_labels={},
            ),
            resolved_max_clusters,
        )

    if n_sites - profile_degeneracy.excluded_count <= 1:
        return (
            build_module_selection_result(
                strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
                selected_module_count=1,
                requested_module_count=None,
                threshold_used=None,
                max_clusters_evaluated=1,
                candidate_scores={},
                reason=(
                    "fewer than two non-degenerate phosphosite profiles remained "
                    "after filtering degenerate rows for correlation scoring"
                )
                + correlation_exclusion_note,
                profile_degeneracy=profile_degeneracy,
                excluded_from_correlation_count=profile_degeneracy.excluded_count,
                candidate_labels={},
            ),
            resolved_max_clusters,
        )

    return None, resolved_max_clusters


def select_best_candidate_count(candidate_scores: dict[int, float]) -> int:
    return max(candidate_scores.items(), key=lambda item: (item[1], -item[0]))[0]


def filter_cluster_candidates(
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    *,
    threshold: float,
) -> dict[int, float]:
    """Return candidate counts whose minimum cluster median passes threshold."""

    return {
        cluster_count: score.mean_median_correlation
        for cluster_count, score in candidate_scores.items()
        if score.min_median_correlation >= threshold
    }


def select_threshold_candidate(
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
    tree_engine: str,
    candidate_scoring_mode: _CandidateScoringMode,
    exact_cluster_tree_built: bool,
    candidate_scoring_evaluated: bool,
    candidate_scoring_skip_reason: str | None,
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None,
) -> _ModuleSelectionComputation | None:
    passing_candidates = filter_cluster_candidates(
        candidate_scores,
        threshold=threshold,
    )
    if not passing_candidates:
        return None
    return build_module_selection_result(
        strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        selected_module_count=select_best_candidate_count(passing_candidates),
        requested_module_count=requested_module_count,
        threshold_used=float(threshold),
        max_clusters_evaluated=int(max_clusters),
        candidate_scores=candidate_scores,
        reason=reason + correlation_exclusion_note + approximation_note,
        profile_degeneracy=profile_degeneracy,
        excluded_from_correlation_count=profile_degeneracy.excluded_count,
        candidate_labels=candidate_labels,
        tree_engine=tree_engine,
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        candidate_scoring_evaluated=candidate_scoring_evaluated,
        candidate_scoring_skip_reason=candidate_scoring_skip_reason,
        candidate_scoring_sampling=candidate_scoring_sampling,
    )


__all__ = [
    "build_module_selection_result",
    "filter_cluster_candidates",
    "resolve_pre_scoring_module_selection",
    "select_best_candidate_count",
    "select_threshold_candidate",
]
