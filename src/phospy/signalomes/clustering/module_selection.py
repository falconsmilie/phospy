"""Module-count selection orchestration collaborator."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phospy.signalomes.clustering.candidate_selection import (
    _ModuleSelectionComputation,
    build_module_selection_result,
    resolve_pre_scoring_module_selection,
    select_threshold_candidate,
)
from phospy.signalomes.clustering.policies import (
    MAX_FULL_CORRELATION_SITE_COUNT,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE,
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_CLUSTERING_SCORING_MODE_EXACT,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
)
from phospy.signalomes.clustering.scale_guards import resolve_max_exact_tree_sites
from phospy.signalomes.clustering.scoring import ModuleScorer, ScorePreconditioner
from phospy.signalomes.clustering.tree_building import ClusterTreeOperations
from phospy.signalomes.clustering.validation import validate_requested_module_count
from phospy.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
)


def _validate_threshold(value: float, *, field_name: str) -> None:
    if not np.isfinite(float(value)):
        raise ValueError(f"{field_name} must be finite")
    if float(value) < 0.0 or float(value) > 1.0:
        raise ValueError(f"{field_name} must be between 0 and 1")


@dataclass(frozen=True, slots=True)
class ModuleSelector:
    """Select the signalome module count and return diagnostics payloads."""

    preconditioner: ScorePreconditioner
    scorer: ModuleScorer

    def select(
        self,
        *,
        scoring_values: np.ndarray,
        requested_module_count: int | None = None,
        primary_threshold: float = 0.5,
        fallback_threshold: float = 0.1,
        max_clusters: int = 10,
        scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
        max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
        max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
        cluster_tree_operations: ClusterTreeOperations | None = None,
    ) -> _ModuleSelectionComputation:
        _validate_threshold(primary_threshold, field_name="primary_threshold")
        _validate_threshold(fallback_threshold, field_name="fallback_threshold")
        if max_clusters < 1:
            raise ValueError("max_clusters must be >= 1")
        if scoring_mode not in {
            SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
            SIGNALOME_CLUSTERING_SCORING_MODE_EXACT,
            SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE,
        }:
            raise ValueError("scoring_mode must be one of: auto, exact, approximate")
        if tree_engine != SIGNALOME_TREE_ENGINE_EXACT:
            raise ValueError("tree_engine must be 'exact'")
        if candidate_scoring_policy not in {
            None,
            SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        }:
            raise ValueError("candidate_scoring_policy must be one of: full, sampled")
        if max_full_candidate_scoring_sites < 1:
            raise ValueError("max_full_candidate_scoring_sites must be >= 1")

        resolved_max_exact_tree_sites = resolve_max_exact_tree_sites(
            max_exact_tree_sites
        )
        scoring_array = np.asarray(scoring_values, dtype=float)
        if scoring_array.ndim != 2:
            raise ValueError("scoring_values must be a 2D array")

        n_sites = int(scoring_array.shape[0])
        requested_module_count = validate_requested_module_count(
            requested_module_count=requested_module_count,
            available_clustering_site_count=n_sites,
            field_name="signalome workflow request config.clustering.module_count",
        )
        profile_degeneracy = self.preconditioner.profile_degeneracy(scoring_array)
        correlation_exclusion_note = self.preconditioner.exclusion_note(
            profile_degeneracy
        )

        early_selection, resolved_max_clusters = resolve_pre_scoring_module_selection(
            requested_module_count=requested_module_count,
            n_sites=n_sites,
            max_clusters=max_clusters,
            profile_degeneracy=profile_degeneracy,
            correlation_exclusion_note=correlation_exclusion_note,
        )
        if early_selection is not None:
            return early_selection

        resolved_candidate_scoring_policy = self.scorer.resolve_policy(
            scoring_mode=scoring_mode,
            candidate_scoring_policy=candidate_scoring_policy,
            n_sites=n_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        )
        clustering_values = self.preconditioner.for_clustering(scoring_array)
        candidate_score_result = self.scorer.score_candidates(
            clustering_values=clustering_values,
            correlation_values=scoring_array,
            candidate_range=range(2, resolved_max_clusters + 1),
            profile_degeneracy=profile_degeneracy,
            n_sites=n_sites,
            scoring_mode=scoring_mode,
            tree_engine=tree_engine,
            candidate_scoring_policy=resolved_candidate_scoring_policy,
            max_exact_tree_sites=resolved_max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
            cluster_tree_operations=cluster_tree_operations,
        )

        primary_selection = select_threshold_candidate(
            candidate_scores=candidate_score_result.candidate_scores,
            candidate_labels=candidate_score_result.candidate_labels,
            max_clusters=resolved_max_clusters,
            threshold=primary_threshold,
            requested_module_count=requested_module_count,
            reason=(
                "selected the highest-scoring candidate that satisfied the primary "
                "within-cluster correlation threshold"
            ),
            profile_degeneracy=profile_degeneracy,
            correlation_exclusion_note=correlation_exclusion_note,
            approximation_note=candidate_score_result.approximation_note,
            candidate_scoring_mode=candidate_score_result.candidate_scoring_mode,
            exact_cluster_tree_built=candidate_score_result.exact_cluster_tree_built,
            candidate_scoring_evaluated=candidate_score_result.candidate_scoring_evaluated,
            candidate_scoring_skip_reason=candidate_score_result.candidate_scoring_skip_reason,
            tree_engine=tree_engine,
            candidate_scoring_sampling=candidate_score_result.candidate_scoring_sampling,
        )
        if primary_selection is not None:
            return primary_selection

        fallback_selection = select_threshold_candidate(
            candidate_scores=candidate_score_result.candidate_scores,
            candidate_labels=candidate_score_result.candidate_labels,
            max_clusters=resolved_max_clusters,
            threshold=fallback_threshold,
            requested_module_count=requested_module_count,
            reason=(
                "no candidate satisfied the primary threshold; selected the "
                "highest-scoring fallback candidate"
            ),
            profile_degeneracy=profile_degeneracy,
            correlation_exclusion_note=correlation_exclusion_note,
            approximation_note=candidate_score_result.approximation_note,
            candidate_scoring_mode=candidate_score_result.candidate_scoring_mode,
            exact_cluster_tree_built=candidate_score_result.exact_cluster_tree_built,
            candidate_scoring_evaluated=candidate_score_result.candidate_scoring_evaluated,
            candidate_scoring_skip_reason=candidate_score_result.candidate_scoring_skip_reason,
            tree_engine=tree_engine,
            candidate_scoring_sampling=candidate_score_result.candidate_scoring_sampling,
        )
        if fallback_selection is not None:
            return fallback_selection

        return build_module_selection_result(
            strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
            selected_module_count=1,
            requested_module_count=requested_module_count,
            threshold_used=None,
            max_clusters_evaluated=resolved_max_clusters,
            candidate_scores=candidate_score_result.candidate_scores,
            reason=(
                "no candidate module count satisfied the configured correlation "
                "thresholds, so the workflow fell back to one module"
            )
            + correlation_exclusion_note
            + candidate_score_result.approximation_note,
            profile_degeneracy=profile_degeneracy,
            excluded_from_correlation_count=profile_degeneracy.excluded_count,
            candidate_labels=candidate_score_result.candidate_labels,
            candidate_scoring_mode=candidate_score_result.candidate_scoring_mode,
            exact_cluster_tree_built=candidate_score_result.exact_cluster_tree_built,
            candidate_scoring_evaluated=candidate_score_result.candidate_scoring_evaluated,
            candidate_scoring_skip_reason=candidate_score_result.candidate_scoring_skip_reason,
            tree_engine=tree_engine,
            candidate_scoring_sampling=candidate_score_result.candidate_scoring_sampling,
        )


__all__ = ["ModuleSelector"]
