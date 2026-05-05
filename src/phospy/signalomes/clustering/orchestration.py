"""Signalome clustering orchestration coordination layer."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.signalomes.clustering.diagnostic_schemas import (
    SignalomeBackendDiagnostics,
    SignalomeCandidateScoringSamplingDiagnostics,
    validate_backend_diagnostics,
    validate_candidate_scoring_sampling_diagnostics,
)
from phospy.signalomes.clustering.diagnostics import (
    approximation_used_from_candidate_mode,
)
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION,
)
from phospy.signalomes.clustering.module_selection import ModuleSelector
from phospy.signalomes.clustering.policies import (
    MAX_FULL_CORRELATION_SITE_COUNT,
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
    _CandidateScoringMode,
)
from phospy.signalomes.clustering.scoring import ModuleScorer, ScorePreconditioner
from phospy.signalomes.clustering.tree_building import (
    ClusterTreeOperations,
    fit_cluster_labels,
)
from phospy.signalomes.clustering.validation import (
    validate_cluster_count_for_site_count,
)
from phospy.signalomes.constants import SITE_CLUSTER_COLUMN
from phospy.signalomes.models import (
    SignalomeModuleSelectionDiagnostics,
)


@dataclass(frozen=True, slots=True)
class ClusterSitesResult:
    """Cluster labels and diagnostics for module-count selection."""

    site_clusters: pd.Series
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics
    tree_engine: str = SIGNALOME_TREE_ENGINE_EXACT
    candidate_scoring_mode: _CandidateScoringMode = (
        SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    )
    exact_cluster_tree_built: bool = False
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None = (
        None
    )
    candidate_scoring_evaluated: bool = False
    candidate_scoring_skip_reason: str | None = None
    backend_name: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    backend_version: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION
    backend_diagnostics: SignalomeBackendDiagnostics | None = None
    approximation_used: bool = False

    def __post_init__(self) -> None:
        if self.candidate_scoring_sampling is not None:
            validate_candidate_scoring_sampling_diagnostics(
                self.candidate_scoring_sampling,
                field_name="cluster_sites_result.candidate_scoring_sampling",
            )
        if self.backend_diagnostics is not None:
            validate_backend_diagnostics(
                self.backend_diagnostics,
                field_name="cluster_sites_result.backend_diagnostics",
            )


@dataclass(frozen=True, slots=True)
class ClusterTreeBuilder:
    """Build final cluster labels for the selected module count."""

    preconditioner: ScorePreconditioner

    def final_labels(
        self,
        *,
        scoring_values: np.ndarray,
        module_count: int,
        cached_labels: np.ndarray | None,
        tree_engine: SignalomeTreeEngine,
        candidate_scoring_policy: SignalomeCandidateScoringPolicy | None,
        max_exact_tree_sites: int | None,
        cluster_tree_operations: ClusterTreeOperations | None,
    ) -> tuple[np.ndarray, bool]:
        n_sites = int(scoring_values.shape[0])
        if module_count == 1:
            return np.ones(n_sites, dtype=int), False
        if cached_labels is not None:
            return cached_labels.astype(int, copy=False) + 1, True

        labels = fit_cluster_labels(
            scoring_values=self.preconditioner.for_clustering(scoring_values),
            cluster_count=module_count,
            tree_engine=tree_engine,
            candidate_scoring_policy=(
                candidate_scoring_policy
                if candidate_scoring_policy is not None
                else SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
            ),
            max_exact_tree_sites=max_exact_tree_sites,
            cluster_tree_operations=cluster_tree_operations,
        )
        return labels + 1, True


_SCORE_PRECONDITIONER = ScorePreconditioner()
_MODULE_SCORER = ModuleScorer()
_MODULE_SELECTOR = ModuleSelector(
    preconditioner=_SCORE_PRECONDITIONER,
    scorer=_MODULE_SCORER,
)
_TREE_BUILDER = ClusterTreeBuilder(preconditioner=_SCORE_PRECONDITIONER)


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> pd.Series:
    """Cluster phosphosites into site clusters."""

    return cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        cluster_tree_operations=cluster_tree_operations,
    ).site_clusters


def cluster_sites_with_diagnostics(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> ClusterSitesResult:
    """Cluster phosphosites and capture module-selection diagnostics."""

    scoring_values = np.asarray(scoring_matrix.to_numpy(dtype=float, copy=False))
    selection = _MODULE_SELECTOR.select(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        cluster_tree_operations=cluster_tree_operations,
    )

    module_count = validate_cluster_count_for_site_count(
        cluster_count=int(selection.diagnostics.selected_module_count),
        available_clustering_site_count=int(scoring_values.shape[0]),
        field_name="selected_module_count",
    )
    labels, built_for_final_assignment = _TREE_BUILDER.final_labels(
        scoring_values=scoring_values,
        module_count=module_count,
        cached_labels=selection.candidate_labels.get(module_count),
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        cluster_tree_operations=cluster_tree_operations,
    )
    exact_cluster_tree_built = bool(
        built_for_final_assignment or selection.exact_cluster_tree_built
    )

    return ClusterSitesResult(
        site_clusters=pd.Series(
            labels,
            index=scoring_matrix.index.copy(),
            dtype=int,
            name=SITE_CLUSTER_COLUMN,
        ),
        module_selection_diagnostics=selection.diagnostics,
        tree_engine=selection.tree_engine,
        candidate_scoring_mode=selection.candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        candidate_scoring_evaluated=selection.candidate_scoring_evaluated,
        candidate_scoring_skip_reason=selection.candidate_scoring_skip_reason,
        candidate_scoring_sampling=selection.candidate_scoring_sampling,
        approximation_used=approximation_used_from_candidate_mode(
            candidate_scoring_mode=str(selection.candidate_scoring_mode),
            candidate_scoring_evaluated=bool(selection.candidate_scoring_evaluated),
        ),
    )


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
) -> int:
    """Select a module count from a scoring matrix."""

    return select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
    ).selected_module_count


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
) -> SignalomeModuleSelectionDiagnostics:
    """Select a module count and return diagnostics."""

    array = (
        scoring_values.to_numpy(dtype=float, copy=False)
        if isinstance(scoring_values, pd.DataFrame)
        else np.asarray(scoring_values, dtype=float)
    )
    return _MODULE_SELECTOR.select(
        scoring_values=array,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
    ).diagnostics


__all__ = [
    "ClusterSitesResult",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "select_module_count",
    "select_module_count_with_diagnostics",
]
