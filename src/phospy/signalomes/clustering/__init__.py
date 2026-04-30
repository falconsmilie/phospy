"""Signalome clustering public facade and backend boundary."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.signalomes.clustering import exact_python as _exact
from phospy.signalomes.clustering import selection as _selection
from phospy.signalomes.clustering.backend_dispatch import (
    available_clustering_engines,
    resolve_clustering_engine,
    run_clustering_engine,
)
from phospy.signalomes.clustering.exact_python import (
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE,
    ClusterSitesResult,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
    derive_protein_modules,
)
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION,
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.signalomes.clustering.protocol import SignalomeClusteringEngine

SIGNALOME_TREE_ENGINE_EXACT = _exact.SIGNALOME_TREE_ENGINE_EXACT
MAX_FULL_CORRELATION_SITE_COUNT = _exact.MAX_FULL_CORRELATION_SITE_COUNT
SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO = (
    _exact.SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO
)
SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS = _exact.SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS
SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES = _exact.SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES
SignalomeClusteringMissingValuePolicy = _exact.SignalomeClusteringMissingValuePolicy
build_cluster_labels_from_tree = _exact.build_cluster_labels_from_tree

# Backward-compatibility re-export: keep internal exact backend constants and
# helpers available through `phospy.signalomes.clustering`.
for _name in dir(_exact):
    if _name.startswith("__"):
        continue
    globals()[_name] = getattr(_exact, _name)


def run_signalome_clustering_engine(
    *,
    scoring_matrix: pd.DataFrame,
    site_to_protein: pd.Series,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    tree_engine: str = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: str | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    clustering_engine: str = SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
) -> SignalomeClusteringEngineResult:
    """Run configured signalome clustering backend."""

    return run_clustering_engine(
        request=SignalomeClusteringEngineRequest(
            scoring_matrix=scoring_matrix,
            site_to_protein=site_to_protein,
            requested_module_count=requested_module_count,
            primary_threshold=primary_threshold,
            fallback_threshold=fallback_threshold,
            max_clusters=max_clusters,
            tree_engine=tree_engine,
            candidate_scoring_policy=candidate_scoring_policy,
            max_exact_tree_sites=max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        ),
        clustering_engine=clustering_engine,
    )


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
) -> pd.Series:
    return _exact.cluster_sites(
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
    )


def cluster_sites_with_diagnostics(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
) -> ClusterSitesResult:
    return _exact.cluster_sites_with_diagnostics(
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
    )


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
) -> int:
    return _selection.select_module_count(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
    )


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
):
    return _selection.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
    )


def fit_cluster_labels(
    scoring_values: np.ndarray,
    cluster_count: int,
    *,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy = SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
) -> np.ndarray:
    return _exact.fit_cluster_labels(
        scoring_values=scoring_values,
        cluster_count=cluster_count,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
    )


__all__ = [
    "ClusterSitesResult",
    "MAX_FULL_CORRELATION_SITE_COUNT",
    "SIGNALOME_CANDIDATE_SCORING_APPLIES_TO",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES",
    "SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON",
    "SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION",
    "SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL",
    "SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION",
    "SIGNALOME_TREE_ENGINE_EXACT",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE",
    "SignalomeClusteringEngine",
    "SignalomeClusteringMissingValuePolicy",
    "SignalomeClusteringEngineRequest",
    "SignalomeClusteringEngineResult",
    "available_clustering_engines",
    "build_cluster_labels_from_tree",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "derive_protein_modules",
    "fit_cluster_labels",
    "resolve_clustering_engine",
    "select_module_count",
    "select_module_count_with_diagnostics",
    "run_clustering_engine",
    "run_signalome_clustering_engine",
]
