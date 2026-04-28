"""Signalome clustering public facade and backend boundary."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.signalomes.clustering import exact_python as _exact
from phospy.signalomes.clustering import selection as _selection
from phospy.signalomes.clustering.backend_dispatch import (
    available_clustering_backends,
    resolve_clustering_backend,
    run_clustering_backend,
)
from phospy.signalomes.clustering.exact_python import (
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE,
    ClusterSitesResult,
    SignalomeCandidateScoringBackend,
    SignalomeClusteringScoringMode,
    SignalomeClusterTreeBackend,
    derive_protein_modules,
)
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION,
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL_VERSION,
    SignalomeClusteringBackendRequest,
    SignalomeClusteringBackendResult,
)
from phospy.signalomes.clustering.protocol import SignalomeClusteringBackend

SIGNALOME_CLUSTER_TREE_BACKEND_EXACT = _exact.SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
MAX_FULL_CORRELATION_SITE_COUNT = _exact.MAX_FULL_CORRELATION_SITE_COUNT
build_cluster_labels_from_tree = _exact.build_cluster_labels_from_tree

# Backward-compatibility re-export: keep internal exact backend constants and
# helpers available through `phospy.signalomes.clustering`.
for _name in dir(_exact):
    if _name.startswith("__"):
        continue
    globals()[_name] = getattr(_exact, _name)


def run_signalome_clustering_backend(
    *,
    scoring_matrix: pd.DataFrame,
    site_to_protein: pd.Series,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    cluster_tree_backend: str = SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
    candidate_scoring_backend: str | None = None,
    max_exact_cluster_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_correlation_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    backend_name: str = SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
) -> SignalomeClusteringBackendResult:
    """Run configured signalome clustering backend."""

    return run_clustering_backend(
        request=SignalomeClusteringBackendRequest(
            scoring_matrix=scoring_matrix,
            site_to_protein=site_to_protein,
            requested_module_count=requested_module_count,
            primary_threshold=primary_threshold,
            fallback_threshold=fallback_threshold,
            max_clusters=max_clusters,
            cluster_tree_backend=cluster_tree_backend,
            candidate_scoring_backend=candidate_scoring_backend,
            max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
            max_full_correlation_sites=max_full_correlation_sites,
        ),
        backend_name=backend_name,
    )


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    cluster_tree_backend: SignalomeClusterTreeBackend = SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
    candidate_scoring_backend: SignalomeCandidateScoringBackend | None = None,
    max_exact_cluster_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_correlation_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
) -> pd.Series:
    return _exact.cluster_sites(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
        max_full_correlation_sites=max_full_correlation_sites,
    )


def cluster_sites_with_diagnostics(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    cluster_tree_backend: SignalomeClusterTreeBackend = SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
    candidate_scoring_backend: SignalomeCandidateScoringBackend | None = None,
    max_exact_cluster_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_correlation_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
) -> ClusterSitesResult:
    return _exact.cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
        max_full_correlation_sites=max_full_correlation_sites,
    )


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_cluster_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
) -> int:
    return _selection.select_module_count(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    )


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_cluster_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
):
    return _selection.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    )


def fit_cluster_labels(
    scoring_values: np.ndarray,
    cluster_count: int,
    *,
    cluster_tree_backend: SignalomeClusterTreeBackend = SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
    candidate_scoring_backend: SignalomeCandidateScoringBackend = SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
    max_exact_cluster_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
) -> np.ndarray:
    return _exact.fit_cluster_labels(
        scoring_values=scoring_values,
        cluster_count=cluster_count,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    )


__all__ = [
    "ClusterSitesResult",
    "MAX_FULL_CORRELATION_SITE_COUNT",
    "SIGNALOME_CANDIDATE_SCORING_APPLIES_TO",
    "SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON",
    "SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION",
    "SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL",
    "SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL_VERSION",
    "SIGNALOME_CLUSTER_TREE_BACKEND_EXACT",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE",
    "SignalomeClusteringBackend",
    "SignalomeClusteringBackendRequest",
    "SignalomeClusteringBackendResult",
    "available_clustering_backends",
    "build_cluster_labels_from_tree",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "derive_protein_modules",
    "fit_cluster_labels",
    "resolve_clustering_backend",
    "select_module_count",
    "select_module_count_with_diagnostics",
    "run_clustering_backend",
    "run_signalome_clustering_backend",
]
