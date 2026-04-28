"""SciPy-backed signalome clustering backend implementation."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import cast

import numpy as np
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import pdist

from phospy.signalomes.clustering import exact_python as _exact
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL_VERSION,
    SignalomeClusteringBackendRequest,
    SignalomeClusteringBackendResult,
)

_SCIPY_LINKAGE_METHOD = "ward"
_SCIPY_DISTANCE_METRIC = "euclidean"


@dataclass(frozen=True, slots=True)
class _ScipyWardClusterTree:
    n_sites: int
    linkage_matrix: np.ndarray


def _build_scipy_cluster_tree(scoring_values: np.ndarray) -> _ScipyWardClusterTree:
    values = np.asarray(scoring_values, dtype=float)
    if values.ndim != 2:
        raise ValueError("scoring_values must be a 2D array")
    n_sites = int(values.shape[0])
    if n_sites <= 1:
        return _ScipyWardClusterTree(
            n_sites=n_sites,
            linkage_matrix=np.zeros((0, 4), dtype=float),
        )

    condensed_distances = np.asarray(
        pdist(values, metric=_SCIPY_DISTANCE_METRIC),
        dtype=float,
    )
    linkage_matrix = np.asarray(
        linkage(condensed_distances, method=_SCIPY_LINKAGE_METHOD),
        dtype=float,
    )
    return _ScipyWardClusterTree(n_sites=n_sites, linkage_matrix=linkage_matrix)


def _canonicalize_labels(raw_labels: np.ndarray) -> np.ndarray:
    labels = np.asarray(raw_labels, dtype=int).reshape(-1)
    if labels.size == 0:
        return np.zeros(0, dtype=int)

    ordered_labels = sorted(
        np.unique(labels).tolist(),
        key=lambda label: (
            int(np.flatnonzero(labels == int(label))[0]),
            int(label),
        ),
    )
    canonical = np.zeros(labels.shape[0], dtype=int)
    for canonical_label, source_label in enumerate(ordered_labels):
        canonical[labels == int(source_label)] = int(canonical_label)
    return canonical


def _build_scipy_cluster_labels_from_tree(
    *,
    cluster_tree: _ScipyWardClusterTree,
    cluster_counts: Iterable[int],
) -> dict[int, np.ndarray]:
    requested_counts = [int(cluster_count) for cluster_count in cluster_counts]
    if not requested_counts:
        return {}
    n_sites = int(cluster_tree.n_sites)
    unique_counts = sorted(
        {max(1, min(int(count), n_sites)) for count in requested_counts}
    )
    labels_by_count: dict[int, np.ndarray] = {}

    if n_sites == 0:
        for count in unique_counts:
            labels_by_count[count] = np.zeros(0, dtype=int)
        return {
            count: labels_by_count[max(1, min(int(count), n_sites))]
            for count in requested_counts
        }

    if n_sites == 1:
        for count in unique_counts:
            labels_by_count[count] = np.zeros(1, dtype=int)
        return {
            count: labels_by_count[max(1, min(int(count), n_sites))]
            for count in requested_counts
        }

    for count in unique_counts:
        if count == 1:
            labels_by_count[count] = np.zeros(n_sites, dtype=int)
            continue
        if count == n_sites:
            labels_by_count[count] = np.arange(n_sites, dtype=int)
            continue
        raw_labels = np.asarray(
            cut_tree(cluster_tree.linkage_matrix, n_clusters=[int(count)]),
            dtype=int,
        ).reshape(-1)
        labels_by_count[count] = _canonicalize_labels(raw_labels)

    return {
        count: labels_by_count[max(1, min(int(count), n_sites))].copy()
        for count in requested_counts
    }


@dataclass(frozen=True, slots=True)
class _ScipyClusterTreeOperations:
    def build_cluster_tree(self, scoring_values: np.ndarray) -> _ScipyWardClusterTree:
        return _build_scipy_cluster_tree(scoring_values)

    def build_cluster_labels_from_tree(
        self,
        *,
        cluster_tree: object,
        cluster_counts: Iterable[int],
    ) -> dict[int, np.ndarray]:
        if not isinstance(cluster_tree, _ScipyWardClusterTree):
            raise TypeError(
                "scipy cluster-tree operations expected a _ScipyWardClusterTree instance"
            )
        return _build_scipy_cluster_labels_from_tree(
            cluster_tree=cluster_tree,
            cluster_counts=cluster_counts,
        )


_SCIPY_CLUSTER_TREE_OPERATIONS = _ScipyClusterTreeOperations()


@dataclass(frozen=True, slots=True)
class ScipyHierarchicalClusteringBackend:
    """SciPy-backed hierarchical clustering backend."""

    name: str = SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL
    version: str = SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL_VERSION

    def run(
        self,
        request: SignalomeClusteringBackendRequest,
    ) -> SignalomeClusteringBackendResult:
        if request.cluster_tree_backend != _exact.SIGNALOME_CLUSTER_TREE_BACKEND_EXACT:
            raise ValueError("cluster_tree_backend must be 'exact'")
        if request.candidate_scoring_backend not in {
            None,
            _exact.SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
            _exact.SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
        }:
            raise ValueError("candidate_scoring_backend must be one of: full, sampled")
        cluster_tree_backend = cast(
            _exact.SignalomeClusterTreeBackend,
            request.cluster_tree_backend,
        )
        candidate_scoring_backend = cast(
            _exact.SignalomeCandidateScoringBackend | None,
            request.candidate_scoring_backend,
        )
        clustering_result = _exact.cluster_sites_with_diagnostics(
            scoring_matrix=request.scoring_matrix,
            requested_module_count=request.requested_module_count,
            primary_threshold=request.primary_threshold,
            fallback_threshold=request.fallback_threshold,
            max_clusters=request.max_clusters,
            cluster_tree_backend=cluster_tree_backend,
            candidate_scoring_backend=candidate_scoring_backend,
            max_exact_cluster_tree_sites=request.max_exact_cluster_tree_sites,
            max_full_correlation_sites=request.max_full_correlation_sites,
            cluster_tree_operations=_SCIPY_CLUSTER_TREE_OPERATIONS,
        )
        protein_modules = _exact.derive_protein_modules(
            site_clusters=clustering_result.site_clusters,
            site_to_protein=request.site_to_protein,
        )

        selected_module_count = int(
            clustering_result.module_selection_diagnostics.selected_module_count
        )
        backend_diagnostics: dict[str, object] = {
            "backend_name": self.name,
            "uses_scipy": True,
            "linkage_method": _SCIPY_LINKAGE_METHOD,
            "distance_metric": _SCIPY_DISTANCE_METRIC,
            "selected_module_count": selected_module_count,
            "input_site_count": int(request.scoring_matrix.shape[0]),
            "exact_tree_path_used": bool(clustering_result.exact_cluster_tree_built),
        }
        return SignalomeClusteringBackendResult(
            site_clusters=clustering_result.site_clusters,
            protein_modules=protein_modules,
            selected_module_count=selected_module_count,
            module_selection_diagnostics=clustering_result.module_selection_diagnostics,
            backend_name=self.name,
            backend_version=self.version,
            approximation_used=bool(clustering_result.approximation_used),
            exact_cluster_tree_built=bool(clustering_result.exact_cluster_tree_built),
            cluster_tree_backend=str(clustering_result.cluster_tree_backend),
            candidate_scoring_mode=str(clustering_result.candidate_scoring_mode),
            candidate_scoring_evaluated=bool(
                clustering_result.candidate_scoring_evaluated
            ),
            candidate_scoring_skip_reason=(
                None
                if clustering_result.candidate_scoring_skip_reason is None
                else str(clustering_result.candidate_scoring_skip_reason)
            ),
            candidate_scoring_sampling=clustering_result.candidate_scoring_sampling,
            backend_diagnostics=backend_diagnostics,
            threshold_metadata={
                "primary_threshold": float(request.primary_threshold),
                "fallback_threshold": float(request.fallback_threshold),
            },
            limit_metadata={
                "max_exact_cluster_tree_sites": (
                    None
                    if request.max_exact_cluster_tree_sites is None
                    else int(request.max_exact_cluster_tree_sites)
                ),
                "max_full_correlation_sites": int(request.max_full_correlation_sites),
                "max_clusters": int(request.max_clusters),
            },
        )


__all__ = ["ScipyHierarchicalClusteringBackend"]
