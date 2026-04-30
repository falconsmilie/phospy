"""SciPy hierarchical Ward tree engine for signalome clustering."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import pdist

from phospy.signalomes.clustering.diagnostic_schemas import (
    SignalomeTreeEngineDiagnostics,
    build_tree_engine_diagnostics,
)
from phospy.signalomes.clustering.validation import (
    validate_cluster_count_for_site_count,
)

_SCIPY_LINKAGE_METHOD = "ward"
_SCIPY_DISTANCE_METRIC = "euclidean"


@dataclass(frozen=True, slots=True)
class ScipyWardClusterTree:
    n_sites: int
    linkage_matrix: np.ndarray


def build_cluster_tree(scoring_values: np.ndarray) -> ScipyWardClusterTree:
    values = np.asarray(scoring_values, dtype=float)
    if values.ndim != 2:
        raise ValueError("scoring_values must be a 2D array")
    n_sites = int(values.shape[0])
    if n_sites <= 1:
        return ScipyWardClusterTree(
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
    return ScipyWardClusterTree(n_sites=n_sites, linkage_matrix=linkage_matrix)


def build_cluster_labels_from_tree(
    *,
    cluster_tree: ScipyWardClusterTree,
    cluster_counts: Iterable[int],
) -> dict[int, np.ndarray]:
    requested_counts = [int(cluster_count) for cluster_count in cluster_counts]
    if not requested_counts:
        return {}
    n_sites = int(cluster_tree.n_sites)
    unique_counts = sorted(
        {
            validate_cluster_count_for_site_count(
                cluster_count=int(count),
                available_clustering_site_count=n_sites,
                field_name="cluster_counts",
            )
            for count in requested_counts
        }
    )
    labels_by_count: dict[int, np.ndarray] = {}

    if n_sites == 0:
        for count in unique_counts:
            labels_by_count[count] = np.zeros(0, dtype=int)
        return {count: labels_by_count[int(count)] for count in requested_counts}

    if n_sites == 1:
        for count in unique_counts:
            labels_by_count[count] = np.zeros(1, dtype=int)
        return {count: labels_by_count[int(count)] for count in requested_counts}

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

    return {count: labels_by_count[int(count)].copy() for count in requested_counts}


def engine_diagnostics() -> SignalomeTreeEngineDiagnostics:
    return build_tree_engine_diagnostics(
        uses_scipy=True,
        linkage_method=_SCIPY_LINKAGE_METHOD,
        distance_metric=_SCIPY_DISTANCE_METRIC,
    )


@dataclass(frozen=True, slots=True)
class ScipyHierarchicalTreeEngine:
    """SciPy Ward hierarchical tree engine."""

    name: str = "scipy_hierarchical_tree"
    version: str = "1"

    def build_tree(self, values: np.ndarray) -> ScipyWardClusterTree:
        return build_cluster_tree(values)

    def labels_for_counts(
        self,
        *,
        tree: object,
        counts: Iterable[int],
    ) -> dict[int, np.ndarray]:
        if not isinstance(tree, ScipyWardClusterTree):
            raise TypeError(
                "scipy tree engine expected a ScipyWardClusterTree instance"
            )
        return build_cluster_labels_from_tree(cluster_tree=tree, cluster_counts=counts)


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


__all__ = [
    "ScipyHierarchicalTreeEngine",
    "ScipyWardClusterTree",
    "build_cluster_labels_from_tree",
    "build_cluster_tree",
    "engine_diagnostics",
]
