"""Backend-independent cluster tree construction helpers."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import cast

import numpy as np

from phospy.signalomes.clustering.backends.exact_python import (
    ExactWardClusterTree,
)
from phospy.signalomes.clustering.backends.exact_python import (
    build_cluster_labels_from_tree as build_exact_cluster_labels_from_tree,
)
from phospy.signalomes.clustering.backends.exact_python import (
    build_cluster_tree as build_exact_cluster_tree,
)
from phospy.signalomes.clustering.contracts import ClusterTreeEngine
from phospy.signalomes.clustering.policies import (
    MAX_FULL_CORRELATION_SITE_COUNT,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeTreeEngine,
)
from phospy.signalomes.clustering.scale_guards import (
    raise_if_exact_tree_limit_exceeded,
    resolve_max_exact_tree_sites,
)
from phospy.signalomes.clustering.validation import (
    validate_cluster_count_for_site_count,
)


@dataclass(frozen=True, slots=True)
class _ClusterTreeOperationsAdapter:
    """Compatibility adapter for legacy cluster-tree operation hooks."""

    engine: ClusterTreeEngine

    def build_cluster_tree(self, scoring_values: np.ndarray) -> object:
        return self.engine.build_tree(scoring_values)

    def build_cluster_labels_from_tree(
        self,
        *,
        cluster_tree: object,
        cluster_counts: Iterable[int],
    ) -> dict[int, np.ndarray]:
        return self.engine.labels_for_counts(tree=cluster_tree, counts=cluster_counts)


@dataclass(frozen=True, slots=True)
class _ExactWardClusterTreeOperations:
    """Default exact-tree operations with legacy monkeypatch hook behavior."""

    def build_cluster_tree(self, scoring_values: np.ndarray) -> object:
        # Keep a runtime indirection through the exact-python module so
        # performance-contract monkeypatching of exact-python helpers remains
        # effective after orchestration extraction.
        from phospy.signalomes.clustering import exact_python as exact_clustering

        return exact_clustering._build_cluster_tree(scoring_values)

    def build_cluster_labels_from_tree(
        self,
        *,
        cluster_tree: object,
        cluster_counts: Iterable[int],
    ) -> dict[int, np.ndarray]:
        from phospy.signalomes.clustering import exact_python as exact_clustering

        return exact_clustering.build_cluster_labels_from_tree(
            cluster_tree=cast(ExactWardClusterTree, cluster_tree),
            cluster_counts=cluster_counts,
        )


ClusterTreeOperations = _ClusterTreeOperationsAdapter
ClusterTreeOperationsAdapter = _ClusterTreeOperationsAdapter

_EXACT_WARD_CLUSTER_TREE_OPERATIONS = _ExactWardClusterTreeOperations()


def resolve_cluster_tree_operations(
    cluster_tree_operations: ClusterTreeOperations | None,
) -> ClusterTreeOperations | _ExactWardClusterTreeOperations:
    if cluster_tree_operations is None:
        return _EXACT_WARD_CLUSTER_TREE_OPERATIONS
    return cluster_tree_operations


def build_cluster_tree(scoring_values: np.ndarray) -> ExactWardClusterTree:
    """Compatibility wrapper for pure-Python exact Ward tree construction."""

    return build_exact_cluster_tree(scoring_values)


def build_cluster_labels_from_tree(
    *,
    cluster_tree: ExactWardClusterTree,
    cluster_counts: Iterable[int],
) -> dict[int, np.ndarray]:
    """Compatibility wrapper for pure-Python exact Ward tree cutting."""

    return build_exact_cluster_labels_from_tree(
        cluster_tree=cluster_tree,
        cluster_counts=cluster_counts,
    )


def build_exact_cluster_tree_with_guard(
    *,
    clustering_values: np.ndarray,
    n_sites: int,
    tree_engine: SignalomeTreeEngine,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy,
    max_exact_tree_sites: int | None,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> object:
    if tree_engine != SIGNALOME_TREE_ENGINE_EXACT:
        raise ValueError("tree_engine must be 'exact'")
    raise_if_exact_tree_limit_exceeded(
        n_sites=n_sites,
        max_exact_tree_sites=max_exact_tree_sites,
        candidate_scoring_policy=candidate_scoring_policy,
    )
    tree_operations = resolve_cluster_tree_operations(cluster_tree_operations)
    return tree_operations.build_cluster_tree(clustering_values)


def fit_cluster_labels(
    scoring_values: np.ndarray,
    cluster_count: int,
    *,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy = (
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    ),
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> np.ndarray:
    """Fit Ward clustering and return 0-indexed labels for one count."""

    values = np.asarray(scoring_values, dtype=float)
    n_sites = int(values.shape[0])
    resolved_cluster_count = validate_cluster_count_for_site_count(
        cluster_count=int(cluster_count),
        available_clustering_site_count=n_sites,
        field_name="cluster_count",
    )
    if resolved_cluster_count == 1:
        return np.zeros(n_sites, dtype=int)
    resolved_max_exact_tree_sites = resolve_max_exact_tree_sites(max_exact_tree_sites)
    tree = build_exact_cluster_tree_with_guard(
        clustering_values=values,
        n_sites=n_sites,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=resolved_max_exact_tree_sites,
        cluster_tree_operations=cluster_tree_operations,
    )
    tree_operations = resolve_cluster_tree_operations(cluster_tree_operations)
    return tree_operations.build_cluster_labels_from_tree(
        cluster_tree=tree,
        cluster_counts=[resolved_cluster_count],
    )[resolved_cluster_count].astype(int, copy=False)


def prepare_scoring_values_for_clustering(scoring_values: np.ndarray) -> np.ndarray:
    values = np.asarray(scoring_values, dtype=float).copy()
    if values.size == 0:
        return values
    values[~np.isfinite(values)] = np.nan
    column_medians = np.nanmedian(values, axis=0)
    column_medians = np.where(np.isfinite(column_medians), column_medians, 0.0)
    row_positions, column_positions = np.where(np.isnan(values))
    if row_positions.size > 0:
        values[row_positions, column_positions] = column_medians[column_positions]
    return values


__all__ = [
    "ClusterTreeOperations",
    "ClusterTreeOperationsAdapter",
    "build_cluster_labels_from_tree",
    "build_cluster_tree",
    "build_exact_cluster_tree_with_guard",
    "fit_cluster_labels",
    "prepare_scoring_values_for_clustering",
    "resolve_cluster_tree_operations",
]
