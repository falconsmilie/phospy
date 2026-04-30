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
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringMissingValuePolicy,
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


@dataclass(frozen=True, slots=True)
class SignalomeClusteringMissingValueDiagnostics:
    policy: SignalomeClusteringMissingValuePolicy
    non_finite_input_value_count: int
    missing_after_non_finite_normalization_count: int
    imputed_value_count: int
    fully_missing_column_count: int


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


def summarize_clustering_missing_value_diagnostics(
    scoring_values: np.ndarray,
    *,
    missing_value_policy: SignalomeClusteringMissingValuePolicy = (
        SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS
    ),
) -> SignalomeClusteringMissingValueDiagnostics:
    _validate_signalome_clustering_missing_value_policy(missing_value_policy)
    values = np.asarray(scoring_values, dtype=float)
    if values.ndim != 2:
        raise ValueError("scoring_values must be a 2D array")
    if values.size == 0:
        return SignalomeClusteringMissingValueDiagnostics(
            policy=missing_value_policy,
            non_finite_input_value_count=0,
            missing_after_non_finite_normalization_count=0,
            imputed_value_count=0,
            fully_missing_column_count=0,
        )
    non_finite_mask = ~np.isfinite(values)
    normalized_values = values.copy()
    normalized_values[non_finite_mask] = np.nan
    missing_mask = np.isnan(normalized_values)
    fully_missing_columns = np.all(missing_mask, axis=0)
    missing_count = int(np.count_nonzero(missing_mask))
    return SignalomeClusteringMissingValueDiagnostics(
        policy=missing_value_policy,
        non_finite_input_value_count=int(np.count_nonzero(non_finite_mask)),
        missing_after_non_finite_normalization_count=missing_count,
        imputed_value_count=missing_count,
        fully_missing_column_count=int(np.count_nonzero(fully_missing_columns)),
    )


def prepare_scoring_values_for_clustering(
    scoring_values: np.ndarray,
    *,
    missing_value_policy: SignalomeClusteringMissingValuePolicy = (
        SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS
    ),
) -> np.ndarray:
    _validate_signalome_clustering_missing_value_policy(missing_value_policy)
    values = np.asarray(scoring_values, dtype=float).copy()
    if values.ndim != 2:
        raise ValueError("scoring_values must be a 2D array")
    if values.size == 0:
        return values
    values[~np.isfinite(values)] = np.nan
    column_medians = _column_medians_with_zero_for_all_missing_columns(values)
    row_positions, column_positions = np.where(np.isnan(values))
    if row_positions.size > 0:
        values[row_positions, column_positions] = column_medians[column_positions]
    return values


def _validate_signalome_clustering_missing_value_policy(
    missing_value_policy: SignalomeClusteringMissingValuePolicy,
) -> None:
    if (
        missing_value_policy
        != SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS
    ):
        raise ValueError(
            "unsupported signalome clustering missing-value policy: "
            f"{missing_value_policy!r}"
        )


def _column_medians_with_zero_for_all_missing_columns(
    values: np.ndarray,
) -> np.ndarray:
    n_columns = int(values.shape[1])
    column_medians = np.zeros(n_columns, dtype=float)
    for column_index in range(n_columns):
        column_values = values[:, column_index]
        finite_values = column_values[np.isfinite(column_values)]
        if finite_values.size == 0:
            column_medians[column_index] = 0.0
            continue
        column_medians[column_index] = float(np.median(finite_values))
    return column_medians


__all__ = [
    "ClusterTreeOperations",
    "ClusterTreeOperationsAdapter",
    "SignalomeClusteringMissingValueDiagnostics",
    "build_cluster_labels_from_tree",
    "build_cluster_tree",
    "build_exact_cluster_tree_with_guard",
    "fit_cluster_labels",
    "prepare_scoring_values_for_clustering",
    "resolve_cluster_tree_operations",
    "summarize_clustering_missing_value_diagnostics",
]
