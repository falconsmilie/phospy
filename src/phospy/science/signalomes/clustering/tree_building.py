"""Backend-independent cluster tree construction helpers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
import pandas as pd

from phospy.provenance.hashing import fingerprint_matrix
from phospy.provenance.models import TableFingerprint
from phospy.science.signalomes.clustering.backends import (
    exact_python as exact_tree_backend,
)
from phospy.science.signalomes.clustering.backends.exact_python import (
    ExactWardClusterTree,
)
from phospy.science.signalomes.clustering.backends.exact_python import (
    build_cluster_labels_from_tree as build_exact_cluster_labels_from_tree,
)
from phospy.science.signalomes.clustering.backends.exact_python import (
    build_cluster_tree as build_exact_cluster_tree,
)
from phospy.science.signalomes.clustering.contracts import ClusterTreeEngine
from phospy.science.signalomes.clustering.policies import (
    MAX_FULL_CORRELATION_SITE_COUNT,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringMissingValuePolicy,
    SignalomeTreeEngine,
)
from phospy.science.signalomes.clustering.scale_guards import (
    raise_if_exact_tree_limit_exceeded,
    resolve_max_exact_tree_sites,
)
from phospy.science.signalomes.clustering.validation import (
    validate_cluster_count_for_site_count,
)
from phospy.science.signalomes.models import SignalomeClusteringPreparationDiagnostics

SIGNALOME_CLUSTERING_PREPARATION_POLICY_ID = SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE
SIGNALOME_CLUSTERING_PREPARED_MATRIX_FINGERPRINT_NAME = (
    "signalome.clustering.prepared_matrix"
)


@dataclass(frozen=True, slots=True)
class _ClusterTreeOperationsAdapter:
    """Adapter that maps ClusterTreeEngine to tree operations."""

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
class ExactWardClusterTreeBuilder:
    """Exact Ward tree builder and label cutter."""

    def build_cluster_tree(self, scoring_values: np.ndarray) -> object:
        return exact_tree_backend.build_cluster_tree(scoring_values)

    def build_cluster_labels_from_tree(
        self,
        *,
        cluster_tree: object,
        cluster_counts: Iterable[int],
    ) -> dict[int, np.ndarray]:
        if not isinstance(cluster_tree, ExactWardClusterTree):
            raise TypeError("cluster_tree must be an ExactWardClusterTree instance")
        return exact_tree_backend.build_cluster_labels_from_tree(
            cluster_tree=cluster_tree,
            cluster_counts=cluster_counts,
        )


SignalomeClusteringMissingValueDiagnostics = SignalomeClusteringPreparationDiagnostics


@dataclass(frozen=True, slots=True)
class PreparedSignalomeClusteringMatrix:
    """Labeled, finite signalome clustering matrix plus preparation audit."""

    prepared_matrix: pd.DataFrame
    retained_column_labels: tuple[str, ...]
    dropped_fully_missing_column_labels: tuple[str, ...]
    imputed_value_count: int
    imputed_value_counts_by_column: Mapping[str, int]
    dropped_fully_missing_cell_count: int
    non_finite_input_value_count: int
    missing_after_non_finite_normalization_count: int
    prepared_matrix_fingerprint: TableFingerprint
    preparation_policy_id: str = SIGNALOME_CLUSTERING_PREPARATION_POLICY_ID

    def __post_init__(self) -> None:
        retained = tuple(str(value) for value in self.retained_column_labels)
        dropped = tuple(
            str(value) for value in self.dropped_fully_missing_column_labels
        )
        imputation_counts = {
            str(key): int(value)
            for key, value in self.imputed_value_counts_by_column.items()
        }
        object.__setattr__(self, "retained_column_labels", retained)
        object.__setattr__(self, "dropped_fully_missing_column_labels", dropped)
        object.__setattr__(
            self,
            "imputed_value_counts_by_column",
            MappingProxyType(imputation_counts),
        )

    @property
    def values(self) -> np.ndarray:
        """Return finite prepared values for tree backends."""

        return self.prepared_matrix.to_numpy(dtype=float, copy=True)

    def to_diagnostics(self) -> SignalomeClusteringPreparationDiagnostics:
        """Return the public diagnostics view of this preparation audit."""

        return SignalomeClusteringPreparationDiagnostics(
            preparation_policy_id=self.preparation_policy_id,
            input_dimension_count=(
                len(self.retained_column_labels)
                + len(self.dropped_fully_missing_column_labels)
            ),
            retained_dimension_count=len(self.retained_column_labels),
            retained_dimension_labels=self.retained_column_labels,
            dropped_fully_missing_dimension_count=(
                len(self.dropped_fully_missing_column_labels)
            ),
            dropped_fully_missing_dimension_labels=(
                self.dropped_fully_missing_column_labels
            ),
            dropped_fully_missing_dimension_preview=(
                self.dropped_fully_missing_column_labels[:5]
            ),
            dropped_fully_missing_value_count=int(
                self.dropped_fully_missing_cell_count
            ),
            non_finite_input_value_count=int(self.non_finite_input_value_count),
            missing_after_non_finite_normalization_count=int(
                self.missing_after_non_finite_normalization_count
            ),
            imputed_value_count=int(self.imputed_value_count),
            imputed_value_counts_by_dimension=dict(self.imputed_value_counts_by_column),
            prepared_matrix_fingerprint=self.prepared_matrix_fingerprint,
        )


ClusterTreeOperations = _ClusterTreeOperationsAdapter
ClusterTreeOperationsAdapter = _ClusterTreeOperationsAdapter

_EXACT_WARD_CLUSTER_TREE_OPERATIONS = ExactWardClusterTreeBuilder()


def resolve_cluster_tree_operations(
    cluster_tree_operations: ClusterTreeOperations | None,
) -> ClusterTreeOperations | ExactWardClusterTreeBuilder:
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
        SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE
    ),
) -> SignalomeClusteringPreparationDiagnostics:
    return prepare_signalome_clustering_matrix(
        scoring_values,
        missing_value_policy=missing_value_policy,
    ).to_diagnostics()


def prepare_scoring_values_for_clustering(
    scoring_values: np.ndarray,
    *,
    missing_value_policy: SignalomeClusteringMissingValuePolicy = (
        SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE
    ),
) -> np.ndarray:
    return prepare_signalome_clustering_matrix(
        scoring_values,
        missing_value_policy=missing_value_policy,
    ).values


def prepare_signalome_clustering_matrix(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    missing_value_policy: SignalomeClusteringMissingValuePolicy = (
        SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE
    ),
) -> PreparedSignalomeClusteringMatrix:
    _validate_signalome_clustering_missing_value_policy(missing_value_policy)
    values = _coerce_scoring_values_to_dataframe(scoring_values)
    if values.ndim != 2:
        raise ValueError("scoring_values must be a 2D matrix")

    numeric_values = values.astype(float).copy(deep=True)
    input_array = numeric_values.to_numpy(dtype=float, copy=True)
    non_finite_mask = ~np.isfinite(input_array)
    non_finite_frame = pd.DataFrame(
        non_finite_mask,
        index=numeric_values.index.copy(),
        columns=numeric_values.columns.copy(),
    )
    normalized = numeric_values.mask(non_finite_frame, other=np.nan)
    missing_mask = normalized.isna()
    fully_missing_mask = missing_mask.all(axis=0)
    dropped_labels = tuple(
        str(label)
        for label in normalized.columns[fully_missing_mask.to_numpy(dtype=bool)]
    )
    retained = normalized.loc[:, ~fully_missing_mask].copy(deep=True)
    if retained.shape[1] == 0:
        raise ValueError(
            "signalome clustering preparation retained no kinase/dimension columns "
            "after dropping fully missing dimensions; dropped_fully_missing_dimension_count="
            f"{len(dropped_labels)}; next_action=provide at least one signalome "
            "score dimension with a finite value before clustering"
        )

    retained_missing_mask = retained.isna()
    per_column_imputation_counts = {
        str(column): int(retained_missing_mask.iloc[:, position].sum())
        for position, column in enumerate(retained.columns)
    }
    column_medians = retained.median(axis=0, skipna=True)
    prepared = retained.fillna(column_medians).astype(float)
    prepared.index = _stable_string_index_for_fingerprint(values.index)
    prepared.columns = pd.Index(
        [str(label) for label in prepared.columns],
        name=values.columns.name,
        dtype=object,
    )
    prepared_array = prepared.to_numpy(dtype=float, copy=True)
    if not np.isfinite(prepared_array).all():
        raise ValueError(
            "signalome clustering preparation produced non-finite values; "
            "next_action=inspect retained partially missing dimensions and input "
            "score values"
        )
    fingerprint = fingerprint_matrix(
        prepared,
        name=SIGNALOME_CLUSTERING_PREPARED_MATRIX_FINGERPRINT_NAME,
    )
    return PreparedSignalomeClusteringMatrix(
        prepared_matrix=prepared,
        retained_column_labels=tuple(str(label) for label in prepared.columns),
        dropped_fully_missing_column_labels=dropped_labels,
        imputed_value_count=int(sum(per_column_imputation_counts.values())),
        imputed_value_counts_by_column=per_column_imputation_counts,
        dropped_fully_missing_cell_count=int(len(dropped_labels) * prepared.shape[0]),
        non_finite_input_value_count=int(np.count_nonzero(non_finite_mask)),
        missing_after_non_finite_normalization_count=int(
            missing_mask.to_numpy(dtype=bool, copy=False).sum()
        ),
        prepared_matrix_fingerprint=fingerprint,
        preparation_policy_id=str(missing_value_policy),
    )


def _validate_signalome_clustering_missing_value_policy(
    missing_value_policy: SignalomeClusteringMissingValuePolicy,
) -> None:
    if (
        missing_value_policy
        != SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE
    ):
        raise ValueError(
            "unsupported signalome clustering missing-value policy: "
            f"{missing_value_policy!r}"
        )


def _coerce_scoring_values_to_dataframe(
    scoring_values: pd.DataFrame | np.ndarray,
) -> pd.DataFrame:
    if isinstance(scoring_values, pd.DataFrame):
        return scoring_values.copy(deep=True)
    values = np.array(scoring_values, dtype=float, copy=True)
    if values.ndim != 2:
        raise ValueError("scoring_values must be a 2D matrix")
    return pd.DataFrame(
        values,
        columns=[f"dimension_{position}" for position in range(values.shape[1])],
        copy=True,
    )


def _stable_string_index_for_fingerprint(index: pd.Index) -> pd.Index:
    if all(isinstance(label, str) for label in index.to_numpy(dtype=object)):
        return pd.Index(
            [str(label) for label in index],
            name=index.name,
            dtype=object,
        )
    copied = index.copy()
    copied.name = index.name
    return copied


__all__ = [
    "ClusterTreeOperations",
    "ClusterTreeOperationsAdapter",
    "PreparedSignalomeClusteringMatrix",
    "SIGNALOME_CLUSTERING_PREPARATION_POLICY_ID",
    "SIGNALOME_CLUSTERING_PREPARED_MATRIX_FINGERPRINT_NAME",
    "SignalomeClusteringMissingValueDiagnostics",
    "build_cluster_labels_from_tree",
    "build_cluster_tree",
    "build_exact_cluster_tree_with_guard",
    "fit_cluster_labels",
    "prepare_signalome_clustering_matrix",
    "prepare_scoring_values_for_clustering",
    "resolve_cluster_tree_operations",
    "summarize_clustering_missing_value_diagnostics",
]
