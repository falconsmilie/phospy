from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, linkage
from sklearn.cluster import AgglomerativeClustering

from ..internal.types import SignalomeModuleSelectionStrategy
from ..validation.values.enums import validate_module_selection_strategy
from ..validation.values.numeric import validate_fraction, validate_positive_int

__all__ = [
    "ClusterCandidateScore",
    "ClusterSitesResult",
    "DEFAULT_SIGNALOME_MODULE_SELECTION_POLICY",
    "SignalomeModuleSelectionDiagnostics",
    "SignalomeModuleSelectionPolicy",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "select_module_count",
    "select_module_count_with_diagnostics",
]

MAX_FULL_CORRELATION_SITE_COUNT = 2000
MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER = 256


@dataclass(frozen=True, slots=True)
class ClusterCandidateScore:
    """Cached score summary for one candidate module count."""

    min_median_correlation: float
    mean_median_correlation: float


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionPolicy:
    """Explicit policy for automatic signalome module-count selection.

    ``strategy`` controls how the automatic selector behaves:

    - ``"correlation_thresholds"`` applies the current PhosPy correlation-based
      heuristic using the configured primary and fallback thresholds.
    - ``"single_module"`` bypasses automatic selection and forces one module
      unless the caller explicitly requests a module count.
    """

    strategy: SignalomeModuleSelectionStrategy = "correlation_thresholds"
    primary_threshold: float = 0.5
    fallback_threshold: float = 0.1
    max_clusters: int = 10

    def __post_init__(self) -> None:
        validate_module_selection_strategy(self.strategy)
        validate_fraction(self.primary_threshold, name="primary_threshold")
        validate_fraction(self.fallback_threshold, name="fallback_threshold")
        validate_positive_int(self.max_clusters, name="max_clusters")

    @classmethod
    def from_value(cls, value: object) -> SignalomeModuleSelectionPolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls(**dict(value))
        msg = "module_selection_policy must be a SignalomeModuleSelectionPolicy or mapping"
        raise TypeError(msg)


DEFAULT_SIGNALOME_MODULE_SELECTION_POLICY = SignalomeModuleSelectionPolicy()


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionDiagnostics:
    """Structured explanation of how a signalome module count was chosen."""

    strategy: SignalomeModuleSelectionStrategy
    selected_module_count: int
    requested_module_count: int | None
    threshold_used: float | None
    max_clusters_evaluated: int
    candidate_scores: dict[int, ClusterCandidateScore]
    reason: str

    @property
    def used_automatic_selection(self) -> bool:
        return self.requested_module_count is None


@dataclass(frozen=True, slots=True)
class ClusterSitesResult:
    """Cluster labels plus module-selection diagnostics."""

    site_clusters: pd.Series
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> pd.Series:
    """Cluster phosphosites into signalome site clusters."""

    return cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        policy=policy,
    ).site_clusters


def cluster_sites_with_diagnostics(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> ClusterSitesResult:
    """Cluster phosphosites and capture how the module count was chosen."""

    n_sites = scoring_matrix.shape[0]
    if n_sites == 1:
        diagnostics = SignalomeModuleSelectionDiagnostics(
            strategy=(
                DEFAULT_SIGNALOME_MODULE_SELECTION_POLICY.strategy
                if policy is None
                else policy.strategy
            ),
            selected_module_count=1,
            requested_module_count=requested_module_count,
            threshold_used=None,
            max_clusters_evaluated=1,
            candidate_scores={},
            reason="single phosphosite input only supports one signalome module",
        )
        return ClusterSitesResult(
            site_clusters=pd.Series(
                [1], index=scoring_matrix.index, dtype=int, name="site_cluster"
            ),
            module_selection_diagnostics=diagnostics,
        )

    scoring_values = scoring_matrix.to_numpy(dtype=float)
    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        policy=policy,
    )
    module_count = max(1, min(diagnostics.selected_module_count, n_sites))

    if module_count == 1:
        labels = np.ones(n_sites, dtype=int)
    else:
        labels = fit_cluster_labels(scoring_values, module_count) + 1

    return ClusterSitesResult(
        site_clusters=pd.Series(
            labels,
            index=scoring_matrix.index,
            dtype=int,
            name="site_cluster",
        ),
        module_selection_diagnostics=diagnostics,
    )


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> int:
    """Choose a signalome module count from the scoring matrix."""

    return select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        policy=policy,
    ).selected_module_count


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> SignalomeModuleSelectionDiagnostics:
    """Choose a module count and explain why that count was selected."""

    resolved_policy = SignalomeModuleSelectionPolicy.from_value(policy)
    scoring_array = np.asarray(scoring_values, dtype=float)
    n_sites = scoring_array.shape[0]
    if n_sites <= 1:
        return SignalomeModuleSelectionDiagnostics(
            strategy=resolved_policy.strategy,
            selected_module_count=1,
            requested_module_count=requested_module_count,
            threshold_used=None,
            max_clusters_evaluated=1,
            candidate_scores={},
            reason="single phosphosite input only supports one signalome module",
        )

    if requested_module_count is not None:
        resolved_count = max(1, min(int(requested_module_count), n_sites))
        return SignalomeModuleSelectionDiagnostics(
            strategy=resolved_policy.strategy,
            selected_module_count=resolved_count,
            requested_module_count=int(requested_module_count),
            threshold_used=None,
            max_clusters_evaluated=min(resolved_policy.max_clusters, n_sites),
            candidate_scores={},
            reason="module_count was provided explicitly by the caller",
        )

    if resolved_policy.strategy == "single_module":
        return SignalomeModuleSelectionDiagnostics(
            strategy=resolved_policy.strategy,
            selected_module_count=1,
            requested_module_count=None,
            threshold_used=None,
            max_clusters_evaluated=1,
            candidate_scores={},
            reason="module_selection_strategy='single_module' forces one module",
        )

    max_clusters = min(resolved_policy.max_clusters, n_sites)
    if max_clusters < 2:
        return SignalomeModuleSelectionDiagnostics(
            strategy=resolved_policy.strategy,
            selected_module_count=1,
            requested_module_count=None,
            threshold_used=None,
            max_clusters_evaluated=max_clusters,
            candidate_scores={},
            reason="fewer than two cluster counts are available for evaluation",
        )

    candidate_range = range(2, max_clusters + 1)
    approximation_note = ""
    if n_sites <= MAX_FULL_CORRELATION_SITE_COUNT:
        site_correlations = np.corrcoef(scoring_array)
        candidate_scores = score_cluster_candidates(
            scoring_values=scoring_array,
            site_correlations=site_correlations,
            cluster_range=candidate_range,
        )
    else:
        candidate_scores = score_cluster_candidates_approximate(
            scoring_values=scoring_array,
            cluster_range=candidate_range,
            max_sites_per_cluster=MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
        )
        approximation_note = (
            " Used sampled within-cluster correlation estimates to avoid "
            "materializing a full site-by-site correlation matrix."
        )
    primary_candidates = filter_cluster_candidates(
        candidate_scores,
        threshold=resolved_policy.primary_threshold,
    )
    if primary_candidates:
        selected_count = max(
            primary_candidates.items(),
            key=lambda item: (item[1], -item[0]),
        )[0]
        return SignalomeModuleSelectionDiagnostics(
            strategy=resolved_policy.strategy,
            selected_module_count=selected_count,
            requested_module_count=None,
            threshold_used=resolved_policy.primary_threshold,
            max_clusters_evaluated=max_clusters,
            candidate_scores=candidate_scores,
            reason=(
                "selected the highest-scoring candidate that satisfied the "
                "primary within-cluster correlation threshold"
            )
            + approximation_note,
        )

    fallback_candidates = filter_cluster_candidates(
        candidate_scores,
        threshold=resolved_policy.fallback_threshold,
    )
    if fallback_candidates:
        selected_count = max(
            fallback_candidates.items(),
            key=lambda item: (item[1], -item[0]),
        )[0]
        return SignalomeModuleSelectionDiagnostics(
            strategy=resolved_policy.strategy,
            selected_module_count=selected_count,
            requested_module_count=None,
            threshold_used=resolved_policy.fallback_threshold,
            max_clusters_evaluated=max_clusters,
            candidate_scores=candidate_scores,
            reason=(
                "no candidate satisfied the primary threshold; selected the "
                "highest-scoring fallback candidate"
            )
            + approximation_note,
        )

    return SignalomeModuleSelectionDiagnostics(
        strategy=resolved_policy.strategy,
        selected_module_count=1,
        requested_module_count=None,
        threshold_used=None,
        max_clusters_evaluated=max_clusters,
        candidate_scores=candidate_scores,
        reason=(
            "no candidate module count satisfied the configured correlation "
            "thresholds, so the workflow fell back to one module"
        )
        + approximation_note,
    )


def filter_cluster_candidates(
    candidate_scores: dict[int, ClusterCandidateScore],
    *,
    threshold: float,
) -> dict[int, float]:
    """Return candidate counts whose cluster medians satisfy a threshold."""

    return {
        cluster_count: score.mean_median_correlation
        for cluster_count, score in candidate_scores.items()
        if score.min_median_correlation >= threshold
    }


def score_cluster_candidates(
    *,
    scoring_values: np.ndarray,
    site_correlations: np.ndarray,
    cluster_range: Iterable[int],
) -> dict[int, ClusterCandidateScore]:
    """Score candidate module counts using one cached Ward hierarchy."""

    cluster_counts = [int(cluster_count) for cluster_count in cluster_range]
    if not cluster_counts:
        return {}

    linkage_matrix = build_cluster_tree(scoring_values)
    candidate_labels = build_cluster_labels_from_tree(
        linkage_matrix=linkage_matrix,
        cluster_counts=cluster_counts,
    )

    candidates: dict[int, ClusterCandidateScore] = {}
    for cluster_count in cluster_counts:
        labels = candidate_labels[cluster_count]
        cluster_medians = [
            cluster_median_correlation(site_correlations, labels, label)
            for label in np.unique(labels)
        ]
        if not cluster_medians:
            continue
        candidates[cluster_count] = ClusterCandidateScore(
            min_median_correlation=float(min(cluster_medians)),
            mean_median_correlation=float(np.mean(cluster_medians)),
        )
    return candidates


def score_cluster_candidates_approximate(
    *,
    scoring_values: np.ndarray,
    cluster_range: Iterable[int],
    max_sites_per_cluster: int,
) -> dict[int, ClusterCandidateScore]:
    """Score candidate counts using sampled cluster-local correlations.

    This avoids materializing a full site-by-site correlation matrix for very
    large phosphosite sets.
    """

    cluster_counts = [int(cluster_count) for cluster_count in cluster_range]
    if not cluster_counts:
        return {}

    linkage_matrix = build_cluster_tree(scoring_values)
    candidate_labels = build_cluster_labels_from_tree(
        linkage_matrix=linkage_matrix,
        cluster_counts=cluster_counts,
    )

    candidates: dict[int, ClusterCandidateScore] = {}
    for cluster_count in cluster_counts:
        labels = candidate_labels[cluster_count]
        cluster_medians = [
            cluster_median_correlation_approximate(
                scoring_values=scoring_values,
                labels=labels,
                label=label,
                max_sites_per_cluster=max_sites_per_cluster,
            )
            for label in np.unique(labels)
        ]
        if not cluster_medians:
            continue
        candidates[cluster_count] = ClusterCandidateScore(
            min_median_correlation=float(min(cluster_medians)),
            mean_median_correlation=float(np.mean(cluster_medians)),
        )
    return candidates


def build_cluster_tree(scoring_values: np.ndarray) -> np.ndarray:
    """Build one Ward hierarchical tree for candidate module evaluation."""

    return linkage(np.asarray(scoring_values, dtype=float), method="ward")


def build_cluster_labels_from_tree(
    *,
    linkage_matrix: np.ndarray,
    cluster_counts: Iterable[int],
) -> dict[int, np.ndarray]:
    """Cut one cached hierarchy into labels for each candidate module count."""

    cluster_count_list = [int(cluster_count) for cluster_count in cluster_counts]
    if not cluster_count_list:
        return {}

    cut_labels = cut_tree(linkage_matrix, n_clusters=cluster_count_list)
    if cut_labels.ndim == 1:
        cut_labels = cut_labels.reshape(-1, 1)

    return {
        cluster_count: cut_labels[:, position].astype(int, copy=False)
        for position, cluster_count in enumerate(cluster_count_list)
    }


def fit_cluster_labels(scoring_values: np.ndarray, cluster_count: int) -> np.ndarray:
    """Fit Ward agglomerative clustering once for one candidate count."""

    return (
        AgglomerativeClustering(
            n_clusters=cluster_count,
            linkage="ward",
        )
        .fit_predict(scoring_values)
        .astype(int)
    )


def cluster_median_correlation(
    site_correlations: np.ndarray,
    labels: np.ndarray,
    label: int,
) -> float:
    """Return the median within-cluster correlation for one cluster label."""

    cluster_positions = np.flatnonzero(labels == label)
    if cluster_positions.size <= 1:
        return 0.0

    cluster_correlations = site_correlations[
        np.ix_(cluster_positions, cluster_positions)
    ]
    cluster_correlations = cluster_correlations.copy()
    np.fill_diagonal(cluster_correlations, np.nan)
    values = cluster_correlations[~np.isnan(cluster_correlations)]
    if values.size == 0:
        return 0.0
    return float(np.median(values))


def cluster_median_correlation_approximate(
    *,
    scoring_values: np.ndarray,
    labels: np.ndarray,
    label: int,
    max_sites_per_cluster: int,
) -> float:
    """Approximate the within-cluster median correlation for one label."""

    cluster_positions = np.flatnonzero(labels == label)
    if cluster_positions.size <= 1:
        return 0.0

    if cluster_positions.size > max_sites_per_cluster:
        sampled_positions = np.linspace(
            0,
            cluster_positions.size - 1,
            num=max_sites_per_cluster,
            dtype=int,
        )
        cluster_positions = cluster_positions[sampled_positions]

    cluster_values = scoring_values[cluster_positions]
    cluster_correlations = np.corrcoef(cluster_values)
    cluster_correlations = np.asarray(cluster_correlations, dtype=float).copy()
    np.fill_diagonal(cluster_correlations, np.nan)
    values = cluster_correlations[~np.isnan(cluster_correlations)]
    if values.size == 0:
        return 0.0
    return float(np.median(values))
