from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, linkage
from sklearn.cluster import AgglomerativeClustering

__all__ = ["cluster_sites"]


@dataclass(frozen=True, slots=True)
class ClusterCandidateScore:
    """Cached score summary for one candidate module count."""

    min_median_correlation: float
    mean_median_correlation: float


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
) -> pd.Series:
    """Cluster phosphosites into signalome site clusters."""

    n_sites = scoring_matrix.shape[0]
    if n_sites == 1:
        return pd.Series(
            [1], index=scoring_matrix.index, dtype=int, name="site_cluster"
        )

    scoring_values = scoring_matrix.to_numpy(dtype=float)
    module_count = (
        requested_module_count
        if requested_module_count is not None
        else select_module_count(scoring_values)
    )
    module_count = max(1, min(module_count, n_sites))

    if module_count == 1:
        labels = np.ones(n_sites, dtype=int)
    else:
        labels = fit_cluster_labels(scoring_values, module_count) + 1

    return pd.Series(labels, index=scoring_matrix.index, dtype=int, name="site_cluster")


def select_module_count(scoring_values: pd.DataFrame | np.ndarray) -> int:
    """Choose a signalome module count from the scoring matrix."""

    scoring_array = np.asarray(scoring_values, dtype=float)
    n_sites = scoring_array.shape[0]
    if n_sites <= 1:
        return 1

    max_clusters = min(10, n_sites)
    if max_clusters < 2:
        return 1

    site_correlations = np.corrcoef(scoring_array)
    candidate_scores = score_cluster_candidates(
        scoring_values=scoring_array,
        site_correlations=site_correlations,
        cluster_range=range(2, max_clusters + 1),
    )

    candidates = filter_cluster_candidates(candidate_scores, threshold=0.5)
    if not candidates:
        candidates = filter_cluster_candidates(candidate_scores, threshold=0.1)
    if not candidates:
        return 1

    return max(candidates.items(), key=lambda item: (item[1], -item[0]))[0]


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
