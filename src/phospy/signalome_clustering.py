from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

__all__ = ["cluster_sites"]


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

    module_count = (
        requested_module_count
        if requested_module_count is not None
        else select_module_count(scoring_matrix)
    )
    module_count = max(1, min(module_count, n_sites))

    if module_count == 1:
        labels = np.ones(n_sites, dtype=int)
    else:
        labels = (
            AgglomerativeClustering(
                n_clusters=module_count,
                linkage="ward",
            )
            .fit_predict(scoring_matrix.to_numpy(dtype=float))
            .astype(int)
            + 1
        )

    return pd.Series(labels, index=scoring_matrix.index, dtype=int, name="site_cluster")


def select_module_count(scoring_matrix: pd.DataFrame) -> int:
    """Choose a signalome module count from the scoring matrix."""

    n_sites = scoring_matrix.shape[0]
    if n_sites <= 1:
        return 1

    max_clusters = min(10, n_sites)
    if max_clusters < 2:
        return 1

    site_correlations = np.corrcoef(scoring_matrix.to_numpy(dtype=float))
    candidates = score_cluster_candidates(
        scoring_matrix=scoring_matrix,
        site_correlations=site_correlations,
        threshold=0.5,
        cluster_range=range(2, max_clusters + 1),
    )
    if not candidates:
        candidates = score_cluster_candidates(
            scoring_matrix=scoring_matrix,
            site_correlations=site_correlations,
            threshold=0.1,
            cluster_range=range(2, max_clusters + 1),
        )
    if not candidates:
        return 1

    return max(candidates.items(), key=lambda item: (item[1], -item[0]))[0]


def score_cluster_candidates(
    *,
    scoring_matrix: pd.DataFrame,
    site_correlations: np.ndarray,
    threshold: float,
    cluster_range: Iterable[int],
) -> dict[int, float]:
    """Score candidate module counts using within-cluster correlations."""

    candidates: dict[int, float] = {}
    for cluster_count in cluster_range:
        labels = (
            AgglomerativeClustering(
                n_clusters=cluster_count,
                linkage="ward",
            )
            .fit_predict(scoring_matrix.to_numpy(dtype=float))
            .astype(int)
        )
        cluster_medians = [
            cluster_median_correlation(site_correlations, labels, label)
            for label in np.unique(labels)
        ]
        if cluster_medians and all(median >= threshold for median in cluster_medians):
            candidates[cluster_count] = float(np.mean(cluster_medians))
    return candidates


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
