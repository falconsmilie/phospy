"""Pure-Python Ward tree engine for signalome clustering."""

from __future__ import annotations

import heapq
from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from phospy.science.signalomes.clustering.validation import (
    validate_cluster_count_for_site_count,
)


@dataclass(frozen=True, slots=True)
class ExactWardClusterTree:
    n_sites: int
    merges: tuple[tuple[int, int], ...]


def build_cluster_tree(scoring_values: np.ndarray) -> ExactWardClusterTree:
    values = np.asarray(scoring_values, dtype=float)
    if values.ndim != 2:
        raise ValueError("scoring_values must be a 2D array")
    n_sites = int(values.shape[0])
    if n_sites <= 1:
        return ExactWardClusterTree(n_sites=n_sites, merges=())

    centroids: dict[int, np.ndarray] = {
        index: values[index].astype(float, copy=False) for index in range(n_sites)
    }
    sizes: dict[int, int] = {index: 1 for index in range(n_sites)}
    versions: dict[int, int] = {index: 0 for index in range(n_sites)}
    active: set[int] = set(range(n_sites))

    heap: list[tuple[float, int, int, int, int]] = []
    for left in range(n_sites):
        for right in range(left + 1, n_sites):
            heapq.heappush(
                heap,
                (
                    _ward_distance(
                        centroids[left],
                        sizes[left],
                        centroids[right],
                        sizes[right],
                    ),
                    left,
                    right,
                    versions[left],
                    versions[right],
                ),
            )

    merges: list[tuple[int, int]] = []
    next_cluster_id = n_sites
    while len(active) > 1:
        left, right = _pop_next_valid_merge(
            heap=heap,
            active=active,
            versions=versions,
        )
        merges.append((left, right))

        left_size = sizes.pop(left)
        right_size = sizes.pop(right)
        left_centroid = centroids.pop(left)
        right_centroid = centroids.pop(right)
        active.remove(left)
        active.remove(right)

        merged_size = left_size + right_size
        merged_centroid = (
            float(left_size) * left_centroid + float(right_size) * right_centroid
        ) / float(merged_size)
        merged_id = next_cluster_id
        next_cluster_id += 1

        sizes[merged_id] = merged_size
        centroids[merged_id] = merged_centroid
        versions[merged_id] = 0
        active.add(merged_id)

        for other in active:
            if other == merged_id:
                continue
            first = min(other, merged_id)
            second = max(other, merged_id)
            heapq.heappush(
                heap,
                (
                    _ward_distance(
                        centroids[other],
                        sizes[other],
                        centroids[merged_id],
                        sizes[merged_id],
                    ),
                    first,
                    second,
                    versions[first],
                    versions[second],
                ),
            )

    return ExactWardClusterTree(n_sites=n_sites, merges=tuple(merges))


def build_cluster_labels_from_tree(
    *,
    cluster_tree: ExactWardClusterTree,
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
        base = np.zeros(1, dtype=int)
        for count in unique_counts:
            labels_by_count[count] = base.copy()
        return {count: labels_by_count[int(count)] for count in requested_counts}

    current_members: dict[int, np.ndarray] = {
        site_id: np.asarray([site_id], dtype=int) for site_id in range(n_sites)
    }
    current_cluster_count = n_sites
    next_cluster_id = n_sites

    if current_cluster_count in unique_counts:
        labels_by_count[current_cluster_count] = _labels_from_members(
            members=current_members,
            n_sites=n_sites,
        )

    for left, right in cluster_tree.merges:
        left_members = current_members.pop(left)
        right_members = current_members.pop(right)
        merged = np.concatenate([left_members, right_members])
        current_members[next_cluster_id] = merged
        next_cluster_id += 1
        current_cluster_count -= 1
        if current_cluster_count in unique_counts:
            labels_by_count[current_cluster_count] = _labels_from_members(
                members=current_members,
                n_sites=n_sites,
            )
        if current_cluster_count <= min(unique_counts):
            break

    return {count: labels_by_count[int(count)].copy() for count in requested_counts}


@dataclass(frozen=True, slots=True)
class ExactPythonTreeEngine:
    """Pure-Python Ward tree engine."""

    name: str = "exact_python_tree"
    version: str = "1"

    def build_tree(self, values: np.ndarray) -> ExactWardClusterTree:
        return build_cluster_tree(values)

    def labels_for_counts(
        self,
        *,
        tree: object,
        counts: Iterable[int],
    ) -> dict[int, np.ndarray]:
        if not isinstance(tree, ExactWardClusterTree):
            raise TypeError(
                "exact tree engine expected an ExactWardClusterTree instance"
            )
        return build_cluster_labels_from_tree(cluster_tree=tree, cluster_counts=counts)


def _ward_distance(
    left_centroid: np.ndarray,
    left_size: int,
    right_centroid: np.ndarray,
    right_size: int,
) -> float:
    delta = np.asarray(left_centroid, dtype=float) - np.asarray(
        right_centroid, dtype=float
    )
    squared_norm = float(np.dot(delta, delta))
    return (
        (float(left_size) * float(right_size))
        / float(left_size + right_size)
        * squared_norm
    )


def _pop_next_valid_merge(
    *,
    heap: list[tuple[float, int, int, int, int]],
    active: set[int],
    versions: dict[int, int],
) -> tuple[int, int]:
    while heap:
        _, left, right, left_version, right_version = heapq.heappop(heap)
        if left not in active or right not in active:
            continue
        if versions[left] != left_version or versions[right] != right_version:
            continue
        return left, right
    raise RuntimeError("failed to resolve a valid merge from ward clustering heap")


def _labels_from_members(
    *,
    members: dict[int, np.ndarray],
    n_sites: int,
) -> np.ndarray:
    labels = np.zeros(n_sites, dtype=int)
    sorted_clusters = sorted(members.items(), key=lambda item: int(item[0]))
    for label, (_, cluster_members) in enumerate(sorted_clusters):
        labels[np.asarray(cluster_members, dtype=int)] = int(label)
    return labels


__all__ = [
    "ExactPythonTreeEngine",
    "ExactWardClusterTree",
    "build_cluster_labels_from_tree",
    "build_cluster_tree",
]
