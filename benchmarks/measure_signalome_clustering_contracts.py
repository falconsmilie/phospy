#!/usr/bin/env python3
"""Benchmark signalome module-selection correlation-path contracts.

Targets:
- `phospy.signalomes.clustering.select_module_count_with_diagnostics`
- `phospy.signalomes.clustering.MAX_FULL_CORRELATION_SITE_COUNT`
- `phospy.signalomes.clustering.MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER`
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

APPROXIMATION_REASON_TOKEN = "Used sampled within-cluster correlation estimates"


def _build_scoring_matrix(*, n_sites: int, n_kinases: int, seed: int) -> pd.DataFrame:
    from tests.support.performance_contracts import deterministic_matrix

    matrix = deterministic_matrix(n_sites=n_sites, n_samples=n_kinases, seed=seed)
    matrix.columns = pd.Index(
        [f"KINASE_{index + 1:03d}" for index in range(matrix.shape[1])],
        name="kinase",
    )
    return matrix


def _patch_cluster_tree_for_contract_measurement(
    clustering: object,
) -> tuple[object, object]:
    original_tree_builder = clustering._build_cluster_tree
    original_label_builder = clustering.build_cluster_labels_from_tree

    def _stub_build_cluster_tree(scoring_values: np.ndarray) -> object:
        n_sites = int(np.asarray(scoring_values, dtype=float).shape[0])
        return clustering._WardClusterTree(n_sites=n_sites, merges=())

    def _stub_build_cluster_labels_from_tree(
        *,
        cluster_tree: object,
        cluster_counts: object,
    ) -> dict[int, np.ndarray]:
        n_sites = int(cluster_tree.n_sites)
        labels_by_count: dict[int, np.ndarray] = {}
        for requested_count in [int(value) for value in cluster_counts]:
            resolved_count = max(1, min(int(requested_count), n_sites))
            if resolved_count == 1:
                labels = np.zeros(n_sites, dtype=int)
            else:
                labels = np.arange(n_sites, dtype=int) % resolved_count
            labels_by_count[requested_count] = labels.astype(int, copy=False)
        return labels_by_count

    clustering._build_cluster_tree = _stub_build_cluster_tree
    clustering.build_cluster_labels_from_tree = _stub_build_cluster_labels_from_tree
    return original_tree_builder, original_label_builder


def _restore_cluster_tree_builders(
    clustering: object,
    original_tree_builder: object,
    original_label_builder: object,
) -> None:
    clustering._build_cluster_tree = original_tree_builder
    clustering.build_cluster_labels_from_tree = original_label_builder


def main() -> None:
    import phospy.signalomes.clustering as clustering
    from tests.support.performance_contracts import measure_runtime_and_peak_mib

    original_tree_builder, original_label_builder = (
        _patch_cluster_tree_for_contract_measurement(clustering)
    )
    try:
        below_threshold_sites = 600
        above_threshold_sites = clustering.MAX_FULL_CORRELATION_SITE_COUNT + 50

        below_threshold_matrix = _build_scoring_matrix(
            n_sites=below_threshold_sites,
            n_kinases=40,
            seed=401,
        )
        above_threshold_matrix = _build_scoring_matrix(
            n_sites=above_threshold_sites,
            n_kinases=12,
            seed=402,
        )

        below_diagnostics, below_runtime_seconds, below_peak_mib = (
            measure_runtime_and_peak_mib(
                lambda: clustering.select_module_count_with_diagnostics(
                    scoring_values=below_threshold_matrix,
                    max_clusters=3,
                ),
                warmup=True,
            )
        )
        above_diagnostics, above_runtime_seconds, above_peak_mib = (
            measure_runtime_and_peak_mib(
                lambda: clustering.select_module_count_with_diagnostics(
                    scoring_values=above_threshold_matrix,
                    max_clusters=3,
                ),
                warmup=True,
            )
        )
    finally:
        _restore_cluster_tree_builders(
            clustering,
            original_tree_builder,
            original_label_builder,
        )

    full_matrix_mib = (above_threshold_sites * above_threshold_sites * 8) / (
        1024 * 1024
    )
    approximation_reported = APPROXIMATION_REASON_TOKEN in above_diagnostics.reason
    below_uses_approximation = APPROXIMATION_REASON_TOKEN in below_diagnostics.reason

    print("cluster_tree_measurement_mode=stubbed_module_selection_contract_scoring")
    print(
        f"max_full_correlation_site_count={clustering.MAX_FULL_CORRELATION_SITE_COUNT}"
    )
    print(
        "max_approx_correlation_samples_per_cluster="
        f"{clustering.MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER}"
    )
    print(f"below_threshold_sites={below_threshold_sites}")
    print(f"below_threshold_runtime_seconds={below_runtime_seconds:.6f}")
    print(f"below_threshold_peak_mib={below_peak_mib:.3f}")
    print(
        f"below_threshold_reports_approximation={str(below_uses_approximation).lower()}"
    )
    print(f"above_threshold_sites={above_threshold_sites}")
    print(f"above_threshold_runtime_seconds={above_runtime_seconds:.6f}")
    print(f"above_threshold_peak_mib={above_peak_mib:.3f}")
    print(
        f"above_threshold_reports_approximation={str(approximation_reported).lower()}"
    )
    print(f"above_threshold_full_matrix_theoretical_mib={full_matrix_mib:.3f}")


if __name__ == "__main__":
    main()
