from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.signalomes.clustering as clustering
from tests.support.performance_contracts import (
    deterministic_matrix,
    measure_runtime_and_peak_mib,
)

pytestmark = pytest.mark.performance

SIGNALOME_CLUSTER_TREE_BENCHMARK_N_SITES = 500
SIGNALOME_CLUSTER_TREE_BENCHMARK_N_KINASES = 8
SIGNALOME_CLUSTER_TREE_RUNTIME_SECONDS_MAX = 15.0
SIGNALOME_CLUSTER_TREE_PEAK_MIB_MAX = 256.0


def _build_signalome_scoring_matrix(*, n_sites: int, n_kinases: int) -> pd.DataFrame:
    matrix = deterministic_matrix(
        n_sites=n_sites,
        n_samples=n_kinases,
        seed=42042,
    )
    matrix.columns = pd.Index(
        [f"KINASE_{index + 1:03d}" for index in range(matrix.shape[1])],
        name="kinase",
    )
    return matrix


def test_signalome_cluster_tree_builder_medium_input_runtime() -> None:
    n_sites = SIGNALOME_CLUSTER_TREE_BENCHMARK_N_SITES
    scoring_matrix = _build_signalome_scoring_matrix(
        n_sites=n_sites,
        n_kinases=SIGNALOME_CLUSTER_TREE_BENCHMARK_N_KINASES,
    )
    scoring_values = scoring_matrix.to_numpy(dtype=float, copy=False)

    cluster_tree, elapsed_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: clustering._build_cluster_tree(scoring_values),
        warmup=True,
    )

    merges = np.asarray(cluster_tree.merges, dtype=int)
    merge_children = merges.reshape(-1)

    assert elapsed_seconds < SIGNALOME_CLUSTER_TREE_RUNTIME_SECONDS_MAX
    assert peak_mib < SIGNALOME_CLUSTER_TREE_PEAK_MIB_MAX
    assert cluster_tree.n_sites == n_sites
    assert len(cluster_tree.merges) == n_sites - 1
    assert merges.shape == (n_sites - 1, 2)
    assert int(merge_children.min()) >= 0
    assert int(merge_children.max()) < (2 * n_sites - 1)
    assert np.unique(merge_children).size == merge_children.size
