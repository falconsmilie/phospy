from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.signalomes.clustering as clustering
from tests.support.performance_contracts import (
    SIGNALOME_CLUSTER_TREE_BENCHMARK_N_KINASES,
    SIGNALOME_CLUSTER_TREE_BENCHMARK_N_SITES,
    SIGNALOME_CLUSTER_TREE_PEAK_MIB_MAX,
    SIGNALOME_CLUSTER_TREE_RUNTIME_SECONDS_MAX,
    deterministic_matrix,
    measure_runtime_and_peak_mib,
)

pytestmark = pytest.mark.performance


def _build_signalome_scoring_matrix(*, n_sites: int, n_kinases: int) -> pd.DataFrame:
    base_matrix = deterministic_matrix(
        n_sites=n_sites,
        n_samples=n_kinases,
        seed=42042,
    )
    values = base_matrix.to_numpy(dtype=float, copy=False)
    rng = np.random.default_rng(42042)

    n_module_profiles = 10
    module_profiles = rng.normal(
        loc=10.0,
        scale=1.6,
        size=(n_module_profiles, int(n_kinases)),
    )
    module_assignments = (np.arange(int(n_sites), dtype=int) * 7) % int(
        n_module_profiles
    )
    realistic_values = module_profiles[module_assignments] + rng.normal(
        loc=0.0,
        scale=0.35,
        size=values.shape,
    )
    realistic_values = np.round(realistic_values, decimals=6)
    matrix = pd.DataFrame(
        realistic_values,
        index=base_matrix.index.copy(),
        columns=base_matrix.columns.copy(),
        dtype=float,
    )
    matrix.columns = pd.Index(
        [f"KINASE_{index + 1:03d}" for index in range(int(matrix.shape[1]))],
        name="kinase",
    )
    return matrix


def test_cluster_tree_builder_medium_input_benchmark() -> None:
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

    # Keep this threshold loose for CI stability while still guarding against
    # accidental severe runtime/memory regressions around exact-tree construction.
    assert elapsed_seconds < SIGNALOME_CLUSTER_TREE_RUNTIME_SECONDS_MAX
    assert peak_mib < SIGNALOME_CLUSTER_TREE_PEAK_MIB_MAX
    assert cluster_tree.n_sites == n_sites
    assert len(cluster_tree.merges) == n_sites - 1
    assert merges.shape == (n_sites - 1, 2)
    assert int(merge_children.min()) >= 0
    assert int(merge_children.max()) < (2 * n_sites - 1)
    assert np.unique(merge_children).size == merge_children.size
