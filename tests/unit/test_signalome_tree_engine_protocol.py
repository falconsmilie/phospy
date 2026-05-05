from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.signalomes.clustering.backends.exact_python import ExactPythonTreeEngine
from phospy.signalomes.clustering.backends.scipy_hierarchical import (
    ScipyHierarchicalTreeEngine,
)
from phospy.signalomes.clustering.contracts import ClusterTreeEngine
from phospy.signalomes.clustering.models import (
    SignalomeClusteringEngineRequest,
)
from phospy.signalomes.clustering.tree_engine_adapter import (
    run_clustering_with_tree_engine,
)


def _scoring_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [
            [1.0, 0.1, 0.2],
            [0.9, 0.2, 0.1],
            [0.0, 1.0, 0.9],
            [0.1, 0.9, 1.0],
        ],
        index=["P1;S1;", "P2;S2;", "P3;S3;", "P4;S4;"],
        columns=["K1", "K2", "K3"],
        dtype=float,
    )


def _site_to_protein() -> pd.Series:
    return pd.Series(
        ["P1", "P2", "P3", "P4"],
        index=pd.Index(["P1;S1;", "P2;S2;", "P3;S3;", "P4;S4;"], name="site_id"),
        name="protein_id",
        dtype=str,
    )


def _same_partition(left: np.ndarray, right: np.ndarray) -> bool:
    left_equal = left[:, None] == left[None, :]
    right_equal = right[:, None] == right[None, :]
    return bool(np.array_equal(left_equal, right_equal))


def test_exact_and_scipy_tree_engines_satisfy_contract_runtime_shape() -> None:
    exact_engine: ClusterTreeEngine = ExactPythonTreeEngine()
    scipy_engine: ClusterTreeEngine = ScipyHierarchicalTreeEngine()
    values = _scoring_matrix().to_numpy(dtype=float)

    exact_tree = exact_engine.build_tree(values)
    scipy_tree = scipy_engine.build_tree(values)
    exact_labels = exact_engine.labels_for_counts(tree=exact_tree, counts=[2, 3])
    scipy_labels = scipy_engine.labels_for_counts(tree=scipy_tree, counts=[2, 3])

    assert set(exact_labels) == {2, 3}
    assert set(scipy_labels) == {2, 3}
    assert exact_labels[2].shape == (4,)
    assert scipy_labels[2].shape == (4,)


def test_exact_and_scipy_tree_engines_produce_equivalent_partitions() -> None:
    exact_engine = ExactPythonTreeEngine()
    scipy_engine = ScipyHierarchicalTreeEngine()
    values = _scoring_matrix().to_numpy(dtype=float)

    exact_tree = exact_engine.build_tree(values)
    scipy_tree = scipy_engine.build_tree(values)
    exact_labels = exact_engine.labels_for_counts(tree=exact_tree, counts=[2])[2]
    scipy_labels = scipy_engine.labels_for_counts(tree=scipy_tree, counts=[2])[2]

    assert _same_partition(exact_labels, scipy_labels)


def test_shared_orchestration_runs_with_either_tree_engine() -> None:
    request = SignalomeClusteringEngineRequest(
        scoring_matrix=_scoring_matrix(),
        site_to_protein=_site_to_protein(),
        requested_module_count=None,
        primary_threshold=0.5,
        fallback_threshold=0.1,
        max_clusters=4,
        candidate_scoring_policy="full",
        max_exact_tree_sites=10,
        max_full_candidate_scoring_sites=10,
    )

    exact_result = run_clustering_with_tree_engine(
        request=request,
        tree_engine=ExactPythonTreeEngine(),
        clustering_engine="exact_python",
        backend_version="1",
        backend_diagnostics={
            "uses_scipy": False,
            "linkage_method": "ward",
            "distance_metric": "euclidean",
        },
    )
    scipy_result = run_clustering_with_tree_engine(
        request=request,
        tree_engine=ScipyHierarchicalTreeEngine(),
        clustering_engine="scipy_hierarchical",
        backend_version="1",
        backend_diagnostics={
            "uses_scipy": True,
            "linkage_method": "ward",
            "distance_metric": "euclidean",
        },
    )

    assert _same_partition(
        exact_result.site_clusters.to_numpy(dtype=int, copy=False),
        scipy_result.site_clusters.to_numpy(dtype=int, copy=False),
    )
    assert _same_partition(
        exact_result.protein_modules.to_numpy(dtype=int, copy=False),
        scipy_result.protein_modules.to_numpy(dtype=int, copy=False),
    )
    assert exact_result.tree_implementation == "exact_python_tree"
    assert scipy_result.tree_implementation == "scipy_hierarchical_tree"
    assert exact_result.backend_diagnostics is not None
    assert scipy_result.backend_diagnostics is not None
    assert (
        exact_result.backend_diagnostics["tree_implementation"]
        == exact_result.tree_implementation
    )
    assert (
        scipy_result.backend_diagnostics["tree_implementation"]
        == scipy_result.tree_implementation
    )
