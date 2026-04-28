"""Tree-engine implementations for signalome clustering orchestration."""

from __future__ import annotations

from phospy.signalomes.clustering.backends.exact_python import (
    ExactPythonTreeEngine,
    ExactWardClusterTree,
)
from phospy.signalomes.clustering.backends.exact_python import (
    build_cluster_labels_from_tree as build_exact_cluster_labels_from_tree,
)
from phospy.signalomes.clustering.backends.exact_python import (
    build_cluster_tree as build_exact_cluster_tree,
)
from phospy.signalomes.clustering.backends.scipy_hierarchical import (
    ScipyHierarchicalTreeEngine,
    ScipyWardClusterTree,
)
from phospy.signalomes.clustering.backends.scipy_hierarchical import (
    build_cluster_labels_from_tree as build_scipy_cluster_labels_from_tree,
)
from phospy.signalomes.clustering.backends.scipy_hierarchical import (
    build_cluster_tree as build_scipy_cluster_tree,
)
from phospy.signalomes.clustering.backends.scipy_hierarchical import (
    engine_diagnostics as scipy_engine_diagnostics,
)

__all__ = [
    "ExactPythonTreeEngine",
    "ExactWardClusterTree",
    "ScipyHierarchicalTreeEngine",
    "ScipyWardClusterTree",
    "build_exact_cluster_labels_from_tree",
    "build_exact_cluster_tree",
    "build_scipy_cluster_labels_from_tree",
    "build_scipy_cluster_tree",
    "scipy_engine_diagnostics",
]
