"""Top-level SciPy hierarchical clustering backend implementation."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.signalomes.clustering.backends.scipy_hierarchical import (
    ScipyHierarchicalTreeEngine,
)
from phospy.signalomes.clustering.backends.scipy_hierarchical import (
    engine_diagnostics as scipy_tree_engine_diagnostics,
)
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION,
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.signalomes.clustering.tree_engine_adapter import (
    run_clustering_with_tree_engine,
)


@dataclass(frozen=True, slots=True)
class ScipyHierarchicalClusteringBackend:
    """SciPy-backed hierarchical clustering backend."""

    name: str = SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    version: str = SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION

    def run(
        self,
        request: SignalomeClusteringEngineRequest,
    ) -> SignalomeClusteringEngineResult:
        return run_clustering_with_tree_engine(
            request=request,
            tree_engine=ScipyHierarchicalTreeEngine(),
            clustering_engine=self.name,
            backend_version=self.version,
            backend_diagnostics=scipy_tree_engine_diagnostics(),
        )


__all__ = ["ScipyHierarchicalClusteringBackend"]
