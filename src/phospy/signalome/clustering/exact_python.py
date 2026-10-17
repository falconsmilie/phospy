"""Internal exact-Python clustering composition."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.signalome.clustering.backends.exact_python import (
    ExactPythonTreeEngine,
)
from phospy.signalome.clustering.diagnostic_schemas import (
    build_tree_engine_diagnostics,
)
from phospy.signalome.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION,
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.signalome.clustering.tree_engine_adapter import (
    run_clustering_with_tree_engine,
)


@dataclass(frozen=True, slots=True)
class ExactPythonClusteringBackend:
    """Top-level exact-Python clustering backend."""

    name: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    version: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION

    def run(
        self,
        request: SignalomeClusteringEngineRequest,
    ) -> SignalomeClusteringEngineResult:
        return run_clustering_with_tree_engine(
            request=request,
            tree_engine=ExactPythonTreeEngine(),
            clustering_engine=self.name,
            backend_version=self.version,
            backend_diagnostics=build_tree_engine_diagnostics(
                uses_scipy=False,
                linkage_method="ward",
                distance_metric="euclidean",
            ),
        )
