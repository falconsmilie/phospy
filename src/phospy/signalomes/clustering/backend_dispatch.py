"""Backend registry/dispatch for signalome clustering."""

from __future__ import annotations

from phospy.signalomes.clustering.exact_python import ExactPythonClusteringBackend
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.signalomes.clustering.protocol import SignalomeClusteringEngine
from phospy.signalomes.clustering.scipy_hierarchical import (
    ScipyHierarchicalClusteringBackend,
)

_BACKENDS: dict[str, SignalomeClusteringEngine] = {
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON: ExactPythonClusteringBackend(),
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL: (
        ScipyHierarchicalClusteringBackend()
    ),
}


def available_clustering_engines() -> tuple[str, ...]:
    """Return supported internal signalome clustering backend names."""

    return tuple(sorted(_BACKENDS))


def resolve_clustering_engine(name: str) -> SignalomeClusteringEngine:
    """Resolve a backend by explicit name."""

    if name in _BACKENDS:
        return _BACKENDS[name]
    available = ", ".join(available_clustering_engines())
    raise ValueError(
        "unsupported signalome clustering backend "
        f"{name!r}; expected one of: {available}"
    )


def run_clustering_engine(
    *,
    request: SignalomeClusteringEngineRequest,
    clustering_engine: str = SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
) -> SignalomeClusteringEngineResult:
    """Run a resolved backend for the provided request."""

    return resolve_clustering_engine(clustering_engine).run(request)


__all__ = [
    "available_clustering_engines",
    "resolve_clustering_engine",
    "run_clustering_engine",
]
