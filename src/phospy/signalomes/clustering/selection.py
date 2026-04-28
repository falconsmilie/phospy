"""Backend selection/dispatch for signalome clustering."""

from __future__ import annotations

from phospy.signalomes.clustering.exact_python import ExactPythonClusteringBackend
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    SignalomeClusteringBackendRequest,
    SignalomeClusteringBackendResult,
)
from phospy.signalomes.clustering.protocol import SignalomeClusteringBackend
from phospy.signalomes.clustering.scipy_hierarchical import (
    ScipyHierarchicalClusteringBackend,
)

_BACKENDS: dict[str, SignalomeClusteringBackend] = {
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON: ExactPythonClusteringBackend(),
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL: (
        ScipyHierarchicalClusteringBackend()
    ),
}


def available_clustering_backends() -> tuple[str, ...]:
    """Return supported internal signalome clustering backend names."""

    return tuple(sorted(_BACKENDS))


def resolve_clustering_backend(name: str) -> SignalomeClusteringBackend:
    """Resolve a backend by explicit name."""

    if name in _BACKENDS:
        return _BACKENDS[name]
    available = ", ".join(available_clustering_backends())
    raise ValueError(
        "unsupported signalome clustering backend "
        f"{name!r}; expected one of: {available}"
    )


def run_clustering_backend(
    *,
    request: SignalomeClusteringBackendRequest,
    backend_name: str = SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
) -> SignalomeClusteringBackendResult:
    """Run a resolved backend for the provided request."""

    return resolve_clustering_backend(backend_name).run(request)


__all__ = [
    "available_clustering_backends",
    "resolve_clustering_backend",
    "run_clustering_backend",
]
