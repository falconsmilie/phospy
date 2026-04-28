"""Signalome clustering public facade and internal backend boundary."""

from __future__ import annotations

from typing import Any

from phospy.signalomes.clustering import exact_python as _exact
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION,
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL_VERSION,
    SignalomeClusteringBackendRequest,
    SignalomeClusteringBackendResult,
)
from phospy.signalomes.clustering.protocol import SignalomeClusteringBackend
from phospy.signalomes.clustering.selection import (
    available_clustering_backends,
    resolve_clustering_backend,
    run_clustering_backend,
)

for _name in dir(_exact):
    if _name.startswith("__"):
        continue
    globals()[_name] = getattr(_exact, _name)

# Explicit bindings keep static analysers aware of runtime-exported constants.
SIGNALOME_CLUSTER_TREE_BACKEND_EXACT = _exact.SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
MAX_FULL_CORRELATION_SITE_COUNT = _exact.MAX_FULL_CORRELATION_SITE_COUNT


_PATCHABLE_EXACT_SYMBOLS = (
    "_build_cluster_tree",
    "build_cluster_labels_from_tree",
    "_resolve_pre_scoring_module_selection",
    "build_correlation_matrix_with_exclusions",
)


def _sync_exact_monkeypatch_hooks() -> None:
    for symbol in _PATCHABLE_EXACT_SYMBOLS:
        if symbol in globals():
            setattr(_exact, symbol, globals()[symbol])


def run_signalome_clustering_backend(
    *,
    scoring_matrix: Any,
    site_to_protein: Any,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    cluster_tree_backend: str = SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
    candidate_scoring_backend: str | None = None,
    max_exact_cluster_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_correlation_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    backend_name: str = SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
) -> SignalomeClusteringBackendResult:
    """Run the internal backend protocol and return typed backend output."""

    _sync_exact_monkeypatch_hooks()
    return run_clustering_backend(
        request=SignalomeClusteringBackendRequest(
            scoring_matrix=scoring_matrix,
            site_to_protein=site_to_protein,
            requested_module_count=requested_module_count,
            primary_threshold=primary_threshold,
            fallback_threshold=fallback_threshold,
            max_clusters=max_clusters,
            cluster_tree_backend=cluster_tree_backend,
            candidate_scoring_backend=candidate_scoring_backend,
            max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
            max_full_correlation_sites=max_full_correlation_sites,
        ),
        backend_name=backend_name,
    )


def cluster_sites(*args: Any, **kwargs: Any) -> Any:
    _sync_exact_monkeypatch_hooks()
    return _exact.cluster_sites(*args, **kwargs)


def cluster_sites_with_diagnostics(*args: Any, **kwargs: Any) -> Any:
    _sync_exact_monkeypatch_hooks()
    return _exact.cluster_sites_with_diagnostics(*args, **kwargs)


def select_module_count(*args: Any, **kwargs: Any) -> Any:
    _sync_exact_monkeypatch_hooks()
    return _exact.select_module_count(*args, **kwargs)


def select_module_count_with_diagnostics(*args: Any, **kwargs: Any) -> Any:
    _sync_exact_monkeypatch_hooks()
    return _exact.select_module_count_with_diagnostics(*args, **kwargs)


def fit_cluster_labels(*args: Any, **kwargs: Any) -> Any:
    _sync_exact_monkeypatch_hooks()
    return _exact.fit_cluster_labels(*args, **kwargs)


__all__ = list(_exact.__all__) + [
    "SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON",
    "SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION",
    "SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL",
    "SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL_VERSION",
    "SignalomeClusteringBackend",
    "SignalomeClusteringBackendRequest",
    "SignalomeClusteringBackendResult",
    "available_clustering_backends",
    "resolve_clustering_backend",
    "run_clustering_backend",
    "run_signalome_clustering_backend",
]
