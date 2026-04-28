"""Internal signalome clustering backend models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.signalomes.models import SignalomeModuleSelectionDiagnostics

SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON = "exact_python"
SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION = "1"


@dataclass(frozen=True, slots=True)
class SignalomeClusteringBackendRequest:
    """Execution request for a clustering backend."""

    scoring_matrix: pd.DataFrame
    site_to_protein: pd.Series
    requested_module_count: int | None
    primary_threshold: float
    fallback_threshold: float
    max_clusters: int
    cluster_tree_backend: str
    candidate_scoring_backend: str | None
    max_exact_cluster_tree_sites: int | None
    max_full_correlation_sites: int


@dataclass(frozen=True, slots=True)
class SignalomeClusteringBackendResult:
    """Typed clustering backend output with diagnostics/provenance metadata."""

    site_clusters: pd.Series
    protein_modules: pd.Series
    selected_module_count: int
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics
    backend_name: str
    backend_version: str
    approximation_used: bool
    exact_cluster_tree_built: bool
    cluster_tree_backend: str
    candidate_scoring_mode: str
    candidate_scoring_evaluated: bool
    candidate_scoring_skip_reason: str | None
    candidate_scoring_sampling: dict[str, object] | None
    threshold_metadata: dict[str, float | None]
    limit_metadata: dict[str, int | None]


__all__ = [
    "SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON",
    "SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION",
    "SignalomeClusteringBackendRequest",
    "SignalomeClusteringBackendResult",
]
