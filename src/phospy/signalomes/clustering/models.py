"""Internal signalome clustering engine models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.signalomes.clustering.diagnostic_schemas import (
    SignalomeBackendDiagnostics,
    SignalomeCandidateScoringSamplingDiagnostics,
    SignalomeClusteringLimitMetadata,
    SignalomeClusteringThresholdMetadata,
    validate_backend_diagnostics,
    validate_candidate_scoring_sampling_diagnostics,
    validate_limit_metadata,
    validate_threshold_metadata,
)
from phospy.signalomes.clustering.policies import SignalomeCandidateScoringPolicy
from phospy.signalomes.models import SignalomeModuleSelectionDiagnostics

SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON = "exact_python"
SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION = "1"
SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL = "scipy_hierarchical"
SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION = "1"


@dataclass(frozen=True, slots=True)
class SignalomeClusteringEngineRequest:
    """Execution request for a clustering engine."""

    scoring_matrix: pd.DataFrame
    site_to_protein: pd.Series
    requested_module_count: int | None
    primary_threshold: float
    fallback_threshold: float
    max_clusters: int
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None
    max_exact_tree_sites: int | None
    max_full_candidate_scoring_sites: int


@dataclass(frozen=True, slots=True)
class SignalomeClusteringEngineResult:
    """Typed clustering engine output with diagnostics/provenance metadata."""

    site_clusters: pd.Series
    protein_modules: pd.Series
    selected_module_count: int
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics
    backend_name: str
    backend_version: str
    approximation_used: bool
    exact_cluster_tree_built: bool
    tree_implementation: str
    candidate_scoring_mode: str
    candidate_scoring_evaluated: bool
    candidate_scoring_skip_reason: str | None
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None
    backend_diagnostics: SignalomeBackendDiagnostics | None
    threshold_metadata: SignalomeClusteringThresholdMetadata
    limit_metadata: SignalomeClusteringLimitMetadata

    def __post_init__(self) -> None:
        if self.candidate_scoring_sampling is not None:
            validate_candidate_scoring_sampling_diagnostics(
                self.candidate_scoring_sampling,
                field_name="clustering_engine_result.candidate_scoring_sampling",
            )
        if self.backend_diagnostics is not None:
            validate_backend_diagnostics(
                self.backend_diagnostics,
                field_name="clustering_engine_result.backend_diagnostics",
            )
        validate_threshold_metadata(
            self.threshold_metadata,
            field_name="clustering_engine_result.threshold_metadata",
        )
        validate_limit_metadata(
            self.limit_metadata,
            field_name="clustering_engine_result.limit_metadata",
        )


__all__ = [
    "SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON",
    "SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION",
    "SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL",
    "SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION",
    "SignalomeClusteringEngineRequest",
    "SignalomeClusteringEngineResult",
]
