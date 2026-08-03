"""Internal signalome clustering engine models."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.science.signalomes.clustering.diagnostic_schemas import (
    SignalomeBackendDiagnostics,
    SignalomeCandidateScoringSamplingDiagnostics,
    SignalomeClusteringLimitMetadata,
    SignalomeClusteringThresholdMetadata,
    validate_backend_diagnostics,
    validate_candidate_scoring_sampling_diagnostics,
    validate_limit_metadata,
    validate_threshold_metadata,
)
from phospy.science.signalomes.clustering.policies import (
    SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES,
    SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS,
    SignalomeCandidateScoringPolicy,
)
from phospy.science.signalomes.models import (
    SignalomeClusteringPreparationDiagnostics,
    SignalomeModuleSelectionDiagnostics,
    default_signalome_clustering_preparation_diagnostics,
)

SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON = "exact_python"
SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION = "1"
SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL = "scipy_hierarchical"
SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION = "1"


@dataclass(frozen=True, slots=True, init=False)
class SignalomeClusteringEngineRequest:
    """Execution request for a clustering engine."""

    scoring_matrix: pd.DataFrame
    site_to_protein_group_id: pd.Series
    requested_module_count: int | None
    primary_threshold: float
    fallback_threshold: float
    max_clusters: int
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None
    max_exact_tree_sites: int | None
    max_full_candidate_scoring_sites: int
    module_selection_stability_perturbations: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS
    )
    module_selection_stability_seed: int | None = None
    module_selection_stability_max_sites: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES
    )

    def __init__(
        self,
        *,
        scoring_matrix: pd.DataFrame,
        site_to_protein_group_id: pd.Series | None = None,
        site_to_protein: pd.Series | None = None,
        requested_module_count: int | None,
        primary_threshold: float,
        fallback_threshold: float,
        max_clusters: int,
        candidate_scoring_policy: SignalomeCandidateScoringPolicy | None,
        max_exact_tree_sites: int | None,
        max_full_candidate_scoring_sites: int,
        module_selection_stability_perturbations: int = (
            SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS
        ),
        module_selection_stability_seed: int | None = None,
        module_selection_stability_max_sites: int = (
            SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES
        ),
    ) -> None:
        """Create a clustering request, accepting legacy site_to_protein alias."""

        if site_to_protein_group_id is None:
            if site_to_protein is None:
                raise TypeError(
                    "SignalomeClusteringEngineRequest requires site_to_protein_group_id"
                )
            site_to_protein_group_id = site_to_protein
        elif site_to_protein is not None and not site_to_protein.equals(
            site_to_protein_group_id
        ):
            raise ValueError(
                "SignalomeClusteringEngineRequest received conflicting "
                "site_to_protein_group_id and legacy site_to_protein mappings"
            )
        object.__setattr__(self, "scoring_matrix", scoring_matrix)
        object.__setattr__(self, "site_to_protein_group_id", site_to_protein_group_id)
        object.__setattr__(self, "requested_module_count", requested_module_count)
        object.__setattr__(self, "primary_threshold", primary_threshold)
        object.__setattr__(self, "fallback_threshold", fallback_threshold)
        object.__setattr__(self, "max_clusters", max_clusters)
        object.__setattr__(self, "candidate_scoring_policy", candidate_scoring_policy)
        object.__setattr__(self, "max_exact_tree_sites", max_exact_tree_sites)
        object.__setattr__(
            self,
            "max_full_candidate_scoring_sites",
            max_full_candidate_scoring_sites,
        )
        object.__setattr__(
            self,
            "module_selection_stability_perturbations",
            module_selection_stability_perturbations,
        )
        object.__setattr__(
            self, "module_selection_stability_seed", module_selection_stability_seed
        )
        object.__setattr__(
            self,
            "module_selection_stability_max_sites",
            module_selection_stability_max_sites,
        )

    @property
    def site_to_protein(self) -> pd.Series:
        """Legacy alias for Signalome site-to-protein-group-id mapping."""

        return self.site_to_protein_group_id


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
    clustering_preparation_diagnostics: SignalomeClusteringPreparationDiagnostics = (
        field(default_factory=default_signalome_clustering_preparation_diagnostics)
    )

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
