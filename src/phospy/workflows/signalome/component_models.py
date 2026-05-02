"""Shared component dataclasses for signalome workflow execution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.prediction.scoring import DownstreamScoreSelectionPolicy
from phospy.signalomes.clustering import (
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    ClusterSitesResult,
)
from phospy.signalomes.clustering.diagnostic_schemas import (
    SignalomeBackendDiagnostics,
    SignalomeCandidateScoringSamplingDiagnostics,
    validate_backend_diagnostics,
    validate_candidate_scoring_sampling_diagnostics,
)
from phospy.signalomes.models import SignalomeNetworkCorrelationDiagnostics


@dataclass(frozen=True, slots=True)
class SignalomeExecutionMetadata:
    prediction_sites: int
    prediction_kinases: int
    downstream_score_sites: int
    downstream_score_kinases: int
    downstream_score_source: str
    downstream_score_selection_policy: DownstreamScoreSelectionPolicy | None = None


@dataclass(frozen=True, slots=True)
class SignalomeScaleGuardDecision:
    site_count: int
    input_protein_count: int
    input_kinase_count: int
    selected_module_count: int
    candidate_module_counts_evaluated: int
    candidate_module_count_upper_bound: int
    clustering_engine: str
    clustering_engine_version: str
    backend_diagnostics: SignalomeBackendDiagnostics | None
    tree_implementation: str
    tree_generation_backend: str
    tree_generation_mode: str
    tree_generation_is_approximate: bool
    tree_generation_scope: str
    tree_generation_guard_triggered: bool
    candidate_scoring_policy: str
    candidate_scoring_requested_policy: str
    candidate_scoring_strategy: str
    candidate_scoring_is_approximate: bool
    candidate_scoring_guard_triggered: bool
    candidate_scoring_sampled_site_total: int | None
    candidate_scoring_sampled_pair_count: int | None
    max_exact_tree_sites: int
    max_full_candidate_scoring_sites: int
    exact_cluster_tree_built: bool
    candidate_scoring_mode: str
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None
    scale_guard_passed: bool
    candidate_scoring_evaluated: bool = False
    candidate_scoring_skip_reason: str | None = None
    candidate_scoring_applies_to: str = SIGNALOME_CANDIDATE_SCORING_APPLIES_TO
    final_module_assignment_backend: str = (
        SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE
    )
    final_module_assignment_uses_candidate_scoring: bool = False

    def __post_init__(self) -> None:
        if self.backend_diagnostics is not None:
            validate_backend_diagnostics(
                self.backend_diagnostics,
                field_name="scale_guard.backend_diagnostics",
            )
        if self.candidate_scoring_sampling is not None:
            validate_candidate_scoring_sampling_diagnostics(
                self.candidate_scoring_sampling,
                field_name="scale_guard.candidate_scoring_sampling",
            )


@dataclass(frozen=True, slots=True)
class SignalomeClusteringRunResult:
    clustering_result: ClusterSitesResult
    protein_modules: pd.Series


@dataclass(frozen=True, slots=True)
class SignalomeSupportSummary:
    kinase_substrates: dict[str, tuple[str, ...]]
    support_counts: dict[str, int]


@dataclass(frozen=True, slots=True)
class SignalomeModuleTableBuildResult:
    module_assignments: pd.DataFrame
    signalome_modules: pd.DataFrame
    module_count: int
    support_summary: SignalomeSupportSummary


@dataclass(frozen=True, slots=True)
class SignalomeNetworkBuildResult:
    edges: pd.DataFrame
    nodes: pd.DataFrame
    candidate_correlations: pd.DataFrame
    correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics


@dataclass(frozen=True, slots=True)
class SignalomeContextTableBuildResult:
    site_membership: pd.DataFrame
    protein_site_context: pd.DataFrame


__all__ = [
    "SignalomeClusteringRunResult",
    "SignalomeContextTableBuildResult",
    "SignalomeExecutionMetadata",
    "SignalomeModuleTableBuildResult",
    "SignalomeNetworkBuildResult",
    "SignalomeScaleGuardDecision",
    "SignalomeSupportSummary",
]
