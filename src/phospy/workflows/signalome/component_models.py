"""Shared component dataclasses for signalome workflow execution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.signalomes.clustering import (
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    ClusterSitesResult,
)
from phospy.signalomes.models import SignalomeNetworkCorrelationDiagnostics


@dataclass(frozen=True, slots=True)
class SignalomeExecutionMetadata:
    prediction_sites: int
    prediction_kinases: int
    downstream_score_sites: int
    downstream_score_kinases: int
    downstream_score_source: str


@dataclass(frozen=True, slots=True)
class SignalomeScaleGuardDecision:
    site_count: int
    cluster_tree_backend: str
    candidate_scoring_backend: str
    max_exact_cluster_tree_sites: int
    max_full_correlation_sites: int
    exact_cluster_tree_built: bool
    candidate_scoring_mode: str
    candidate_scoring_sampling: dict[str, object] | None
    scale_guard_passed: bool
    candidate_scoring_evaluated: bool = False
    candidate_scoring_skip_reason: str | None = None
    candidate_scoring_applies_to: str = SIGNALOME_CANDIDATE_SCORING_APPLIES_TO
    final_module_assignment_backend: str = (
        SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE
    )
    final_module_assignment_uses_candidate_scoring: bool = False


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
