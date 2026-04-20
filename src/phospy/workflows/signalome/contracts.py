"""Signalome workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.api.configs import SignalomeAssignmentPolicy, SignalomeKinaseNetworkPolicy
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset


@dataclass(frozen=True, slots=True)
class ResolvedSignalomeExecutionConfig:
    """Execution-ready signalome config resolved by the interpreter."""

    substrate_support_cutoff: float
    network_correlation_threshold: float
    network_policy: SignalomeKinaseNetworkPolicy
    assignment_policy: SignalomeAssignmentPolicy
    module_selection_primary_threshold: float
    module_selection_fallback_threshold: float
    module_selection_max_clusters: int
    requested_module_count: int | None


@dataclass(frozen=True, slots=True)
class ResolvedSignalomeWorkflowRequest:
    """Interpreter output for signalome workflow execution.

    ``site_to_protein`` must provide a non-empty explicit protein identifier from
    ``dataset.site_metadata.protein_id`` for every site in ``prediction_matrix.index``.
    ``downstream_score_matrix`` is the same authoritative matrix lane that drove
    upstream kinase prediction, after interpreter preconditioning of unsupported
    all-missing score rows.
    """

    dataset: AnalysisReadyPhosphoDataset
    kinase_result: KinaseWorkflowResult
    execution_config: ResolvedSignalomeExecutionConfig
    downstream_score_matrix: pd.DataFrame
    downstream_score_source: str
    prediction_matrix: pd.DataFrame
    site_to_protein: pd.Series


class SignalomeWorkflowValidatorContract(Protocol):
    """Internal contract for signalome workflow request validation."""

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowRequest:
        """Validate the workflow request and return the same request."""


class SignalomeWorkflowInterpreterContract(Protocol):
    """Internal contract for signalome workflow request interpretation."""

    def run(
        self, request: SignalomeWorkflowRequest
    ) -> ResolvedSignalomeWorkflowRequest:
        """Resolve runtime defaults into execution-ready signalome inputs."""


class SignalomeWorkflowExecutorContract(Protocol):
    """Internal contract for signalome workflow execution."""

    def run(self, request: ResolvedSignalomeWorkflowRequest) -> SignalomeWorkflowResult:
        """Execute signalome domain logic and assemble public results."""


__all__ = [
    "ResolvedSignalomeExecutionConfig",
    "ResolvedSignalomeWorkflowRequest",
    "SignalomeWorkflowExecutorContract",
    "SignalomeWorkflowInterpreterContract",
    "SignalomeWorkflowValidatorContract",
]
