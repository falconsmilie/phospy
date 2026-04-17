"""Signalome workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.api.configs import SignalomeConfig
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import SignalomeWorkflowResult, SimpleKinaseWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset


@dataclass(frozen=True, slots=True)
class ResolvedSignalomeWorkflowRequest:
    """Interpreter output for signalome workflow execution."""

    dataset: AnalysisReadyPhosphoDataset
    kinase_result: SimpleKinaseWorkflowResult
    config: SignalomeConfig
    score_matrix: pd.DataFrame
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
    "ResolvedSignalomeWorkflowRequest",
    "SignalomeWorkflowExecutorContract",
    "SignalomeWorkflowInterpreterContract",
    "SignalomeWorkflowValidatorContract",
]
