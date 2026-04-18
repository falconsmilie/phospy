"""Kinase workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
)
from phospy.api.requests import KinaseWorkflowRequest
from phospy.api.results import KinaseWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.references.models import ReferenceBundle

if TYPE_CHECKING:
    from phospy.references.resolution import ReferenceResolverContract


@dataclass(frozen=True, slots=True)
class ResolvedKinaseWorkflowRequest:
    """Interpreter output for kinase workflow execution."""

    dataset: AnalysisReadyPhosphoDataset
    references: ReferenceBundle
    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
    scoring_site_index: pd.Index
    scoring_config: KinaseScoringConfig
    prediction_config: KinasePredictionConfig
    activity_config: KinaseActivityConfig | None


class KinaseWorkflowValidatorContract(Protocol):
    """Internal contract for kinase workflow request validation."""

    def run(self, request: KinaseWorkflowRequest) -> KinaseWorkflowRequest:
        """Validate the workflow request and return the same request."""


class KinaseWorkflowInterpreterContract(Protocol):
    """Internal contract for kinase workflow request interpretation."""

    _reference_resolver: ReferenceResolverContract

    def run(self, request: KinaseWorkflowRequest) -> ResolvedKinaseWorkflowRequest:
        """Resolve references and runtime defaults for execution."""


class KinaseWorkflowExecutorContract(Protocol):
    """Internal contract for kinase workflow execution."""

    def run(self, request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
        """Execute workflow domain logic and assemble public results."""


__all__ = [
    "KinaseWorkflowExecutorContract",
    "KinaseWorkflowInterpreterContract",
    "KinaseWorkflowValidatorContract",
    "ResolvedKinaseWorkflowRequest",
]
