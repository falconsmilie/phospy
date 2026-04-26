"""Signalome workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

import pandas as pd

from phospy.api.configs import (
    SignalomeAssignmentPolicy,
    SignalomeClusteringBackend,
    SignalomeKinaseNetworkPolicy,
    SignalomeScorePreconditioningPolicy,
)
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.signalomes.models import (
    SignalomeScorePreconditioningDiagnostics,
    default_signalome_score_preconditioning_diagnostics,
)
from phospy.tables.kinase import KinasePredictionMatrix, KinaseScoreMatrix


@dataclass(frozen=True, slots=True)
class ResolvedSignalomeExecutionConfig:
    """Execution-ready signalome config resolved by the interpreter."""

    substrate_support_cutoff: float
    network_correlation_threshold: float
    network_policy: SignalomeKinaseNetworkPolicy
    assignment_policy: SignalomeAssignmentPolicy
    score_preconditioning_policy: SignalomeScorePreconditioningPolicy
    module_selection_primary_threshold: float
    module_selection_fallback_threshold: float
    module_selection_max_clusters: int
    clustering_backend: SignalomeClusteringBackend
    max_exact_clustering_sites: int
    requested_module_count: int | None


@dataclass(frozen=True, slots=True)
class ResolvedSignalomeWorkflowRequest:
    """Interpreter output for signalome workflow execution.

    ``site_to_protein`` must provide a non-empty explicit protein identifier from
    ``dataset.site_metadata.protein_id`` for every site in ``prediction_matrix.index``.
    ``downstream_score_matrix`` is the same authoritative matrix lane that drove
    upstream kinase prediction, after interpreter preconditioning of unsupported
    all-missing score rows. ``score_preconditioning_diagnostics`` surfaces the
    aligned input row count, dropped all-missing row count, retained row count,
    and active `SignalomeConfig.score_preconditioning_policy`.
    """

    dataset: AnalysisReadyPhosphoDataset
    kinase_result: KinaseWorkflowResult
    execution_config: ResolvedSignalomeExecutionConfig
    downstream_score_matrix: pd.DataFrame
    downstream_score_source: str
    prediction_matrix: pd.DataFrame
    site_to_protein: pd.Series
    score_preconditioning_diagnostics: SignalomeScorePreconditioningDiagnostics = field(
        default_factory=default_signalome_score_preconditioning_diagnostics
    )
    _downstream_score_table: KinaseScoreMatrix = field(
        init=False,
        repr=False,
        compare=False,
    )
    _prediction_table: KinasePredictionMatrix = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        downstream_score_table = KinaseScoreMatrix(
            frame=self.downstream_score_matrix,
            field_name="signalome_request.downstream_score_matrix",
            _assume_owned=True,
        )
        prediction_table = KinasePredictionMatrix(
            frame=self.prediction_matrix,
            field_name="signalome_request.prediction_matrix",
            _assume_owned=True,
        )
        if not isinstance(self.site_to_protein, pd.Series):
            raise WorkflowBoundaryError(
                "signalome workflow boundary validation failed at seam="
                "signalome.contracts.site_to_protein_type; "
                "site_to_protein must be a pandas Series; "
                "next_action=ensure signalome interpreter resolves an explicit "
                "site-to-protein mapping series"
            )
        object.__setattr__(
            self, "downstream_score_matrix", downstream_score_table.frame
        )
        object.__setattr__(self, "prediction_matrix", prediction_table.frame)
        object.__setattr__(self, "_downstream_score_table", downstream_score_table)
        object.__setattr__(self, "_prediction_table", prediction_table)

    @property
    def downstream_score_table(self) -> KinaseScoreMatrix:
        return self._downstream_score_table

    @property
    def prediction_table(self) -> KinasePredictionMatrix:
        return self._prediction_table


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
