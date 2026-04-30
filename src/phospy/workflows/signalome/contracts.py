"""Signalome workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

import pandas as pd

from phospy.api.configs import (
    SignalomeAssignmentPolicy,
    SignalomeCandidateScoringPolicy,
    SignalomeKinaseNetworkPolicy,
    SignalomeScorePreconditioningPolicy,
    SignalomeTreeEngine,
)
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
)
from phospy.signalomes.models import (
    SignalomeAlignmentDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
    default_signalome_alignment_diagnostics,
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
    tree_engine: SignalomeTreeEngine
    candidate_scoring_policy: SignalomeCandidateScoringPolicy
    max_exact_tree_sites: int
    max_full_candidate_scoring_sites: int
    requested_module_count: int | None
    clustering_engine: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON


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
    ``alignment_diagnostics`` reports provided/retained/dropped counts (and
    exclusion reasons) for scientific input alignment across sites, kinases,
    and protein identifiers.
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
    alignment_diagnostics: SignalomeAlignmentDiagnostics = field(
        default_factory=default_signalome_alignment_diagnostics
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
        downstream_site_index = downstream_score_table.frame.index
        prediction_site_index = prediction_table.frame.index
        site_to_protein_index = self.site_to_protein.index
        if not prediction_site_index.equals(
            downstream_site_index
        ) or not site_to_protein_index.equals(downstream_site_index):
            raise WorkflowBoundaryError(
                "signalome workflow boundary validation failed at seam="
                "signalome.contracts.site_index_alignment; "
                "downstream_score_matrix.index, prediction_matrix.index, and "
                "site_to_protein.index must match exactly; "
                f"downstream_score_sites={int(downstream_site_index.size)}; "
                f"prediction_sites={int(prediction_site_index.size)}; "
                f"site_to_protein_sites={int(site_to_protein_index.size)}; "
                "next_action=ensure interpreter aligns prediction_matrix and "
                "site_to_protein to retained downstream score sites after "
                "score preconditioning"
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

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowRequest: ...


class SignalomeWorkflowInterpreterContract(Protocol):
    """Internal contract for signalome workflow request interpretation."""

    def run(
        self, request: SignalomeWorkflowRequest
    ) -> ResolvedSignalomeWorkflowRequest: ...


class SignalomeWorkflowExecutorContract(Protocol):
    """Internal contract for signalome workflow execution."""

    def run(
        self, request: ResolvedSignalomeWorkflowRequest
    ) -> SignalomeWorkflowResult: ...


__all__ = [
    "ResolvedSignalomeExecutionConfig",
    "ResolvedSignalomeWorkflowRequest",
    "SignalomeWorkflowExecutorContract",
    "SignalomeWorkflowInterpreterContract",
    "SignalomeWorkflowValidatorContract",
]
