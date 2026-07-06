"""Signalome workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

import pandas as pd

from phospy.contracts.configs import (
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT,
    LocalisationRequirement,
    SignalomeAssignmentPolicy,
    SignalomeCandidateScoringPolicy,
    SignalomeKinaseNetworkPolicy,
    SignalomeScorePreconditioningPolicy,
)
from phospy.contracts.requests import SignalomeWorkflowRequest
from phospy.contracts.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.scoring import DownstreamScoreSelectionPolicy
from phospy.science.scoring.policy_models import DownstreamScoreSource
from phospy.science.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
)
from phospy.science.signalomes.clustering.policies import (
    SignalomeCandidateScoringPolicyDefinition,
)
from phospy.science.signalomes.models import (
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
    allow_mixed_total_protein_quantitative_meaning: bool
    module_selection_primary_threshold: float
    module_selection_fallback_threshold: float
    module_selection_max_clusters: int
    candidate_scoring_policy: SignalomeCandidateScoringPolicy
    max_exact_tree_sites: int
    max_full_candidate_scoring_sites: int
    requested_module_count: int | None
    clustering_engine: str = SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    network_min_paired_finite_observations: int = (
        SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT
    )
    candidate_scoring_policy_definition: (
        SignalomeCandidateScoringPolicyDefinition | None
    ) = None
    score_preconditioning_policy_definition: ScientificPolicyRecord | None = None
    localisation_requirement: LocalisationRequirement = field(
        default_factory=LocalisationRequirement
    )


@dataclass(frozen=True, slots=True)
class ResolvedSignalomeWorkflowRequest:
    """Interpreter output for signalome workflow execution.

    ``site_to_protein`` must provide a non-empty signalome protein grouping
    label from ``dataset.site_metadata.protein_id`` for every site in
    ``prediction_matrix.index``.
    ``downstream_score_matrix`` is the same authoritative matrix lane that drove
    upstream kinase prediction, after interpreter preconditioning of unsupported
    all-missing score rows. ``score_preconditioning_diagnostics`` surfaces the
    aligned input row count, dropped all-missing row count, retained row count,
    and active ``SignalomeConfig.validation.score_preconditioning_policy``.
    ``alignment_diagnostics`` reports provided/retained/dropped counts (and
    exclusion reasons) for scientific input alignment across sites, kinases,
    and protein identifiers.
    """

    dataset: AnalysisReadyPhosphoDataset
    kinase_result: KinaseWorkflowResult
    execution_config: ResolvedSignalomeExecutionConfig
    downstream_score_matrix: pd.DataFrame
    downstream_score_source: DownstreamScoreSource
    prediction_matrix: pd.DataFrame
    site_to_protein: pd.Series
    score_preconditioning_diagnostics: SignalomeScorePreconditioningDiagnostics = field(
        default_factory=default_signalome_score_preconditioning_diagnostics
    )
    alignment_diagnostics: SignalomeAlignmentDiagnostics = field(
        default_factory=default_signalome_alignment_diagnostics
    )
    downstream_score_selection_policy: DownstreamScoreSelectionPolicy | None = None
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
        downstream_score_source = DownstreamScoreSource.parse(
            self.downstream_score_source,
            field_name="signalome_request.downstream_score_source",
        )
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
        object.__setattr__(self, "downstream_score_source", downstream_score_source)
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
