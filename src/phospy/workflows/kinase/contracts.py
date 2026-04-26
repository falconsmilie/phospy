"""Kinase workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.api.configs import (
    KinaseAdaptivePolicy,
    KinasePredictionMode,
    KinaseProfileMissingValueStrategy,
)
from phospy.api.requests import KinaseWorkflowRequest
from phospy.api.results import KinaseWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.references.models import ReferenceBundle
from phospy.tables.datasets import PhosphoIntensityMatrix
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference

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
    activity_phospho_matrix: pd.DataFrame
    execution_config: ResolvedKinaseExecutionConfig
    _kinase_substrate_reference: KinaseSubstrateReference = field(
        init=False,
        repr=False,
        compare=False,
    )
    _site_sequence_reference: SiteSequenceReference = field(
        init=False,
        repr=False,
        compare=False,
    )
    _activity_phospho_table: PhosphoIntensityMatrix = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        kinase_substrate_reference = KinaseSubstrateReference(
            frame=self.kinase_substrate_map,
            _assume_owned=True,
        )
        site_sequence_reference = SiteSequenceReference(
            frame=self.site_sequences,
            _assume_owned=True,
        )
        activity_phospho_table = PhosphoIntensityMatrix(
            frame=self.activity_phospho_matrix,
            allow_missing=True,
            _assume_owned=True,
        )
        if not isinstance(self.scoring_site_index, pd.Index):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.scoring_site_index_type; "
                "scoring_site_index must be a pandas Index; "
                "next_action=ensure interpreter passes a pandas Index for "
                "resolved scoring-site alignment"
            )
        if not activity_phospho_table.frame.index.equals(self.scoring_site_index):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.activity_site_alignment; "
                "activity_phospho_matrix.index must exactly match scoring_site_index; "
                "next_action=ensure interpreted activity phospho rows are aligned "
                "to the resolved scoring-site index"
            )
        object.__setattr__(
            self, "kinase_substrate_map", kinase_substrate_reference.frame
        )
        object.__setattr__(self, "site_sequences", site_sequence_reference.frame)
        object.__setattr__(
            self, "activity_phospho_matrix", activity_phospho_table.frame
        )
        object.__setattr__(
            self,
            "_kinase_substrate_reference",
            kinase_substrate_reference,
        )
        object.__setattr__(
            self,
            "_site_sequence_reference",
            site_sequence_reference,
        )
        object.__setattr__(
            self,
            "_activity_phospho_table",
            activity_phospho_table,
        )

    @property
    def kinase_substrate_reference(self) -> KinaseSubstrateReference:
        return self._kinase_substrate_reference

    @property
    def site_sequence_reference(self) -> SiteSequenceReference:
        return self._site_sequence_reference

    @property
    def activity_phospho_table(self) -> PhosphoIntensityMatrix:
        return self._activity_phospho_table


@dataclass(frozen=True, slots=True)
class ResolvedKinaseActivityExecutionConfig:
    """Execution-ready kinase activity-stage config."""

    threshold: float
    min_substrates: int
    top_n_substrates: int


@dataclass(frozen=True, slots=True)
class ResolvedKinaseExecutionConfig:
    """Execution-ready kinase workflow config resolved by the interpreter."""

    scoring_min_substrates: int
    include_diagnostic_scoring_tables: bool
    profile_missing_value_strategy: KinaseProfileMissingValueStrategy
    prediction_top_k: int
    prediction_deterministic_max_selected_kinases: int
    prediction_adaptive_ensemble_runs: int
    prediction_mode: KinasePredictionMode
    prediction_adaptive_policy: KinaseAdaptivePolicy
    prediction_n_iterations: int
    prediction_random_state: int | None
    activity: ResolvedKinaseActivityExecutionConfig | None


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
    "ResolvedKinaseActivityExecutionConfig",
    "ResolvedKinaseExecutionConfig",
    "ResolvedKinaseWorkflowRequest",
]
