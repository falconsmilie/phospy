"""Signalome domain models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import pandas as pd

from phospy.contracts.configs import (
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    SignalomeScorePreconditioningPolicy,
)
from phospy.errors.validation import PhosPyValidationError
from phospy.frames.ownership import (
    export_dataframe,
    export_optional_dataframe,
    own_dataframe,
    own_optional_dataframe,
)

SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS = "correlation_thresholds"
SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT = "explicit_module_count"
SignalomeModuleSelectionStrategy = Literal[
    "correlation_thresholds",
    "explicit_module_count",
]
SIGNALOME_CORRELATION_STATUS_FINITE = "finite"
SIGNALOME_CORRELATION_STATUS_CONSTANT_PROFILE = "constant_profile"
SIGNALOME_CORRELATION_STATUS_INSUFFICIENT_OBSERVATIONS = "insufficient_observations"
SIGNALOME_CORRELATION_STATUS_MISSING_VALUES = "missing_values"
SIGNALOME_CORRELATION_STATUS_NON_FINITE_VALUES = "non_finite_values"
SIGNALOME_CORRELATION_STATUS_UNDEFINED = "undefined"
SignalomeCorrelationStatus = Literal[
    "finite",
    "constant_profile",
    "insufficient_observations",
    "missing_values",
    "non_finite_values",
    "undefined",
]


@dataclass(frozen=True, slots=True)
class SignalomeClusterCandidateScore:
    """One candidate module-count score summary."""

    min_median_correlation: float
    mean_median_correlation: float


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionDiagnostics:
    """Structured module-selection diagnostics."""

    strategy: SignalomeModuleSelectionStrategy
    selected_module_count: int
    requested_module_count: int | None
    threshold_used: float | None
    max_clusters_evaluated: int
    candidate_scores: dict[int, SignalomeClusterCandidateScore]
    reason: str
    zero_variance_profile_count: int = 0
    near_constant_profile_count: int = 0
    excluded_from_correlation_count: int = 0

    @property
    def used_automatic_selection(self) -> bool:
        return self.requested_module_count is None


@dataclass(frozen=True, slots=True)
class SignalomeScorePreconditioningDiagnostics:
    """Structured score preconditioning diagnostics for signalome inputs."""

    input_row_count: int
    dropped_all_missing_row_count: int
    retained_row_count: int
    policy: SignalomeScorePreconditioningPolicy = (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )


@dataclass(frozen=True, slots=True)
class SignalomeAlignmentInputDiagnostics:
    """Structured provided/retained/dropped counts for one aligned input lane."""

    provided_count: int
    retained_count: int
    dropped_count: int
    dropped_reasons: dict[str, int] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class SignalomeAlignmentDiagnostics:
    """Structured alignment diagnostics across signalome scientific inputs."""

    dataset_sites: SignalomeAlignmentInputDiagnostics
    prediction_score_sites: SignalomeAlignmentInputDiagnostics
    downstream_score_sites: SignalomeAlignmentInputDiagnostics
    kinases: SignalomeAlignmentInputDiagnostics
    protein_identifiers: SignalomeAlignmentInputDiagnostics


@dataclass(frozen=True, slots=True)
class SignalomeNetworkCorrelationDiagnostics:
    """Structured diagnostics for kinase-network correlation eligibility."""

    total_candidate_correlations: int
    finite_correlations: int
    undefined_correlations: int
    constant_profile_correlations: int
    insufficient_observation_correlations: int
    missing_value_correlations: int
    non_finite_value_correlations: int
    edges_created: int
    edges_skipped_non_finite_correlation: int


def default_signalome_module_selection_diagnostics() -> (
    SignalomeModuleSelectionDiagnostics
):
    """Return a stable placeholder diagnostics payload."""

    return SignalomeModuleSelectionDiagnostics(
        strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        selected_module_count=1,
        requested_module_count=None,
        threshold_used=None,
        max_clusters_evaluated=1,
        candidate_scores={},
        reason="module selection diagnostics were not captured",
    )


def default_signalome_score_preconditioning_diagnostics() -> (
    SignalomeScorePreconditioningDiagnostics
):
    """Return a stable placeholder score-preconditioning diagnostics payload."""

    return SignalomeScorePreconditioningDiagnostics(
        input_row_count=0,
        dropped_all_missing_row_count=0,
        retained_row_count=0,
        policy=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    )


def default_signalome_alignment_diagnostics() -> SignalomeAlignmentDiagnostics:
    """Return stable placeholder alignment diagnostics payload."""

    def _empty() -> SignalomeAlignmentInputDiagnostics:
        return SignalomeAlignmentInputDiagnostics(
            provided_count=0,
            retained_count=0,
            dropped_count=0,
            dropped_reasons={},
        )

    return SignalomeAlignmentDiagnostics(
        dataset_sites=_empty(),
        prediction_score_sites=_empty(),
        downstream_score_sites=_empty(),
        kinases=_empty(),
        protein_identifiers=_empty(),
    )


def default_signalome_network_correlation_diagnostics() -> (
    SignalomeNetworkCorrelationDiagnostics
):
    """Return a stable placeholder network-correlation diagnostics payload."""

    return SignalomeNetworkCorrelationDiagnostics(
        total_candidate_correlations=0,
        finite_correlations=0,
        undefined_correlations=0,
        constant_profile_correlations=0,
        insufficient_observation_correlations=0,
        missing_value_correlations=0,
        non_finite_value_correlations=0,
        edges_created=0,
        edges_skipped_non_finite_correlation=0,
    )


@dataclass(frozen=True, slots=True, init=False)
class SignalomeAssignments:
    """Score-derived signalome module assignment table."""

    _table: pd.DataFrame = field(init=False, repr=False)

    def __init__(self, table: pd.DataFrame, _assume_owned: bool = False) -> None:
        from phospy.tables.signalome import SignalomeAssignmentsTable

        table = own_dataframe(
            table,
            field_name="signalome_result.module_assignments.table",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(
            self,
            "_table",
            SignalomeAssignmentsTable(frame=table, _assume_owned=True).frame,
        )

    @classmethod
    def _from_owned(cls, *, table: pd.DataFrame) -> SignalomeAssignments:
        return cls(table=table, _assume_owned=True)

    def to_pandas(self) -> pd.DataFrame:
        """Return an assignments snapshot isolated from this object."""

        return export_dataframe(self._table)

    @property
    def table(self) -> pd.DataFrame:
        return export_dataframe(self._table)


@dataclass(frozen=True, slots=True, init=False)
class SignalomeModules:
    """Candidate kinase-supported module summary table."""

    _table: pd.DataFrame = field(init=False, repr=False)

    def __init__(self, table: pd.DataFrame, _assume_owned: bool = False) -> None:
        from phospy.tables.signalome import SignalomeModulesTable

        table = own_dataframe(
            table,
            field_name="signalome_result.signalome_modules.table",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(
            self,
            "_table",
            SignalomeModulesTable(frame=table, _assume_owned=True).frame,
        )

    @classmethod
    def _from_owned(cls, *, table: pd.DataFrame) -> SignalomeModules:
        return cls(table=table, _assume_owned=True)

    def to_pandas(self) -> pd.DataFrame:
        """Return a modules snapshot isolated from this object."""

        return export_dataframe(self._table)

    @property
    def table(self) -> pd.DataFrame:
        return export_dataframe(self._table)


@dataclass(frozen=True, slots=True, init=False)
class KinaseNetwork:
    """Kinase score-profile association tables derived from signalome analysis."""

    correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics = field(
        default_factory=default_signalome_network_correlation_diagnostics
    )
    _edges: pd.DataFrame = field(init=False, repr=False)
    _nodes: pd.DataFrame | None = field(init=False, repr=False)
    _candidate_correlations: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        edges: pd.DataFrame,
        nodes: pd.DataFrame | None = None,
        candidate_correlations: pd.DataFrame | None = None,
        correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics | None = None,
        _assume_owned: bool = False,
    ) -> None:
        from phospy.tables.signalome import (
            KinaseNetworkCandidateCorrelationsTable,
            KinaseNetworkEdgesTable,
            KinaseNetworkNodesTable,
        )

        object.__setattr__(
            self,
            "correlation_diagnostics",
            (
                default_signalome_network_correlation_diagnostics()
                if correlation_diagnostics is None
                else correlation_diagnostics
            ),
        )

        edges = own_dataframe(
            edges,
            field_name="signalome_result.kinase_network.edges",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        nodes = own_optional_dataframe(
            nodes,
            field_name="signalome_result.kinase_network.nodes",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        candidate_correlations = own_optional_dataframe(
            candidate_correlations,
            field_name="signalome_result.kinase_network.candidate_correlations",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        if not isinstance(
            self.correlation_diagnostics,
            SignalomeNetworkCorrelationDiagnostics,
        ):
            raise PhosPyValidationError(
                "signalome_result.kinase_network.correlation_diagnostics must be "
                "SignalomeNetworkCorrelationDiagnostics"
            )
        edges = KinaseNetworkEdgesTable(frame=edges, _assume_owned=True).frame
        if nodes is not None:
            nodes = KinaseNetworkNodesTable(frame=nodes, _assume_owned=True).frame
        if candidate_correlations is not None:
            candidate_correlations = KinaseNetworkCandidateCorrelationsTable(
                frame=candidate_correlations,
                _assume_owned=True,
            ).frame
        object.__setattr__(self, "_edges", edges)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_candidate_correlations", candidate_correlations)

    @property
    def edges(self) -> pd.DataFrame:
        return export_dataframe(self._edges)

    @property
    def nodes(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._nodes)

    @property
    def candidate_correlations(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._candidate_correlations)

    @classmethod
    def _from_owned(
        cls,
        *,
        edges: pd.DataFrame,
        nodes: pd.DataFrame | None = None,
        candidate_correlations: pd.DataFrame | None = None,
        correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics | None = None,
    ) -> KinaseNetwork:
        return cls(
            edges=edges,
            nodes=nodes,
            candidate_correlations=candidate_correlations,
            correlation_diagnostics=correlation_diagnostics,
            _assume_owned=True,
        )

    def to_pandas(self) -> pd.DataFrame:
        """Return an edges snapshot isolated from this network."""

        return export_dataframe(self._edges)

    def nodes_dataframe(self) -> pd.DataFrame | None:
        """Return an optional nodes snapshot isolated from this network."""

        return export_optional_dataframe(self._nodes)

    def candidate_correlations_dataframe(self) -> pd.DataFrame | None:
        """Return optional candidate correlations isolated from this network."""

        return export_optional_dataframe(self._candidate_correlations)


__all__ = [
    "KinaseNetwork",
    "SignalomeAlignmentDiagnostics",
    "SignalomeAlignmentInputDiagnostics",
    "SIGNALOME_CORRELATION_STATUS_CONSTANT_PROFILE",
    "SIGNALOME_CORRELATION_STATUS_FINITE",
    "SIGNALOME_CORRELATION_STATUS_INSUFFICIENT_OBSERVATIONS",
    "SIGNALOME_CORRELATION_STATUS_MISSING_VALUES",
    "SIGNALOME_CORRELATION_STATUS_NON_FINITE_VALUES",
    "SIGNALOME_CORRELATION_STATUS_UNDEFINED",
    "SignalomeAssignments",
    "SignalomeClusterCandidateScore",
    "SignalomeCorrelationStatus",
    "SignalomeModuleSelectionDiagnostics",
    "SignalomeModuleSelectionStrategy",
    "SignalomeNetworkCorrelationDiagnostics",
    "SignalomeScorePreconditioningDiagnostics",
    "SignalomeScorePreconditioningPolicy",
    "SignalomeModules",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP",
    "default_signalome_alignment_diagnostics",
    "default_signalome_score_preconditioning_diagnostics",
    "SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS",
    "SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT",
    "default_signalome_network_correlation_diagnostics",
    "default_signalome_module_selection_diagnostics",
]
