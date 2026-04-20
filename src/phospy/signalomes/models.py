"""Signalome domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass
from typing import Literal

import pandas as pd

from phospy._frame_ownership import own_dataframe, own_optional_dataframe
from phospy.errors.validation import PhosPyValidationError

SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS = "correlation_thresholds"
SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT = "explicit_module_count"
SignalomeModuleSelectionStrategy = Literal[
    "correlation_thresholds",
    "explicit_module_count",
]
SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT = "allow_and_report"
SignalomeScorePreconditioningPolicy = Literal["allow_and_report"]


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


@dataclass(frozen=True, slots=True)
class SignalomeAssignments:
    """Signalome module assignment table."""

    table: pd.DataFrame
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        object.__setattr__(
            self,
            "table",
            own_dataframe(
                self.table,
                field_name="signalome_result.module_assignments.table",
                error_type=PhosPyValidationError,
                assume_owned=_assume_owned,
            ),
        )

    @classmethod
    def _from_owned(cls, *, table: pd.DataFrame) -> SignalomeAssignments:
        return cls(table=table, _assume_owned=True)


@dataclass(frozen=True, slots=True)
class SignalomeModules:
    """Signalome module table."""

    table: pd.DataFrame
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        object.__setattr__(
            self,
            "table",
            own_dataframe(
                self.table,
                field_name="signalome_result.signalome_modules.table",
                error_type=PhosPyValidationError,
                assume_owned=_assume_owned,
            ),
        )

    @classmethod
    def _from_owned(cls, *, table: pd.DataFrame) -> SignalomeModules:
        return cls(table=table, _assume_owned=True)


@dataclass(frozen=True, slots=True)
class KinaseNetwork:
    """Kinase network tables derived from signalome analysis."""

    edges: pd.DataFrame
    nodes: pd.DataFrame | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        edges = own_dataframe(
            self.edges,
            field_name="signalome_result.kinase_network.edges",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        nodes = own_optional_dataframe(
            self.nodes,
            field_name="signalome_result.kinase_network.nodes",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(self, "edges", edges)
        object.__setattr__(self, "nodes", nodes)

    @classmethod
    def _from_owned(
        cls,
        *,
        edges: pd.DataFrame,
        nodes: pd.DataFrame | None = None,
    ) -> KinaseNetwork:
        return cls(edges=edges, nodes=nodes, _assume_owned=True)


__all__ = [
    "KinaseNetwork",
    "SignalomeAssignments",
    "SignalomeClusterCandidateScore",
    "SignalomeModuleSelectionDiagnostics",
    "SignalomeModuleSelectionStrategy",
    "SignalomeScorePreconditioningDiagnostics",
    "SignalomeScorePreconditioningPolicy",
    "SignalomeModules",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT",
    "default_signalome_score_preconditioning_diagnostics",
    "SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS",
    "SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT",
    "default_signalome_module_selection_diagnostics",
]
