"""Activity method diagnostics value models and assembly."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.frames.comparison import optional_dataframe_equals
from phospy.frames.ownership import export_optional_dataframe
from phospy.science.activities.method_models import (
    SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD,
    SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD,
    ActivityMethodMetadata,
    ActivityMethodSummary,
)
from phospy.science.activities.threshold_membership import (
    ActivityThresholdMembershipDiagnostics,
)
from phospy.science.tables.activity import ActivityStatisticsTable


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ActivityMethodDiagnostics:
    """Typed method diagnostics carried alongside activity result matrices.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit diagnostics-content comparison.
    """

    __hash__ = object.__hash__

    method_summary: ActivityMethodSummary | None
    threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics | None
    _statistics_table: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        *,
        method_summary: ActivityMethodSummary | None = None,
        threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics
        | None = None,
        statistics_table: pd.DataFrame | None = None,
        _assume_owned: bool = False,
    ) -> None:
        _validate_optional_method_summary(method_summary)
        _validate_optional_threshold_membership_diagnostics(
            threshold_membership_diagnostics
        )
        if statistics_table is not None:
            statistics_table = ActivityStatisticsTable(
                frame=statistics_table,
                _assume_owned=_assume_owned,
            ).frame
        object.__setattr__(self, "method_summary", method_summary)
        object.__setattr__(
            self,
            "threshold_membership_diagnostics",
            threshold_membership_diagnostics,
        )
        object.__setattr__(self, "_statistics_table", statistics_table)

    @property
    def statistics_table(self) -> pd.DataFrame | None:
        """Return an optional statistics-table snapshot, not a report export."""

        return export_optional_dataframe(self._statistics_table)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another diagnostics object has the same content."""

        if not isinstance(other, ActivityMethodDiagnostics):
            return False
        return (
            self.method_summary == other.method_summary
            and self.threshold_membership_diagnostics
            == other.threshold_membership_diagnostics
            and optional_dataframe_equals(
                self._statistics_table,
                other._statistics_table,
            )
        )


class WeightedSubstrateActivityDiagnostics(ActivityMethodDiagnostics):
    """Diagnostics for simplified weighted substrate activity-like scores."""

    __slots__ = ()


class KseaZScoreActivityDiagnostics(ActivityMethodDiagnostics):
    """Diagnostics for KSEA-style kinase activity score outputs."""

    __slots__ = ()


class SsgseaSubstrateEnrichmentActivityDiagnostics(ActivityMethodDiagnostics):
    """Diagnostics for ssGSEA-style substrate enrichment activity-like scores."""

    __slots__ = ()


def _validate_optional_method_summary(
    method_summary: ActivityMethodSummary | None,
) -> None:
    if method_summary is not None and not isinstance(
        method_summary,
        ActivityMethodSummary,
    ):
        raise WorkflowBoundaryError(
            "activity_result.method_summary must be ActivityMethodSummary or None"
        )


def _validate_optional_threshold_membership_diagnostics(
    threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics | None,
) -> None:
    if threshold_membership_diagnostics is not None and not isinstance(
        threshold_membership_diagnostics,
        ActivityThresholdMembershipDiagnostics,
    ):
        raise WorkflowBoundaryError(
            "activity_result.threshold_membership_diagnostics must be "
            "ActivityThresholdMembershipDiagnostics or None"
        )


def _build_activity_method_diagnostics(
    *,
    activity_method: ActivityMethodMetadata,
    method_summary: ActivityMethodSummary | None,
    threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics | None,
    statistics_table: pd.DataFrame | None,
) -> ActivityMethodDiagnostics:
    diagnostics_cls: type[ActivityMethodDiagnostics]
    if activity_method.is_ksea:
        diagnostics_cls = KseaZScoreActivityDiagnostics
    elif (
        activity_method.activity_method_id
        == SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD.activity_method_id
    ):
        diagnostics_cls = SsgseaSubstrateEnrichmentActivityDiagnostics
    elif (
        activity_method.activity_method_id
        == SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD.activity_method_id
    ):
        diagnostics_cls = WeightedSubstrateActivityDiagnostics
    else:
        diagnostics_cls = ActivityMethodDiagnostics
    return diagnostics_cls(
        method_summary=method_summary,
        threshold_membership_diagnostics=threshold_membership_diagnostics,
        statistics_table=statistics_table,
        _assume_owned=True,
    )


__all__ = [
    "ActivityMethodDiagnostics",
    "KseaZScoreActivityDiagnostics",
    "SsgseaSubstrateEnrichmentActivityDiagnostics",
    "WeightedSubstrateActivityDiagnostics",
]
