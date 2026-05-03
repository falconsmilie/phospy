"""Activity stage result models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import pandas as pd

from phospy._frame_ownership import (
    export_dataframe,
    export_optional_dataframe,
    export_series,
)
from phospy.errors.validation import PhosPyValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.tables.activity import (
    ActivityCountSeries,
    ActivityMatrix,
    ActivityStatisticsTable,
    ActivityTargetTable,
)
from phospy.tables.kinase import KinasePredictionMatrix


@dataclass(frozen=True, slots=True)
class ActivityMethodMetadata:
    """Stable scientific identity metadata for an activity scoring method."""

    activity_method_id: str
    activity_method_family: str
    activity_method_label: str
    is_ksea: bool
    is_phosr_kinase_activity_equivalent: bool

    def to_payload(self) -> dict[str, object]:
        return {
            "activity_method_id": self.activity_method_id,
            "activity_method_family": self.activity_method_family,
            "activity_method_label": self.activity_method_label,
            "is_ksea": bool(self.is_ksea),
            "is_phosr_kinase_activity_equivalent": bool(
                self.is_phosr_kinase_activity_equivalent
            ),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMethodMetadata:
        method_id = str(payload.get("activity_method_id", "")).strip()
        method_family = str(payload.get("activity_method_family", "")).strip()
        method_label = str(payload.get("activity_method_label", "")).strip()
        is_ksea = payload.get("is_ksea")
        is_phosr_equivalent = payload.get("is_phosr_kinase_activity_equivalent")
        if not method_id:
            raise ValueError("activity_method_id must be a non-empty string")
        if not method_family:
            raise ValueError("activity_method_family must be a non-empty string")
        if not method_label:
            raise ValueError("activity_method_label must be a non-empty string")
        if not isinstance(is_ksea, bool):
            raise ValueError("is_ksea must be a bool")
        if not isinstance(is_phosr_equivalent, bool):
            raise ValueError("is_phosr_kinase_activity_equivalent must be a bool")
        return cls(
            activity_method_id=method_id,
            activity_method_family=method_family,
            activity_method_label=method_label,
            is_ksea=is_ksea,
            is_phosr_kinase_activity_equivalent=is_phosr_equivalent,
        )


SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD = ActivityMethodMetadata(
    activity_method_id="simplified_weighted_substrate_activity_v1",
    activity_method_family="heuristic_weighted_substrate_score",
    activity_method_label="simplified weighted substrate activity",
    is_ksea=False,
    is_phosr_kinase_activity_equivalent=False,
)

KSEA_ZSCORE_ACTIVITY_METHOD = ActivityMethodMetadata(
    activity_method_id="ksea_zscore_v1",
    activity_method_family="substrate_set_enrichment",
    activity_method_label="ksea z-score activity",
    is_ksea=True,
    is_phosr_kinase_activity_equivalent=False,
)


@dataclass(frozen=True, slots=True)
class ActivityMethodSummary:
    """Method-level score computability counters."""

    kinases_evaluated: int
    kinase_condition_pairs_evaluated: int
    kinase_condition_pairs_computed: int
    kinase_condition_pairs_insufficient_substrates: int
    kinase_condition_pairs_invalid_background_variance: int
    kinase_condition_pairs_no_finite_background_values: int
    kinase_condition_pairs_no_finite_substrate_values: int

    def to_payload(self) -> dict[str, int]:
        return {
            "kinases_evaluated": int(self.kinases_evaluated),
            "kinase_condition_pairs_evaluated": int(
                self.kinase_condition_pairs_evaluated
            ),
            "kinase_condition_pairs_computed": int(
                self.kinase_condition_pairs_computed
            ),
            "kinase_condition_pairs_insufficient_substrates": int(
                self.kinase_condition_pairs_insufficient_substrates
            ),
            "kinase_condition_pairs_invalid_background_variance": int(
                self.kinase_condition_pairs_invalid_background_variance
            ),
            "kinase_condition_pairs_no_finite_background_values": int(
                self.kinase_condition_pairs_no_finite_background_values
            ),
            "kinase_condition_pairs_no_finite_substrate_values": int(
                self.kinase_condition_pairs_no_finite_substrate_values
            ),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMethodSummary:
        return cls(
            kinases_evaluated=int(payload.get("kinases_evaluated", 0)),
            kinase_condition_pairs_evaluated=int(
                payload.get("kinase_condition_pairs_evaluated", 0)
            ),
            kinase_condition_pairs_computed=int(
                payload.get("kinase_condition_pairs_computed", 0)
            ),
            kinase_condition_pairs_insufficient_substrates=int(
                payload.get("kinase_condition_pairs_insufficient_substrates", 0)
            ),
            kinase_condition_pairs_invalid_background_variance=int(
                payload.get("kinase_condition_pairs_invalid_background_variance", 0)
            ),
            kinase_condition_pairs_no_finite_background_values=int(
                payload.get("kinase_condition_pairs_no_finite_background_values", 0)
            ),
            kinase_condition_pairs_no_finite_substrate_values=int(
                payload.get("kinase_condition_pairs_no_finite_substrate_values", 0)
            ),
        )


@dataclass(frozen=True, slots=True)
class PredMatOverlapSummary:
    """Resolved overlap diagnostics between prediction and phospho matrices."""

    overlap_count: int
    pred_mat_rows: int
    phospho_rows: int


@dataclass(frozen=True, slots=True)
class KinaseActivityInputs:
    """Trusted activity-stage inputs resolved by workflow validation."""

    pred_mat: pd.DataFrame
    phospho_matrix: pd.DataFrame
    threshold: float
    min_substrates: int
    top_n_substrates: int
    overlap_summary: PredMatOverlapSummary

    def __post_init__(self) -> None:
        try:
            pred_mat = KinasePredictionMatrix(
                frame=self.pred_mat,
                field_name="prediction_result.pred_mat",
                _assume_owned=True,
            ).frame
            phospho_matrix = ActivityMatrix(
                frame=self.phospho_matrix,
                field_name="dataset.phospho",
                _assume_owned=True,
            ).frame
        except PhosPyValidationError as exc:
            raise WorkflowBoundaryError(
                seam="kinase.activity.input_schema",
                next_action=(
                    "ensure prediction_result.pred_mat and dataset.phospho satisfy "
                    "activity-stage input table schema requirements"
                ),
                details={"schema_error": str(exc)},
                message_prefix="kinase workflow boundary validation failed",
            ) from exc
        if not isinstance(self.overlap_summary, PredMatOverlapSummary):
            raise WorkflowBoundaryError(
                "activity input overlap_summary must be PredMatOverlapSummary"
            )
        object.__setattr__(self, "pred_mat", pred_mat)
        object.__setattr__(self, "phospho_matrix", phospho_matrix)


@dataclass(frozen=True, slots=True, init=False)
class KinaseActivityResult:
    """Activity-stage outputs.

    Outputs are deliberately table-first and mirror the historical baseline stage:

    - ``weighted_activity``: weighted kinase activity matrix
    - ``thresholded_substrate_mean_activity``: sample-by-kinase mean phospho
      signal over predicted substrates above threshold
    - ``thresholded_substrate_counts``: number of selected thresholded predicted
      substrates per kinase
    - ``target_counts``: thresholded predicted target counts per kinase
    - ``target_table``: thresholded kinase-target edge table
    - ``activity_method``: stable method identity metadata for these outputs
    """

    activity_method: ActivityMethodMetadata
    _weighted_activity: pd.DataFrame = field(init=False, repr=False)
    _thresholded_substrate_mean_activity: pd.DataFrame = field(init=False, repr=False)
    _thresholded_substrate_counts: pd.Series = field(init=False, repr=False)
    _target_counts: pd.Series = field(init=False, repr=False)
    _target_table: pd.DataFrame = field(init=False, repr=False)
    _statistics_table: pd.DataFrame | None = field(init=False, repr=False)
    method_summary: ActivityMethodSummary | None

    def __init__(
        self,
        weighted_activity: pd.DataFrame,
        thresholded_substrate_mean_activity: pd.DataFrame,
        thresholded_substrate_counts: pd.Series,
        target_counts: pd.Series,
        target_table: pd.DataFrame,
        statistics_table: pd.DataFrame | None = None,
        method_summary: ActivityMethodSummary | None = None,
        activity_method: ActivityMethodMetadata = (
            SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD
        ),
        _assume_owned: bool = False,
    ) -> None:
        if not isinstance(activity_method, ActivityMethodMetadata):
            raise WorkflowBoundaryError(
                "activity_result.activity_method must be ActivityMethodMetadata"
            )
        weighted_activity = ActivityMatrix(
            frame=weighted_activity,
            field_name="activity_result.weighted_activity",
            _assume_owned=_assume_owned,
        ).frame
        thresholded_substrate_mean_activity = ActivityMatrix(
            frame=thresholded_substrate_mean_activity,
            field_name="activity_result.thresholded_substrate_mean_activity",
            _assume_owned=_assume_owned,
        ).frame
        thresholded_substrate_counts = ActivityCountSeries(
            series=thresholded_substrate_counts,
            field_name="activity_result.thresholded_substrate_counts",
            _assume_owned=_assume_owned,
        ).series
        target_counts = ActivityCountSeries(
            series=target_counts,
            field_name="activity_result.target_counts",
            _assume_owned=_assume_owned,
        ).series
        target_table = ActivityTargetTable(
            frame=target_table,
            _assume_owned=_assume_owned,
        ).frame
        if statistics_table is not None:
            statistics_table = ActivityStatisticsTable(
                frame=statistics_table,
                _assume_owned=_assume_owned,
            ).frame
        if method_summary is not None and not isinstance(
            method_summary, ActivityMethodSummary
        ):
            raise WorkflowBoundaryError(
                "activity_result.method_summary must be ActivityMethodSummary or None"
            )
        object.__setattr__(self, "_weighted_activity", weighted_activity)
        object.__setattr__(
            self,
            "_thresholded_substrate_mean_activity",
            thresholded_substrate_mean_activity,
        )
        object.__setattr__(
            self,
            "_thresholded_substrate_counts",
            thresholded_substrate_counts,
        )
        object.__setattr__(self, "_target_counts", target_counts)
        object.__setattr__(self, "activity_method", activity_method)
        object.__setattr__(self, "_target_table", target_table)
        object.__setattr__(self, "_statistics_table", statistics_table)
        object.__setattr__(self, "method_summary", method_summary)

    @property
    def weighted_activity(self) -> pd.DataFrame:
        return export_dataframe(self._weighted_activity)

    @property
    def thresholded_substrate_mean_activity(self) -> pd.DataFrame:
        return export_dataframe(self._thresholded_substrate_mean_activity)

    @property
    def thresholded_substrate_counts(self) -> pd.Series:
        return export_series(self._thresholded_substrate_counts)

    @property
    def target_counts(self) -> pd.Series:
        return export_series(self._target_counts)

    @property
    def target_table(self) -> pd.DataFrame:
        return export_dataframe(self._target_table)

    @property
    def statistics_table(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._statistics_table)

    @classmethod
    def _from_owned(
        cls,
        *,
        weighted_activity: pd.DataFrame,
        thresholded_substrate_mean_activity: pd.DataFrame,
        thresholded_substrate_counts: pd.Series,
        target_counts: pd.Series,
        target_table: pd.DataFrame,
        statistics_table: pd.DataFrame | None = None,
        method_summary: ActivityMethodSummary | None = None,
        activity_method: ActivityMethodMetadata = (
            SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD
        ),
    ) -> KinaseActivityResult:
        return cls(
            weighted_activity=weighted_activity,
            thresholded_substrate_mean_activity=thresholded_substrate_mean_activity,
            thresholded_substrate_counts=thresholded_substrate_counts,
            target_counts=target_counts,
            target_table=target_table,
            statistics_table=statistics_table,
            method_summary=method_summary,
            activity_method=activity_method,
            _assume_owned=True,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Return a weighted-activity snapshot isolated from this result."""

        return export_dataframe(self._weighted_activity)

    def thresholded_substrate_mean_activity_dataframe(self) -> pd.DataFrame:
        """Return a thresholded-mean snapshot isolated from this result."""

        return export_dataframe(self._thresholded_substrate_mean_activity)

    def target_table_dataframe(self) -> pd.DataFrame:
        """Return a target-table snapshot isolated from this result."""

        return export_dataframe(self._target_table)

    def statistics_table_dataframe(self) -> pd.DataFrame | None:
        """Return a statistics-table snapshot isolated from this result."""

        return export_optional_dataframe(self._statistics_table)
