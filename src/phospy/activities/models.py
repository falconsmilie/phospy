"""Activity stage result models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import export_dataframe
from phospy.errors.validation import PhosPyValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.tables.activity import (
    ActivityCountSeries,
    ActivityMatrix,
    ActivityTargetTable,
)
from phospy.tables.kinase import KinasePredictionMatrix


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


@dataclass(frozen=True, slots=True)
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
    """

    weighted_activity: pd.DataFrame
    thresholded_substrate_mean_activity: pd.DataFrame
    thresholded_substrate_counts: pd.Series
    target_counts: pd.Series
    target_table: pd.DataFrame
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        weighted_activity = ActivityMatrix(
            frame=self.weighted_activity,
            field_name="activity_result.weighted_activity",
            _assume_owned=_assume_owned,
        ).frame
        thresholded_substrate_mean_activity = ActivityMatrix(
            frame=self.thresholded_substrate_mean_activity,
            field_name="activity_result.thresholded_substrate_mean_activity",
            _assume_owned=_assume_owned,
        ).frame
        thresholded_substrate_counts = ActivityCountSeries(
            series=self.thresholded_substrate_counts,
            field_name="activity_result.thresholded_substrate_counts",
            _assume_owned=_assume_owned,
        ).series
        target_counts = ActivityCountSeries(
            series=self.target_counts,
            field_name="activity_result.target_counts",
            _assume_owned=_assume_owned,
        ).series
        target_table = ActivityTargetTable(
            frame=self.target_table,
            _assume_owned=_assume_owned,
        ).frame
        object.__setattr__(self, "weighted_activity", weighted_activity)
        object.__setattr__(
            self,
            "thresholded_substrate_mean_activity",
            thresholded_substrate_mean_activity,
        )
        object.__setattr__(
            self,
            "thresholded_substrate_counts",
            thresholded_substrate_counts,
        )
        object.__setattr__(self, "target_counts", target_counts)
        object.__setattr__(self, "target_table", target_table)

    @classmethod
    def _from_owned(
        cls,
        *,
        weighted_activity: pd.DataFrame,
        thresholded_substrate_mean_activity: pd.DataFrame,
        thresholded_substrate_counts: pd.Series,
        target_counts: pd.Series,
        target_table: pd.DataFrame,
    ) -> KinaseActivityResult:
        return cls(
            weighted_activity=weighted_activity,
            thresholded_substrate_mean_activity=thresholded_substrate_mean_activity,
            thresholded_substrate_counts=thresholded_substrate_counts,
            target_counts=target_counts,
            target_table=target_table,
            _assume_owned=True,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Return a weighted-activity snapshot isolated from this result."""

        return export_dataframe(self.weighted_activity)

    def thresholded_substrate_mean_activity_dataframe(self) -> pd.DataFrame:
        """Return a thresholded-mean snapshot isolated from this result."""

        return export_dataframe(self.thresholded_substrate_mean_activity)

    def target_table_dataframe(self) -> pd.DataFrame:
        """Return a target-table snapshot isolated from this result."""

        return export_dataframe(self.target_table)
