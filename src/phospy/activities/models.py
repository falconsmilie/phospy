"""Activity stage result models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import own_dataframe, own_series


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
        if not isinstance(self.pred_mat, pd.DataFrame):
            raise TypeError("activity input pred_mat must be a pandas DataFrame")
        if not isinstance(self.phospho_matrix, pd.DataFrame):
            raise TypeError("activity input phospho_matrix must be a pandas DataFrame")
        if not isinstance(self.overlap_summary, PredMatOverlapSummary):
            raise TypeError(
                "activity input overlap_summary must be PredMatOverlapSummary"
            )


@dataclass(frozen=True, slots=True)
class KinaseActivityResult:
    """Activity-stage outputs.

    Outputs are deliberately table-first and mirror the legacy scientific stage:

    - ``weighted_activity``: weighted kinase activity matrix
    - ``ksea_scores``: KSEA-style sample-by-kinase score matrix
    - ``ksea_counts``: number of selected KSEA substrates per kinase
    - ``target_counts``: thresholded predicted target counts per kinase
    - ``target_table``: thresholded kinase-target edge table
    """

    weighted_activity: pd.DataFrame
    ksea_scores: pd.DataFrame
    ksea_counts: pd.Series
    target_counts: pd.Series
    target_table: pd.DataFrame
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        weighted_activity = own_dataframe(
            self.weighted_activity,
            field_name="activity_result.weighted_activity",
            assume_owned=_assume_owned,
        )
        ksea_scores = own_dataframe(
            self.ksea_scores,
            field_name="activity_result.ksea_scores",
            assume_owned=_assume_owned,
        )
        ksea_counts = own_series(
            self.ksea_counts,
            field_name="activity_result.ksea_counts",
            assume_owned=_assume_owned,
        )
        target_counts = own_series(
            self.target_counts,
            field_name="activity_result.target_counts",
            assume_owned=_assume_owned,
        )
        target_table = own_dataframe(
            self.target_table,
            field_name="activity_result.target_table",
            assume_owned=_assume_owned,
        )
        object.__setattr__(self, "weighted_activity", weighted_activity)
        object.__setattr__(self, "ksea_scores", ksea_scores)
        object.__setattr__(self, "ksea_counts", ksea_counts)
        object.__setattr__(self, "target_counts", target_counts)
        object.__setattr__(self, "target_table", target_table)

    @classmethod
    def _from_owned(
        cls,
        *,
        weighted_activity: pd.DataFrame,
        ksea_scores: pd.DataFrame,
        ksea_counts: pd.Series,
        target_counts: pd.Series,
        target_table: pd.DataFrame,
    ) -> KinaseActivityResult:
        return cls(
            weighted_activity=weighted_activity,
            ksea_scores=ksea_scores,
            ksea_counts=ksea_counts,
            target_counts=target_counts,
            target_table=target_table,
            _assume_owned=True,
        )
