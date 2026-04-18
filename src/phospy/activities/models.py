"""Activity stage result models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


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
        object.__setattr__(self, "pred_mat", self.pred_mat.copy(deep=True))
        object.__setattr__(self, "phospho_matrix", self.phospho_matrix.copy(deep=True))


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

    def __post_init__(self) -> None:
        if not isinstance(self.weighted_activity, pd.DataFrame):
            raise TypeError(
                "activity_result.weighted_activity must be a pandas DataFrame"
            )
        if not isinstance(self.ksea_scores, pd.DataFrame):
            raise TypeError("activity_result.ksea_scores must be a pandas DataFrame")
        if not isinstance(self.ksea_counts, pd.Series):
            raise TypeError("activity_result.ksea_counts must be a pandas Series")
        if not isinstance(self.target_counts, pd.Series):
            raise TypeError("activity_result.target_counts must be a pandas Series")
        if not isinstance(self.target_table, pd.DataFrame):
            raise TypeError("activity_result.target_table must be a pandas DataFrame")
        object.__setattr__(
            self,
            "weighted_activity",
            self.weighted_activity.copy(deep=True),
        )
        object.__setattr__(
            self,
            "ksea_scores",
            self.ksea_scores.copy(deep=True),
        )
        object.__setattr__(
            self,
            "ksea_counts",
            self.ksea_counts.copy(deep=True),
        )
        object.__setattr__(
            self,
            "target_counts",
            self.target_counts.copy(deep=True),
        )
        object.__setattr__(
            self,
            "target_table",
            self.target_table.copy(deep=True),
        )
