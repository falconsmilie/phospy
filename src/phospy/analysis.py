from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .activities import (
    build_kinase_target_table,
    compute_ksea_scores,
    compute_weighted_kinase_activity,
    count_predicted_targets,
)
from .io import load_pred_mat
from .validation.compatibility import validate_pred_mat_overlap
from .validation.errors import PhospyValidationError
from .validation.tables import PredMatSchema, SiteMatrixSchema


@dataclass(slots=True)
class KinaseActivityResult:
    weighted_activity: pd.DataFrame
    ksea_scores: pd.DataFrame
    ksea_counts: pd.Series
    target_counts: pd.Series
    target_table: pd.DataFrame


class KinaseActivityAnalyzer:
    """Owns downstream kinase summaries derived from a prediction matrix."""

    def __init__(self, pred_mat: pd.DataFrame) -> None:
        self.pred_mat = PredMatSchema.validate(pred_mat, context="pred_mat")

    @classmethod
    def from_csv(cls, pred_mat_path: str | Path) -> KinaseActivityAnalyzer:
        pred_mat = load_pred_mat(pred_mat_path)
        return cls(pred_mat=pred_mat)

    def build_target_table(self, threshold: float = 0.6) -> pd.DataFrame:
        _validate_probability_threshold(threshold, name="threshold")
        return build_kinase_target_table(self.pred_mat, threshold=threshold)

    def count_predicted_targets(self, threshold: float = 0.6) -> pd.Series:
        _validate_probability_threshold(threshold, name="threshold")
        return count_predicted_targets(self.pred_mat, threshold=threshold)

    def compute_weighted_activity(
        self,
        phospho_matrix: pd.DataFrame,
        top_n_substrates: int = 20,
        min_substrates: int = 3,
    ) -> pd.DataFrame:
        _validate_positive_int(top_n_substrates, name="top_n_substrates")
        _validate_positive_int(min_substrates, name="min_substrates")
        validated_matrix = _validate_site_matrix(phospho_matrix)
        return compute_weighted_kinase_activity(
            pred_mat=self.pred_mat,
            phospho_matrix=validated_matrix,
            top_n_substrates=top_n_substrates,
            min_substrates=min_substrates,
        )

    def compute_ksea_scores(
        self,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
    ) -> tuple[pd.DataFrame, pd.Series]:
        _validate_probability_threshold(threshold, name="threshold")
        _validate_positive_int(min_substrates, name="min_substrates")
        validated_matrix = _validate_site_matrix(phospho_matrix)
        return compute_ksea_scores(
            pred_mat=self.pred_mat,
            phospho_matrix=validated_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
        )

    def analyze(
        self,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
    ) -> KinaseActivityResult:
        validated_matrix = _validate_site_matrix(phospho_matrix)
        validate_pred_mat_overlap(self.pred_mat, validated_matrix)
        weighted_activity = self.compute_weighted_activity(
            phospho_matrix=validated_matrix,
            top_n_substrates=top_n_substrates,
            min_substrates=min_substrates,
        )
        ksea_scores, ksea_counts = self.compute_ksea_scores(
            phospho_matrix=validated_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
        )
        target_counts = self.count_predicted_targets(threshold=threshold)
        target_table = self.build_target_table(threshold=threshold)
        return KinaseActivityResult(
            weighted_activity=weighted_activity,
            ksea_scores=ksea_scores,
            ksea_counts=ksea_counts,
            target_counts=target_counts,
            target_table=target_table,
        )

    @staticmethod
    def write_outputs(result: KinaseActivityResult, outdir: str | Path) -> None:
        from .writers import KinaseActivityWriter

        KinaseActivityWriter.write(result, outdir)


def _validate_probability_threshold(value: float, *, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhospyValidationError(f"{name} must be a number between 0.0 and 1.0")
    if not 0.0 <= float(value) <= 1.0:
        raise PhospyValidationError(f"{name} must be between 0.0 and 1.0")


def _validate_positive_int(value: int, *, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhospyValidationError(f"{name} must be an integer")
    if value < 1:
        raise PhospyValidationError(f"{name} must be at least 1")


def _validate_site_matrix(phospho_matrix: pd.DataFrame) -> pd.DataFrame:
    return SiteMatrixSchema.validate(phospho_matrix, context="phospho_matrix")
