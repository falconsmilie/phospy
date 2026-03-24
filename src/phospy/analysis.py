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
from .validation.requests import KinaseActivityRequest
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

    @classmethod
    def from_validated_pred_mat(
        cls,
        pred_mat: pd.DataFrame,
    ) -> KinaseActivityAnalyzer:
        instance = cls.__new__(cls)
        instance.pred_mat = pred_mat.copy()
        return instance

    def build_target_table(self, threshold: float = 0.6) -> pd.DataFrame:
        return build_kinase_target_table(self.pred_mat, threshold=threshold)

    def count_predicted_targets(self, threshold: float = 0.6) -> pd.Series:
        return count_predicted_targets(self.pred_mat, threshold=threshold)

    def compute_weighted_activity(
        self,
        phospho_matrix: pd.DataFrame,
        top_n_substrates: int = 20,
        min_substrates: int = 3,
    ) -> pd.DataFrame:
        return compute_weighted_kinase_activity(
            pred_mat=self.pred_mat,
            phospho_matrix=phospho_matrix,
            top_n_substrates=top_n_substrates,
            min_substrates=min_substrates,
        )

    def compute_ksea_scores(
        self,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
    ) -> tuple[pd.DataFrame, pd.Series]:
        return compute_ksea_scores(
            pred_mat=self.pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
        )

    def analyze(
        self,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
        request: KinaseActivityRequest | None = None,
    ) -> KinaseActivityResult:
        validated_matrix = SiteMatrixSchema.validate(
            phospho_matrix,
            context="phospho_matrix",
        )
        activity_request = request or KinaseActivityRequest.validate_request(
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )
        validate_pred_mat_overlap(self.pred_mat, validated_matrix)
        weighted_activity = self.compute_weighted_activity(
            phospho_matrix=validated_matrix,
            top_n_substrates=activity_request.top_n_substrates,
            min_substrates=activity_request.min_substrates,
        )
        ksea_scores, ksea_counts = self.compute_ksea_scores(
            phospho_matrix=validated_matrix,
            threshold=activity_request.threshold,
            min_substrates=activity_request.min_substrates,
        )
        target_counts = self.count_predicted_targets(
            threshold=activity_request.threshold
        )
        target_table = self.build_target_table(threshold=activity_request.threshold)
        return KinaseActivityResult(
            weighted_activity=weighted_activity,
            ksea_scores=ksea_scores,
            ksea_counts=ksea_counts,
            target_counts=target_counts,
            target_table=target_table,
        )

    @staticmethod
    def write_outputs(result: KinaseActivityResult, outdir: str | Path) -> None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        result.weighted_activity.to_csv(outdir / "kinase_activity_matrix.csv")
        result.ksea_scores.to_csv(outdir / "ksea_scores.csv")
        result.ksea_counts.to_csv(outdir / "ksea_counts.csv")
        result.target_counts.to_csv(outdir / "kinase_target_counts.csv")
        result.target_table.to_csv(outdir / "kinase_target_table.csv", index=False)
