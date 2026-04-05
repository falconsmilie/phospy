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
from .validation.analysis import (
    KinaseActivityRequest,
    ValidatedAnalysisRequest,
    validate_analysis_request,
)
from .validation.tables import SiteMatrixSchema
from .writers import KinaseActivityResultWriter, KinaseActivityWriter

__all__ = ["KinaseActivityAnalyzer", "KinaseActivityResult"]


@dataclass(slots=True)
class KinaseActivityResult:
    """Detached snapshot bundle for downstream kinase-activity outputs.

    The tables and series stored here are produced result tables, not live views
    into analyzer request inputs or dataset workspace state. Mutating them only
    affects this result instance.
    """

    weighted_activity: pd.DataFrame
    ksea_scores: pd.DataFrame
    ksea_counts: pd.Series
    target_counts: pd.Series
    target_table: pd.DataFrame


@dataclass(slots=True)
class KinaseActivityAnalyzer:
    """Application service for downstream kinase activity analysis."""

    result_writer: type[KinaseActivityResultWriter] = KinaseActivityWriter

    @staticmethod
    def _load_pred_mat(pred_mat_path: str | Path) -> pd.DataFrame:
        """Load and validate a kinase prediction matrix from disk."""

        return load_pred_mat(pred_mat_path)

    def _validate_request(
        self,
        *,
        pred_mat: pd.DataFrame,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
    ) -> ValidatedAnalysisRequest:
        return validate_analysis_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )

    @staticmethod
    def analyze(
        pred_mat: pd.DataFrame,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
    ) -> KinaseActivityResult:
        """Compute downstream kinase summaries from raw public inputs."""

        analyzer = KinaseActivityAnalyzer()
        request = analyzer._validate_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )
        return analyzer._analyze_validated_request(request=request)

    @staticmethod
    def _analyze_validated_request(
        *, request: ValidatedAnalysisRequest
    ) -> KinaseActivityResult:
        if not isinstance(request, ValidatedAnalysisRequest):
            msg = (
                "_analyze_validated_request requires a ValidatedAnalysisRequest. "
                "Call _validate_request(...) first."
            )
            raise TypeError(msg)
        weighted_activity = compute_weighted_kinase_activity(
            pred_mat=request.pred_mat,
            phospho_matrix=request.phospho_matrix,
            top_n_substrates=request.request.top_n_substrates,
            min_substrates=request.request.min_substrates,
        )
        ksea_scores, ksea_counts = compute_ksea_scores(
            pred_mat=request.pred_mat,
            phospho_matrix=request.phospho_matrix,
            threshold=request.request.threshold,
            min_substrates=request.request.min_substrates,
        )
        target_counts = count_predicted_targets(
            request.pred_mat,
            threshold=request.request.threshold,
        )
        target_table = build_kinase_target_table(
            request.pred_mat,
            threshold=request.request.threshold,
        )

        return KinaseActivityResult(
            weighted_activity=weighted_activity,
            ksea_scores=ksea_scores,
            ksea_counts=ksea_counts,
            target_counts=target_counts,
            target_table=target_table,
        )

    def _load_and_analyze(
        self,
        pred_mat_path: str | Path,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
    ) -> KinaseActivityResult:
        """Load a prediction matrix from disk and compute downstream kinase summaries."""

        validated_pred_mat = self._load_pred_mat(pred_mat_path)
        validated_matrix = SiteMatrixSchema.validate(
            phospho_matrix,
            context="phospho_matrix",
        )
        raw_request = KinaseActivityRequest.validate_request(
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )
        request = ValidatedAnalysisRequest.from_trusted_inputs(
            request=raw_request,
            pred_mat=validated_pred_mat,
            phospho_matrix=validated_matrix,
            pred_context="pred_mat",
            matrix_context="phospho_matrix",
            min_overlap=1,
            min_fraction=0.1,
        )
        return self._analyze_validated_request(request=request)

    def write_outputs(self, result: KinaseActivityResult, outdir: str | Path) -> None:
        self.result_writer.write(result, outdir)
