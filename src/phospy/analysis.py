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
    ValidatedKinaseActivityInputs,
    build_kinase_activity_inputs,
    build_loaded_kinase_activity_inputs,
)
from .writers import KinaseActivityResultWriter, KinaseActivityWriter


@dataclass(slots=True)
class KinaseActivityResult:
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
    def load_pred_mat(pred_mat_path: str | Path) -> pd.DataFrame:
        """Load and validate a kinase prediction matrix from disk."""

        return load_pred_mat(pred_mat_path)

    @staticmethod
    def analyze(
        pred_mat: pd.DataFrame,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
    ) -> KinaseActivityResult:
        """Compute downstream kinase summaries from a validated prediction matrix."""

        request = KinaseActivityRequest.validate_request(
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )
        return KinaseActivityAnalyzer.analyze_request(
            request=request,
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
        )

    @staticmethod
    def analyze_request(
        *,
        request: KinaseActivityRequest,
        pred_mat: pd.DataFrame,
        phospho_matrix: pd.DataFrame,
    ) -> KinaseActivityResult:
        validated_inputs = build_kinase_activity_inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
        )
        return KinaseActivityAnalyzer._analyze_validated_request(
            request=request,
            validated_inputs=validated_inputs,
        )

    @staticmethod
    def _analyze_validated_request(
        *,
        request: KinaseActivityRequest,
        validated_inputs: ValidatedKinaseActivityInputs,
    ) -> KinaseActivityResult:
        weighted_activity = compute_weighted_kinase_activity(
            pred_mat=validated_inputs.pred_mat,
            phospho_matrix=validated_inputs.phospho_matrix,
            top_n_substrates=request.top_n_substrates,
            min_substrates=request.min_substrates,
        )
        ksea_scores, ksea_counts = compute_ksea_scores(
            pred_mat=validated_inputs.pred_mat,
            phospho_matrix=validated_inputs.phospho_matrix,
            threshold=request.threshold,
            min_substrates=request.min_substrates,
        )
        target_counts = count_predicted_targets(
            validated_inputs.pred_mat,
            threshold=request.threshold,
        )
        target_table = build_kinase_target_table(
            validated_inputs.pred_mat,
            threshold=request.threshold,
        )

        return KinaseActivityResult(
            weighted_activity=weighted_activity,
            ksea_scores=ksea_scores,
            ksea_counts=ksea_counts,
            target_counts=target_counts,
            target_table=target_table,
        )

    def load_and_analyze(
        self,
        pred_mat_path: str | Path,
        phospho_matrix: pd.DataFrame,
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
    ) -> KinaseActivityResult:
        """Load a prediction matrix from disk and compute downstream kinase summaries."""

        validated_pred_mat = self.load_pred_mat(pred_mat_path)
        request = KinaseActivityRequest.validate_request(
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )
        validated_inputs = build_loaded_kinase_activity_inputs(
            validated_pred_mat=validated_pred_mat,
            phospho_matrix=phospho_matrix,
            pred_context="pred_mat",
            matrix_context="phospho_matrix",
            min_overlap=1,
            min_fraction=0.1,
        )
        return self._analyze_validated_request(
            request=request,
            validated_inputs=validated_inputs,
        )

    def write_outputs(self, result: KinaseActivityResult, outdir: str | Path) -> None:
        self.result_writer.write(result, outdir)
