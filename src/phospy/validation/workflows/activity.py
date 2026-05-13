"""Validation helpers for kinase activity-stage execution inputs."""

from __future__ import annotations

import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.activities.models import (
    KinaseActivityInputs,
    PredMatOverlapSummary,
)
from phospy.validation.common.dataframes import (
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.common.missing_values import MissingValuePolicy
from phospy.validation.common.numbers import require_int_at_least, require_real_between
from phospy.validation.common.numeric_frames import require_numeric_matrix

DEFAULT_MIN_PRED_MAT_OVERLAP = 1
DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION = 0.5


class KinaseActivityInputValidator:
    """Validate and normalize trusted kinase activity stage inputs."""

    def run(
        self,
        *,
        pred_mat: object,
        phospho_matrix: object,
        threshold: object,
        min_substrates: object,
        top_n_substrates: object,
        min_overlap: int = DEFAULT_MIN_PRED_MAT_OVERLAP,
        min_fraction: float = DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION,
    ) -> KinaseActivityInputs:
        normalized_threshold = require_real_between(
            threshold,
            field_name="activity_config.threshold",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowBoundaryError,
        )
        normalized_min_substrates = require_int_at_least(
            min_substrates,
            field_name="activity_config.min_substrates",
            minimum=1,
            error_type=WorkflowBoundaryError,
        )
        normalized_top_n_substrates = require_int_at_least(
            top_n_substrates,
            field_name="activity_config.top_n_substrates",
            minimum=1,
            error_type=WorkflowBoundaryError,
        )
        require_int_at_least(
            min_overlap,
            field_name="activity stage min_overlap",
            minimum=1,
            error_type=WorkflowBoundaryError,
        )
        normalized_min_fraction = require_real_between(
            min_fraction,
            field_name="activity stage min_fraction",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowBoundaryError,
        )

        validated_pred_mat = require_numeric_matrix(
            pred_mat,
            field_name="prediction_result.pred_mat",
            allow_empty=False,
            missing_value_policy=MissingValuePolicy.ALLOW,
            error_type=WorkflowBoundaryError,
        )
        validated_pred_mat = require_unique_index(
            validated_pred_mat,
            field_name="prediction_result.pred_mat",
            error_type=WorkflowBoundaryError,
        )
        validated_pred_mat = require_unique_columns(
            validated_pred_mat,
            field_name="prediction_result.pred_mat",
            error_type=WorkflowBoundaryError,
        )

        validated_phospho = require_numeric_matrix(
            phospho_matrix,
            field_name="dataset.phospho",
            allow_empty=False,
            missing_value_policy=MissingValuePolicy.ALLOW,
            error_type=WorkflowBoundaryError,
        )
        validated_phospho = require_unique_index(
            validated_phospho,
            field_name="dataset.phospho",
            error_type=WorkflowBoundaryError,
        )
        validated_phospho = require_unique_columns(
            validated_phospho,
            field_name="dataset.phospho",
            error_type=WorkflowBoundaryError,
        )

        overlap_summary = self._validate_overlap(
            pred_mat=validated_pred_mat,
            phospho_matrix=validated_phospho,
            min_overlap=min_overlap,
            min_fraction=normalized_min_fraction,
        )
        return KinaseActivityInputs(
            pred_mat=validated_pred_mat,
            phospho_matrix=validated_phospho,
            threshold=normalized_threshold,
            min_substrates=normalized_min_substrates,
            top_n_substrates=normalized_top_n_substrates,
            overlap_summary=overlap_summary,
        )

    def _validate_overlap(
        self,
        *,
        pred_mat: pd.DataFrame,
        phospho_matrix: pd.DataFrame,
        min_overlap: int,
        min_fraction: float,
    ) -> PredMatOverlapSummary:
        overlap_count = int(pred_mat.index.intersection(phospho_matrix.index).size)
        pred_mat_rows = int(pred_mat.index.size)
        phospho_rows = int(phospho_matrix.index.size)
        overlap_fraction = overlap_count / max(phospho_rows, 1)
        if overlap_count == 0:
            self._raise_overlap_error(
                overlap_count=overlap_count,
                pred_mat_rows=pred_mat_rows,
                phospho_rows=phospho_rows,
                min_overlap=min_overlap,
                min_fraction=min_fraction,
            )
        if overlap_count < min_overlap or overlap_fraction < min_fraction:
            self._raise_overlap_error(
                overlap_count=overlap_count,
                pred_mat_rows=pred_mat_rows,
                phospho_rows=phospho_rows,
                min_overlap=min_overlap,
                min_fraction=min_fraction,
            )
        return PredMatOverlapSummary(
            overlap_count=overlap_count,
            pred_mat_rows=pred_mat_rows,
            phospho_rows=phospho_rows,
        )

    @staticmethod
    def _raise_overlap_error(
        *,
        overlap_count: int,
        pred_mat_rows: int,
        phospho_rows: int,
        min_overlap: int,
        min_fraction: float,
    ) -> None:
        raise WorkflowBoundaryError(
            "kinase workflow boundary validation failed at "
            "seam=kinase.activity.input_overlap; "
            f"overlap_sites={overlap_count}, pred_mat_sites={pred_mat_rows}, "
            f"phospho_sites={phospho_rows}, min_overlap={min_overlap}, "
            f"min_fraction={min_fraction}; next_action=ensure "
            "prediction_result.pred_mat and dataset.phospho share phosphosite "
            "IDs and originate from the same workflow run"
        )


__all__ = [
    "DEFAULT_MIN_PRED_MAT_OVERLAP",
    "DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION",
    "KinaseActivityInputValidator",
]
