from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, ValidationError

from ...errors import (
    NoCandidateKinasesError,
    RequestValidationError,
    format_empty_prediction_matrix_message,
)
from ...internal.defaults import (
    DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES,
    DEFAULT_KINASE_ACTIVITY_THRESHOLD,
    DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES,
)
from ..compatibility import (
    DEFAULT_MIN_PRED_MAT_OVERLAP,
    DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION,
    PredMatOverlapSummary,
    validate_pred_mat_overlap,
)
from ..schema.tables import ActivitySiteMatrixSchema, PredMatSchema

if TYPE_CHECKING:
    from ...prediction.results import PredMatResult


class KinaseActivityRequest(BaseModel):
    """Raw boundary options for downstream kinase activity analysis."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    threshold: float = Field(default=DEFAULT_KINASE_ACTIVITY_THRESHOLD, ge=0.0, le=1.0)
    min_substrates: int = Field(default=DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES, ge=1)
    top_n_substrates: int = Field(
        default=DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES, ge=1
    )

    @classmethod
    def validate_request(cls, **data: object) -> KinaseActivityRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid kinase activity request",
                error=error,
            ) from error


@dataclass(slots=True)
class AnalysisInputs:
    """Trusted analysis inputs owned by the activity-analysis boundary."""

    pred_mat: pd.DataFrame
    phospho_matrix: pd.DataFrame
    threshold: float
    min_substrates: int
    top_n_substrates: int
    overlap_summary: PredMatOverlapSummary

    @classmethod
    def from_trusted_inputs(
        cls,
        *,
        pred_mat: pd.DataFrame,
        phospho_matrix: pd.DataFrame,
        threshold: float,
        min_substrates: int,
        top_n_substrates: int,
        pred_context: str = "pred_mat",
        matrix_context: str = "phospho_matrix",
        min_overlap: int = DEFAULT_MIN_PRED_MAT_OVERLAP,
        min_fraction: float = DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION,
    ) -> AnalysisInputs:
        """Build trusted analysis inputs from already-owned validated matrices."""

        overlap_summary = validate_pred_mat_overlap(
            pred_mat,
            phospho_matrix,
            pred_context=pred_context,
            matrix_context=matrix_context,
            min_overlap=min_overlap,
            min_fraction=min_fraction,
        )
        return cls(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
            overlap_summary=overlap_summary,
        )


def validate_analysis_request(
    *,
    pred_mat: pd.DataFrame | PredMatResult,
    phospho_matrix: pd.DataFrame,
    threshold: float = DEFAULT_KINASE_ACTIVITY_THRESHOLD,
    min_substrates: int = DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES,
    top_n_substrates: int = DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = DEFAULT_MIN_PRED_MAT_OVERLAP,
    min_fraction: float = DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION,
) -> AnalysisInputs:
    """Validate raw analysis inputs and return trusted analysis inputs."""

    request = KinaseActivityRequest.validate_request(
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )
    from ...prediction.results import PredMatResult

    normalized_pred_mat = (
        pred_mat.to_owned_frame() if isinstance(pred_mat, PredMatResult) else pred_mat
    )
    if normalized_pred_mat is None:
        msg = f"{pred_context} must be provided"
        raise RequestValidationError(msg)
    if normalized_pred_mat.shape[1] == 0:
        msg = format_empty_prediction_matrix_message(
            context=pred_context,
            phosphosite_rows=int(normalized_pred_mat.shape[0]),
            source_hint="analysis input",
        )
        raise NoCandidateKinasesError(msg)
    validated_pred_mat = PredMatSchema.validate(
        normalized_pred_mat,
        context=pred_context,
    )
    validated_matrix = ActivitySiteMatrixSchema.validate(
        phospho_matrix,
        context=matrix_context,
    )
    return AnalysisInputs.from_trusted_inputs(
        pred_mat=validated_pred_mat,
        phospho_matrix=validated_matrix,
        threshold=request.threshold,
        min_substrates=request.min_substrates,
        top_n_substrates=request.top_n_substrates,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )
