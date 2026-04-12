from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd
from pydantic import Field, ValidationError

from ...errors import NoCandidateKinasesError, RequestValidationError
from ..compatibility import validate_pred_mat_overlap
from ..schema.tables import ActivitySiteMatrixSchema, PredMatSchema
from .shared import PhospyRequestModel, normalize_pred_mat_input

if TYPE_CHECKING:
    from ...prediction.results import PredMatResult


class KinaseActivityRequest(PhospyRequestModel):
    """Raw boundary options for downstream kinase activity analysis."""

    threshold: float = Field(default=0.6, ge=0.0, le=1.0)
    min_substrates: int = Field(default=3, ge=1)
    top_n_substrates: int = Field(default=20, ge=1)

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
        min_overlap: int = 1,
        min_fraction: float = 0.1,
    ) -> AnalysisInputs:
        """Build trusted analysis inputs from already-owned validated matrices."""

        validate_pred_mat_overlap(
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
        )


def validate_analysis_request(
    *,
    pred_mat: pd.DataFrame | PredMatResult,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> AnalysisInputs:
    """Validate raw analysis inputs and return trusted analysis inputs."""

    request = KinaseActivityRequest.validate_request(
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )
    normalized_pred_mat = normalize_pred_mat_input(pred_mat)
    if normalized_pred_mat is None:
        msg = f"{pred_context} must be provided"
        raise RequestValidationError(msg)
    if normalized_pred_mat.shape[1] == 0:
        msg = (
            f"{pred_context} does not contain any kinase columns because no "
            "candidate kinases qualified for prediction. Regenerate predMat with "
            "less restrictive top, score_threshold, or inclusion settings."
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
