from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from pydantic import Field, ValidationError

from ._models import PhospyRequestModel
from .compatibility import validate_pred_mat_overlap
from .errors import RequestValidationError
from .tables import PredMatSchema, SiteMatrixSchema


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


@dataclass(frozen=True, slots=True)
class ValidatedAnalysisRequest:
    """Trusted boundary request for the public :class:`phospy.KinaseActivityAnalyzer` API."""

    request: KinaseActivityRequest
    pred_mat: pd.DataFrame
    phospho_matrix: pd.DataFrame


# Backward-compatible alias for older internal names.
ValidatedKinaseActivityInputs = ValidatedAnalysisRequest


def build_validated_analysis_request(
    *,
    request: KinaseActivityRequest,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> ValidatedAnalysisRequest:
    validated_pred_mat = PredMatSchema.validate(pred_mat, context=pred_context)
    return build_loaded_analysis_request(
        request=request,
        validated_pred_mat=validated_pred_mat,
        phospho_matrix=phospho_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )


def build_loaded_analysis_request(
    *,
    request: KinaseActivityRequest,
    validated_pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> ValidatedAnalysisRequest:
    validated_matrix = SiteMatrixSchema.validate(phospho_matrix, context=matrix_context)
    validate_pred_mat_overlap(
        validated_pred_mat,
        validated_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )
    return ValidatedAnalysisRequest(
        request=request,
        pred_mat=validated_pred_mat,
        phospho_matrix=validated_matrix,
    )


def validate_analysis_request(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> ValidatedAnalysisRequest:
    request = KinaseActivityRequest.validate_request(
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )
    return build_validated_analysis_request(
        request=request,
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )


def build_kinase_activity_inputs(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> ValidatedAnalysisRequest:
    """Compatibility wrapper returning validated activity inputs without options."""

    request = KinaseActivityRequest.validate_request()
    return build_validated_analysis_request(
        request=request,
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )


def build_loaded_kinase_activity_inputs(
    *,
    validated_pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> ValidatedAnalysisRequest:
    request = KinaseActivityRequest.validate_request()
    return build_loaded_analysis_request(
        request=request,
        validated_pred_mat=validated_pred_mat,
        phospho_matrix=phospho_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )


def validate_kinase_activity_inputs(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    validated_inputs = build_kinase_activity_inputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )
    return validated_inputs.pred_mat, validated_inputs.phospho_matrix


__all__ = [
    "KinaseActivityRequest",
    "ValidatedAnalysisRequest",
    "ValidatedKinaseActivityInputs",
    "build_kinase_activity_inputs",
    "build_loaded_analysis_request",
    "build_loaded_kinase_activity_inputs",
    "build_validated_analysis_request",
    "validate_analysis_request",
    "validate_kinase_activity_inputs",
]
