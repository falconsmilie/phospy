from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd
from pydantic import Field, ValidationError, field_validator

from ..prediction.models import KinasePredictionResult, PredMatResult
from ..scoring import KinaseScoringResult
from ._models import PhospyRequestModel
from .errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    RequestValidationError,
)
from .tables import PredictionScoreMatrixSchema, PredMatSchema, SiteMatrixSchema

__all__ = [
    "SignalomeRequest",
    "ValidatedSignalomeRequest",
    "validate_signalome_request",
]


class SignalomeRequest(PhospyRequestModel):
    """Raw boundary request for public signalome construction."""

    kinases_of_interest: tuple[str, ...]
    kinase_network_threshold: float = Field(default=0.9, ge=0.0, le=1.0)
    signalome_cutoff: float = Field(default=0.5, ge=0.0, le=1.0)
    module_count: int | None = Field(default=None, ge=1)
    min_kinase_module_share_percent: float = Field(default=1.0, ge=0.0)

    @field_validator("kinases_of_interest", mode="before")
    @classmethod
    def normalize_kinases_of_interest(
        cls,
        value: Sequence[str],
    ) -> tuple[str, ...]:
        if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
            msg = "kinases_of_interest must be provided as a sequence of kinase names"
            raise ValueError(msg)

        normalized = tuple(dict.fromkeys(str(kinase) for kinase in value))
        if not normalized:
            msg = "kinases_of_interest must contain at least one kinase name"
            raise ValueError(msg)
        return normalized

    @classmethod
    def validate_request(cls, **data: object) -> SignalomeRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid signalome request",
                error=error,
            ) from error


@dataclass(slots=True)
class ValidatedSignalomeRequest:
    """Trusted aligned inputs for signalome construction."""

    request: SignalomeRequest
    scoring_matrix: pd.DataFrame
    pred_mat: pd.DataFrame
    expression_matrix: pd.DataFrame


def validate_signalome_request(
    *,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult | PredMatResult,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    kinase_network_threshold: float = 0.9,
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
) -> ValidatedSignalomeRequest:
    """Validate raw signalome inputs and return a trusted aligned bundle."""

    request = SignalomeRequest.validate_request(
        kinases_of_interest=kinases_of_interest,
        kinase_network_threshold=kinase_network_threshold,
        signalome_cutoff=signalome_cutoff,
        module_count=module_count,
        min_kinase_module_share_percent=min_kinase_module_share_percent,
    )

    scoring_matrix = PredictionScoreMatrixSchema.validate(
        _resolve_scoring_matrix(scoring_result),
        context="scoring_result",
    )
    pred_mat = _validate_prediction_result_pred_mat(prediction_result)
    validated_expression_matrix = SiteMatrixSchema.validate(
        expression_matrix,
        context="expression_matrix",
    )

    common_sites = [
        site_id
        for site_id in scoring_matrix.index.astype(str)
        if site_id in pred_mat.index and site_id in validated_expression_matrix.index
    ]
    if not common_sites:
        msg = (
            "scoring_result, prediction_result, and expression_matrix must share "
            "at least one phosphosite row"
        )
        raise InputCompatibilityError(msg)

    common_kinases = [
        kinase
        for kinase in scoring_matrix.columns.astype(str)
        if kinase in pred_mat.columns
    ]
    if not common_kinases:
        msg = (
            "scoring_result and prediction_result must share at least one kinase column"
        )
        raise InputCompatibilityError(msg)

    missing_koi = [
        kinase for kinase in request.kinases_of_interest if kinase not in common_kinases
    ]
    if missing_koi:
        missing = ", ".join(missing_koi)
        msg = f"kinases_of_interest are not available in the aligned signalome inputs: {missing}"
        raise InputCompatibilityError(msg)

    if request.module_count is not None and request.module_count > len(common_sites):
        msg = "module_count cannot exceed the number of aligned phosphosite rows"
        raise InputCompatibilityError(msg)

    return ValidatedSignalomeRequest(
        request=request,
        scoring_matrix=scoring_matrix.loc[common_sites, common_kinases],
        pred_mat=pred_mat.loc[common_sites, common_kinases],
        expression_matrix=validated_expression_matrix.loc[common_sites],
    )


def _resolve_scoring_matrix(scoring_result: KinaseScoringResult) -> pd.DataFrame:
    if not isinstance(scoring_result, KinaseScoringResult):
        msg = "scoring_result must be a KinaseScoringResult"
        raise RequestValidationError(msg)

    if scoring_result.combined_scores is not None:
        return scoring_result.combined_scores
    return scoring_result.profile_scores


def _validate_prediction_result_pred_mat(
    prediction_result: KinasePredictionResult | PredMatResult,
) -> pd.DataFrame:
    pred_mat = _resolve_pred_mat(prediction_result)
    if pred_mat.shape[1] == 0:
        msg = (
            "prediction_result does not contain any kinase columns because no "
            "candidate kinases qualified for prediction. Regenerate predMat with "
            "less restrictive top, score_threshold, or inclusion settings."
        )
        raise NoCandidateKinasesError(msg)
    return PredMatSchema.validate(pred_mat, context="prediction_result")


def _resolve_pred_mat(
    prediction_result: KinasePredictionResult | PredMatResult,
) -> pd.DataFrame:
    if isinstance(prediction_result, KinasePredictionResult):
        return prediction_result.pred_mat_result.to_frame(copy=False)
    if isinstance(prediction_result, PredMatResult):
        return prediction_result.to_frame(copy=False)
    msg = "prediction_result must be a KinasePredictionResult or PredMatResult"
    raise RequestValidationError(msg)
