from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd
from pydantic import Field, ValidationError, field_validator

from ..compatibility import validate_signalome_alignment
from ..domain import (
    resolve_scoring_matrix,
    validate_prediction_result_pred_mat,
    validate_signalome_site_grouping,
)
from ..errors import RequestValidationError
from ..values.collections import (
    normalize_site_to_protein_mapping,
    normalize_string_sequence,
)
from .shared import PhospyRequestModel

if TYPE_CHECKING:
    from ...prediction.results import KinasePredictionResult, PredMatResult
    from ...prediction.scoring import KinaseScoringResult


class SignalomeRequest(PhospyRequestModel):
    """Raw boundary request for public signalome construction."""

    kinases_of_interest: tuple[str, ...]
    site_to_protein: dict[str, str] | None = None
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
        return normalize_string_sequence(
            value,
            field_name="kinases_of_interest",
            empty_message="kinases_of_interest must contain at least one kinase name",
            invalid_message=(
                "kinases_of_interest must be provided as a sequence of kinase names"
            ),
            deduplicate=True,
        )

    @field_validator("site_to_protein", mode="before")
    @classmethod
    def normalize_site_to_protein(
        cls,
        value: object,
    ) -> dict[str, str] | None:
        return normalize_site_to_protein_mapping(value)

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
    site_to_protein: pd.Series


def validate_signalome_request(
    *,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult | PredMatResult,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    site_to_protein: Mapping[str, str] | None = None,
    kinase_network_threshold: float = 0.9,
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
) -> ValidatedSignalomeRequest:
    """Validate raw signalome inputs and return a trusted aligned bundle."""

    request = SignalomeRequest.validate_request(
        kinases_of_interest=kinases_of_interest,
        site_to_protein=site_to_protein,
        kinase_network_threshold=kinase_network_threshold,
        signalome_cutoff=signalome_cutoff,
        module_count=module_count,
        min_kinase_module_share_percent=min_kinase_module_share_percent,
    )

    return _build_validated_signalome_request(
        request=request,
        scoring_matrix=resolve_scoring_matrix(scoring_result),
        pred_mat=validate_prediction_result_pred_mat(prediction_result),
        expression_matrix=expression_matrix,
        scoring_context="scoring_result",
        pred_mat_context="prediction_result",
        expression_context="expression_matrix",
    )


def _build_validated_signalome_request(
    *,
    request: SignalomeRequest,
    scoring_matrix: pd.DataFrame,
    pred_mat: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    scoring_context: str,
    pred_mat_context: str,
    expression_context: str,
) -> ValidatedSignalomeRequest:
    (
        validated_scoring_matrix,
        validated_pred_mat,
        validated_expression_matrix,
        common_sites,
    ) = validate_signalome_alignment(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=request.kinases_of_interest,
        module_count=request.module_count,
        scoring_context=scoring_context,
        pred_mat_context=pred_mat_context,
        expression_context=expression_context,
    )

    validated_site_to_protein = validate_signalome_site_grouping(
        site_ids=common_sites,
        site_to_protein=request.site_to_protein,
    )

    return ValidatedSignalomeRequest(
        request=request,
        scoring_matrix=validated_scoring_matrix,
        pred_mat=validated_pred_mat,
        expression_matrix=validated_expression_matrix,
        site_to_protein=validated_site_to_protein,
    )
