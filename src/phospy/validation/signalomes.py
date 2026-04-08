from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
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
        if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
            msg = "kinases_of_interest must be provided as a sequence of kinase names"
            raise ValueError(msg)

        normalized = tuple(dict.fromkeys(str(kinase) for kinase in value))
        if not normalized:
            msg = "kinases_of_interest must contain at least one kinase name"
            raise ValueError(msg)
        return normalized

    @field_validator("site_to_protein", mode="before")
    @classmethod
    def normalize_site_to_protein(
        cls,
        value: object,
    ) -> dict[str, str] | None:
        if value is None:
            return None
        if not isinstance(value, dict):
            try:
                value = dict(value)  # type: ignore[arg-type]
            except (TypeError, ValueError) as error:
                msg = "site_to_protein must be provided as a mapping of site IDs to protein IDs"
                raise ValueError(msg) from error

        normalized: dict[str, str] = {}
        for raw_site_id, raw_protein_id in value.items():
            site_id = str(raw_site_id)
            protein_id = str(raw_protein_id).strip()
            if not site_id:
                msg = "site_to_protein keys must be non-empty site IDs"
                raise ValueError(msg)
            if not protein_id:
                msg = "site_to_protein values must be non-empty protein IDs"
                raise ValueError(msg)
            normalized[site_id] = protein_id
        if not normalized:
            msg = "site_to_protein must contain at least one site-to-protein mapping"
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
        scoring_matrix=_resolve_scoring_matrix(scoring_result),
        pred_mat=_validate_prediction_result_pred_mat(prediction_result),
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
    validated_scoring_matrix = PredictionScoreMatrixSchema.validate(
        scoring_matrix,
        context=scoring_context,
    )
    validated_pred_mat = PredMatSchema.validate(
        pred_mat,
        context=pred_mat_context,
    )
    _ensure_finite_pred_mat(
        validated_pred_mat,
        context=pred_mat_context,
    )
    validated_expression_matrix = SiteMatrixSchema.validate(
        expression_matrix,
        context=expression_context,
    )

    common_sites = [
        site_id
        for site_id in validated_scoring_matrix.index.astype(str)
        if site_id in validated_pred_mat.index
        and site_id in validated_expression_matrix.index
    ]
    if not common_sites:
        msg = (
            f"{scoring_context}, {pred_mat_context}, and {expression_context} "
            "must share at least one phosphosite row"
        )
        raise InputCompatibilityError(msg)

    common_kinases = [
        kinase
        for kinase in validated_scoring_matrix.columns.astype(str)
        if kinase in validated_pred_mat.columns
    ]
    if not common_kinases:
        msg = (
            f"{scoring_context} and {pred_mat_context} must share at least one "
            "kinase column"
        )
        raise InputCompatibilityError(msg)

    missing_koi = [
        kinase for kinase in request.kinases_of_interest if kinase not in common_kinases
    ]
    if missing_koi:
        missing = ", ".join(missing_koi)
        msg = (
            "kinases_of_interest are not available in the aligned signalome "
            f"inputs: {missing}"
        )
        raise InputCompatibilityError(msg)

    if request.module_count is not None and request.module_count > len(common_sites):
        msg = "module_count cannot exceed the number of aligned phosphosite rows"
        raise InputCompatibilityError(msg)

    validated_site_to_protein = _validate_signalome_site_grouping(
        site_ids=common_sites,
        site_to_protein=request.site_to_protein,
    )

    return ValidatedSignalomeRequest(
        request=request,
        scoring_matrix=validated_scoring_matrix.loc[common_sites, common_kinases],
        pred_mat=validated_pred_mat.loc[common_sites, common_kinases],
        expression_matrix=validated_expression_matrix.loc[common_sites],
        site_to_protein=validated_site_to_protein,
    )


def _ensure_finite_pred_mat(
    frame: pd.DataFrame,
    *,
    context: str,
) -> None:
    failures: list[str] = []
    for column in frame.columns.astype(str):
        series = frame.loc[:, column]
        invalid_mask = ~np.isfinite(series.to_numpy(dtype=float))
        if invalid_mask.any():
            sample_values = series.loc[invalid_mask].astype(str).unique()[:3]
            sample_preview = ", ".join(str(value) for value in sample_values)
            failures.append(f"{column} ({sample_preview})")
    if failures:
        failures_str = "; ".join(failures)
        msg = f"{context} contains non-finite values in numeric columns: {failures_str}"
        raise InputCompatibilityError(msg)


def _validate_signalome_site_grouping(
    *,
    site_ids: Sequence[str],
    site_to_protein: Mapping[str, str] | None,
) -> pd.Series:
    if site_to_protein is not None:
        return _validate_explicit_site_to_protein_mapping(
            site_ids=site_ids,
            site_to_protein=site_to_protein,
        )
    return _validate_supported_signalome_site_ids(site_ids)


def _validate_explicit_site_to_protein_mapping(
    *,
    site_ids: Sequence[str],
    site_to_protein: Mapping[str, str],
) -> pd.Series:
    missing_site_ids = [
        site_id for site_id in site_ids if site_id not in site_to_protein
    ]
    if missing_site_ids:
        preview = ", ".join(missing_site_ids[:3])
        msg = (
            "site_to_protein must define a protein ID for every aligned phosphosite "
            f"row. Missing mappings for: {preview}"
        )
        if len(missing_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    protein_ids = [str(site_to_protein[site_id]).strip() for site_id in site_ids]
    invalid_site_ids = [
        site_id
        for site_id, protein_id in zip(site_ids, protein_ids, strict=True)
        if not protein_id
    ]
    if invalid_site_ids:
        preview = ", ".join(invalid_site_ids[:3])
        msg = (
            "site_to_protein must map aligned phosphosite rows to non-empty protein "
            f"IDs. Invalid mappings for: {preview}"
        )
        if len(invalid_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    series = pd.Series(
        protein_ids, index=pd.Index(site_ids, dtype=object), dtype=object
    )
    series.index.name = "site_id"
    series.name = "protein_id"
    return series


def _validate_supported_signalome_site_ids(site_ids: Sequence[str]) -> pd.Series:
    protein_ids: list[str] = []
    invalid_site_ids: list[str] = []
    for site_id in site_ids:
        protein_id = _protein_id_from_supported_site_id(site_id)
        if protein_id is None:
            invalid_site_ids.append(site_id)
            continue
        protein_ids.append(protein_id)

    if invalid_site_ids:
        preview = ", ".join(invalid_site_ids[:3])
        msg = (
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in the supported 'PROTEIN;SITE;...' "
            f"format. Invalid aligned site IDs: {preview}"
        )
        if len(invalid_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    series = pd.Series(
        protein_ids, index=pd.Index(site_ids, dtype=object), dtype=object
    )
    series.index.name = "site_id"
    series.name = "protein_id"
    return series


def _protein_id_from_supported_site_id(site_id: str) -> str | None:
    parts = [part.strip() for part in str(site_id).split(";")]
    if len(parts) < 3:
        return None
    protein_id, residue = parts[0], parts[1]
    if not protein_id or not residue:
        return None
    return protein_id


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
