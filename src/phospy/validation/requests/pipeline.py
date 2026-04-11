from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from pydantic import Field, ValidationError, field_validator, model_validator

from ...constants import ComparisonSpec
from ...datasets.schema import DatasetSchema
from ...preprocessing.core import (
    CorePreprocessingConfig,
    resolve_core_preprocessing_config,
)
from ..domain import validate_dataset_comparisons
from ..errors import InputCompatibilityError, RequestValidationError
from ..schema.files import validate_existing_file_path
from ..schema.tables import PredMatSchema
from .analysis import KinaseActivityRequest, ValidatedAnalysisRequest
from .shared import PhospyRequestModel, normalize_pred_mat_input

if TYPE_CHECKING:
    from ...datasets.models import PhosphoDataset
    from ...prediction.models import PredMatResult


class CorePipelineRequest(PhospyRequestModel):
    """Validated file-backed boundary request for pipeline construction."""

    total_path: Path
    phospho_path: Path
    pred_mat_path: Path | None = None
    phospho_encoding: str | None = None
    dataset_schema: DatasetSchema = Field(default_factory=DatasetSchema, alias="schema")
    comparisons: tuple[ComparisonSpec, ...] | None = None
    localization_threshold: float = Field(default=0.75, ge=0.0, le=1.0)
    min_observed: int = Field(default=4, ge=1)
    total_sentinel: float = 10.0
    phospho_sentinel: float = 12.0
    max_unmatched_fraction: float = Field(default=0.0, ge=0.0, le=1.0)
    kinase_activity_threshold: float = Field(default=0.6, ge=0.0, le=1.0)
    kinase_activity_min_substrates: int = Field(default=3, ge=1)
    kinase_activity_top_n_substrates: int = Field(default=20, ge=1)

    @field_validator("total_path", "phospho_path", "pred_mat_path")
    @classmethod
    def validate_existing_path(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        try:
            return validate_existing_file_path(value, context="core pipeline file path")
        except RequestValidationError as error:
            raise ValueError(str(error)) from error

    @model_validator(mode="after")
    def validate_comparisons(self) -> CorePipelineRequest:
        try:
            validated = validate_dataset_comparisons(
                schema=self.dataset_schema,
                comparisons=self.comparisons,
                context="Core pipeline request",
            )
        except (InputCompatibilityError, TypeError, ValueError) as error:
            raise ValueError(str(error)) from error
        object.__setattr__(self, "comparisons", validated)
        return self

    @classmethod
    def validate_request(cls, **data: object) -> CorePipelineRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid core pipeline request",
                error=error,
            ) from error


@dataclass(slots=True)
class ValidatedPipelineRequest:
    """Trusted pipeline inputs owned by the pipeline boundary."""

    dataset: PhosphoDataset
    pred_mat: pd.DataFrame | None
    preprocessing_config: CorePreprocessingConfig
    kinase_activity_request: KinaseActivityRequest | None


def build_pipeline_request(
    *,
    dataset: PhosphoDataset,
    validated_pred_mat: pd.DataFrame | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> ValidatedPipelineRequest:
    """Build a trusted pipeline request from already-owned inputs."""

    from ...datasets.models import PhosphoDataset

    if not isinstance(dataset, PhosphoDataset):
        msg = (
            "Invalid pipeline construction request: dataset must be a "
            "PhosphoDataset instance"
        )
        raise RequestValidationError(msg)

    try:
        resolved_config = resolve_core_preprocessing_config(
            config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            context="Invalid pipeline construction request",
            config_param_name="preprocessing_config",
        )
    except (TypeError, ValueError) as error:
        raise RequestValidationError(str(error)) from error

    if validated_pred_mat is not None and not isinstance(
        validated_pred_mat, pd.DataFrame
    ):
        msg = (
            "Invalid pipeline construction request: pred_mat must be a "
            "pandas DataFrame when provided"
        )
        raise RequestValidationError(msg)

    kinase_activity_request = None
    if validated_pred_mat is not None:
        kinase_activity_request = KinaseActivityRequest.validate_request(
            threshold=kinase_activity_threshold,
            min_substrates=kinase_activity_min_substrates,
            top_n_substrates=kinase_activity_top_n_substrates,
        )

    return ValidatedPipelineRequest(
        dataset=dataset,
        pred_mat=validated_pred_mat,
        preprocessing_config=resolved_config,
        kinase_activity_request=kinase_activity_request,
    )


def validate_pipeline_construction_request(
    *,
    dataset: PhosphoDataset,
    pred_mat: pd.DataFrame | PredMatResult | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> ValidatedPipelineRequest:
    """Validate raw in-memory inputs for pipeline construction only."""

    normalized_pred_mat = normalize_pred_mat_input(pred_mat)
    validated_pred_mat = None
    if normalized_pred_mat is not None:
        validated_pred_mat = PredMatSchema.validate(
            normalized_pred_mat,
            context="pipeline pred_mat",
        )

    return build_pipeline_request(
        dataset=dataset,
        validated_pred_mat=validated_pred_mat,
        preprocessing_config=preprocessing_config,
        localization_threshold=localization_threshold,
        min_observed=min_observed,
        max_unmatched_fraction=max_unmatched_fraction,
        total_sentinel=total_sentinel,
        phospho_sentinel=phospho_sentinel,
        kinase_activity_threshold=kinase_activity_threshold,
        kinase_activity_min_substrates=kinase_activity_min_substrates,
        kinase_activity_top_n_substrates=kinase_activity_top_n_substrates,
    )


def validate_pipeline_runtime_compatibility(
    *,
    request: ValidatedPipelineRequest,
    site_matrix: pd.DataFrame,
) -> ValidatedAnalysisRequest | None:
    """Validate post-preprocessing overlap before kinase analysis runs."""

    if request.pred_mat is None or request.kinase_activity_request is None:
        return None

    try:
        return ValidatedAnalysisRequest.from_trusted_inputs(
            request=request.kinase_activity_request,
            pred_mat=request.pred_mat,
            phospho_matrix=site_matrix,
            pred_context="pipeline pred_mat",
            matrix_context="preprocessed site matrix",
        )
    except InputCompatibilityError as error:
        raise InputCompatibilityError(
            f"Pipeline runtime compatibility failed after preprocessing: {error}"
        ) from error
