from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from pydantic import Field, ValidationError, field_validator, model_validator

from ..constants import ComparisonSpec
from ..core_processing import (
    CorePreprocessingConfig,
    resolve_core_preprocessing_config,
)
from ..dataset_schema import DatasetSchema
from ._models import PhospyRequestModel
from .errors import InputCompatibilityError, RequestValidationError
from .paths import validate_existing_file_path
from .tables import PredMatSchema

if TYPE_CHECKING:
    from ..dataset import PhosphoDataset


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
            validated = self.dataset_schema.validate_comparisons(
                self.comparisons,
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
    """Trusted validated bundle for the public :class:`phospy.PhosRPipeline` boundary.

    The dataset workspace is already owned application state. Optional
    ``pred_mat`` tables arrive here only after the raw pipeline boundary has
    schema-validated and copied them once. Trusted pipeline orchestration then
    reuses those owned objects without copying again by default.
    """

    dataset: PhosphoDataset
    pred_mat: pd.DataFrame | None
    preprocessing_config: CorePreprocessingConfig


def build_validated_pipeline_request(
    *,
    dataset: PhosphoDataset,
    validated_pred_mat: pd.DataFrame | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
) -> ValidatedPipelineRequest:
    """Build a trusted pipeline request from already-owned validated inputs.

    This helper expects a dataset workspace and optional validated ``pred_mat``
    table that are already owned trusted state. It does not copy them again.
    Use :func:`validate_pipeline_request` at raw public boundaries.
    """

    from ..dataset import PhosphoDataset

    if not isinstance(dataset, PhosphoDataset):
        msg = "Invalid pipeline request: dataset must be a PhosphoDataset instance"
        raise RequestValidationError(msg)

    try:
        resolved_config = resolve_core_preprocessing_config(
            config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            context="Invalid pipeline request",
            config_param_name="preprocessing_config",
        )
    except (TypeError, ValueError) as error:
        raise RequestValidationError(str(error)) from error

    if validated_pred_mat is not None and not isinstance(
        validated_pred_mat, pd.DataFrame
    ):
        msg = "Invalid pipeline request: pred_mat must be a pandas DataFrame when provided"
        raise RequestValidationError(msg)

    return ValidatedPipelineRequest(
        dataset=dataset,
        pred_mat=validated_pred_mat,
        preprocessing_config=resolved_config,
    )


def validate_pipeline_request(
    *,
    dataset: PhosphoDataset,
    pred_mat: pd.DataFrame | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
) -> ValidatedPipelineRequest:
    """Validate raw in-memory pipeline inputs for the public pipeline boundary.

    The dataset argument is already an owned workspace. When a raw ``pred_mat``
    is supplied, this boundary takes ownership by schema-validating and copying
    it once before handing trusted state downstream.
    """

    validated_pred_mat = None
    if pred_mat is not None:
        validated_pred_mat = PredMatSchema.validate(pred_mat, context="pred_mat")

    return build_validated_pipeline_request(
        dataset=dataset,
        validated_pred_mat=validated_pred_mat,
        preprocessing_config=preprocessing_config,
        localization_threshold=localization_threshold,
        min_observed=min_observed,
        max_unmatched_fraction=max_unmatched_fraction,
        total_sentinel=total_sentinel,
        phospho_sentinel=phospho_sentinel,
    )


__all__ = [
    "CorePipelineRequest",
    "build_validated_pipeline_request",
    "ValidatedPipelineRequest",
    "validate_pipeline_request",
]
