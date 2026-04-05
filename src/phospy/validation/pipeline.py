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
from .analysis import (
    KinaseActivityRequest,
    ValidatedAnalysisRequest,
)
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
class ValidatedPipelineConstructionRequest:
    """Trusted pipeline inputs after raw construction validation.

    The dataset workspace is already owned application state. Optional
    ``pred_mat`` tables arrive here only after the raw pipeline construction
    boundary has schema-validated and copied them once. Runtime compatibility
    with the eventual preprocessed site matrix is validated later, after core
    preprocessing produces that matrix.
    """

    dataset: PhosphoDataset
    pred_mat: pd.DataFrame | None
    preprocessing_config: CorePreprocessingConfig
    kinase_activity_request: KinaseActivityRequest | None


def build_validated_pipeline_construction_request(
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
) -> ValidatedPipelineConstructionRequest:
    """Build trusted pipeline construction inputs from already-owned state.

    This helper expects a dataset workspace and optional validated ``pred_mat``
    table that are already owned trusted state. It does not copy them again.
    Use :func:`validate_pipeline_construction_request` at raw public boundaries.
    Runtime compatibility with the preprocessed site matrix belongs to
    :func:`validate_pipeline_runtime_compatibility`.
    """

    from ..dataset import PhosphoDataset

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

    return ValidatedPipelineConstructionRequest(
        dataset=dataset,
        pred_mat=validated_pred_mat,
        preprocessing_config=resolved_config,
        kinase_activity_request=kinase_activity_request,
    )


def validate_pipeline_construction_request(
    *,
    dataset: PhosphoDataset,
    pred_mat: pd.DataFrame | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> ValidatedPipelineConstructionRequest:
    """Validate raw in-memory pipeline inputs for construction only.

    This raw boundary takes ownership by copying the caller-managed ``pred_mat``
    once during schema validation. It intentionally does not validate runtime
    compatibility against the eventual preprocessed site matrix, because that
    matrix does not exist until core preprocessing has run.
    """

    validated_pred_mat = None
    if pred_mat is not None:
        validated_pred_mat = PredMatSchema.validate(
            pred_mat,
            context="pipeline pred_mat",
        )

    return build_validated_pipeline_construction_request(
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
    request: ValidatedPipelineConstructionRequest,
    site_matrix: pd.DataFrame,
) -> ValidatedAnalysisRequest | None:
    """Validate post-preprocessing runtime compatibility for pipeline analysis.

    This second validation phase runs only after preprocessing has produced the
    site matrix needed to check overlap against ``pred_mat``. It returns a
    trusted downstream analysis request when kinase activity analysis should run.
    """

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


__all__ = [
    "CorePipelineRequest",
    "ValidatedPipelineConstructionRequest",
    "build_validated_pipeline_construction_request",
    "validate_pipeline_construction_request",
    "validate_pipeline_runtime_compatibility",
]
