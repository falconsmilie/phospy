from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    field_validator,
    model_validator,
)

from ...datasets.models import PhosphoDataset
from ...datasets.schema import DatasetSchema
from ...errors import InputCompatibilityError, RequestValidationError
from ...internal.constants import ComparisonSpec
from ...matrices import SiteMatrixPolicy
from ...preprocessing.core import (
    CorePreprocessingConfig,
    resolve_core_preprocessing_config,
)
from ..domain import validate_dataset_comparisons
from ..schema.files import validate_existing_file_path
from ..schema.tables import PredMatSchema
from .analysis import AnalysisInputs, KinaseActivityRequest

if TYPE_CHECKING:
    from ...prediction.results import PredMatResult


class CorePipelineRequest(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

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
    site_matrix_policy: SiteMatrixPolicy = Field(default_factory=SiteMatrixPolicy)
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

    @field_validator("site_matrix_policy", mode="before")
    @classmethod
    def validate_site_matrix_policy(cls, value: object) -> SiteMatrixPolicy:
        try:
            return SiteMatrixPolicy.from_value(value)
        except (TypeError, ValueError) as error:
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
class PipelineInputs:
    """Trusted pipeline inputs owned by the pipeline boundary."""

    dataset: PhosphoDataset
    pred_mat: pd.DataFrame | None
    preprocessing_config: CorePreprocessingConfig
    kinase_activity_threshold: float | None
    kinase_activity_min_substrates: int | None
    kinase_activity_top_n_substrates: int | None


def build_pipeline_inputs(
    *,
    dataset: PhosphoDataset,
    pred_mat: pd.DataFrame | None = None,
    preprocessing_config: CorePreprocessingConfig | None = None,
    localization_threshold: float = 0.75,
    min_observed: int = 4,
    max_unmatched_fraction: float = 0.0,
    total_sentinel: float = 10.0,
    phospho_sentinel: float = 12.0,
    site_matrix_policy: SiteMatrixPolicy | Mapping[str, object] | None = None,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> PipelineInputs:
    """Build trusted pipeline inputs from already-owned dataset state."""

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
            site_matrix_policy=site_matrix_policy,
            context="Invalid pipeline construction request",
            config_param_name="preprocessing_config",
        )
    except (TypeError, ValueError) as error:
        raise RequestValidationError(str(error)) from error

    if pred_mat is not None and not isinstance(pred_mat, pd.DataFrame):
        msg = (
            "Invalid pipeline construction request: pred_mat must be a "
            "pandas DataFrame when provided"
        )
        raise RequestValidationError(msg)

    resolved_threshold: float | None = None
    resolved_min_substrates: int | None = None
    resolved_top_n_substrates: int | None = None
    if pred_mat is not None:
        kinase_activity_request = KinaseActivityRequest.validate_request(
            threshold=kinase_activity_threshold,
            min_substrates=kinase_activity_min_substrates,
            top_n_substrates=kinase_activity_top_n_substrates,
        )
        resolved_threshold = kinase_activity_request.threshold
        resolved_min_substrates = kinase_activity_request.min_substrates
        resolved_top_n_substrates = kinase_activity_request.top_n_substrates

    return PipelineInputs(
        dataset=dataset,
        pred_mat=pred_mat,
        preprocessing_config=resolved_config,
        kinase_activity_threshold=resolved_threshold,
        kinase_activity_min_substrates=resolved_min_substrates,
        kinase_activity_top_n_substrates=resolved_top_n_substrates,
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
    site_matrix_policy: SiteMatrixPolicy | Mapping[str, object] | None = None,
    kinase_activity_threshold: float = 0.6,
    kinase_activity_min_substrates: int = 3,
    kinase_activity_top_n_substrates: int = 20,
) -> PipelineInputs:
    """Validate raw in-memory inputs for pipeline construction."""

    from ...prediction.results import PredMatResult

    normalized_pred_mat = (
        pred_mat.to_frame(copy=False)
        if isinstance(pred_mat, PredMatResult)
        else pred_mat
    )
    validated_pred_mat = None
    if normalized_pred_mat is not None:
        validated_pred_mat = PredMatSchema.validate(
            normalized_pred_mat,
            context="pipeline pred_mat",
        )

    return build_pipeline_inputs(
        dataset=dataset,
        pred_mat=validated_pred_mat,
        preprocessing_config=preprocessing_config,
        localization_threshold=localization_threshold,
        min_observed=min_observed,
        max_unmatched_fraction=max_unmatched_fraction,
        total_sentinel=total_sentinel,
        phospho_sentinel=phospho_sentinel,
        site_matrix_policy=site_matrix_policy,
        kinase_activity_threshold=kinase_activity_threshold,
        kinase_activity_min_substrates=kinase_activity_min_substrates,
        kinase_activity_top_n_substrates=kinase_activity_top_n_substrates,
    )


def validate_pipeline_runtime_compatibility(
    *,
    request: PipelineInputs,
    site_matrix: pd.DataFrame,
) -> AnalysisInputs | None:
    """Validate post-preprocessing overlap before kinase analysis runs."""

    if (
        request.pred_mat is None
        or request.kinase_activity_threshold is None
        or request.kinase_activity_min_substrates is None
        or request.kinase_activity_top_n_substrates is None
    ):
        return None

    try:
        return AnalysisInputs.from_trusted_inputs(
            pred_mat=request.pred_mat,
            phospho_matrix=site_matrix,
            threshold=request.kinase_activity_threshold,
            min_substrates=request.kinase_activity_min_substrates,
            top_n_substrates=request.kinase_activity_top_n_substrates,
            pred_context="pipeline pred_mat",
            matrix_context="preprocessed site matrix",
        )
    except InputCompatibilityError as error:
        raise InputCompatibilityError(
            f"Pipeline runtime compatibility failed after preprocessing: {error}"
        ) from error
