from __future__ import annotations

from pathlib import Path

from pydantic import Field, ValidationError, field_validator, model_validator

from ..constants import ComparisonSpec
from ..dataset_schema import DatasetSchema
from ._models import PhospyRequestModel
from .errors import InputCompatibilityError, RequestValidationError
from .paths import validate_existing_file_path


class CorePipelineRequest(PhospyRequestModel):
    """Validated boundary request for file-backed pipeline construction."""

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


__all__ = ["CorePipelineRequest"]
