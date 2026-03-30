from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    field_validator,
    model_validator,
)

from phospy.constants import DEFAULT_TOTAL_COLS, ComparisonSpec
from phospy.types import DegenerateProbabilityPolicy, PredictionSvmMode

from .errors import RequestValidationError

_VALID_COMPARISON_GROUPS = frozenset(DEFAULT_TOTAL_COLS)


class PhospyRequestModel(BaseModel):
    """Base request model with shared validation behaviour."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)


class CorePipelineRequest(PhospyRequestModel):
    total_path: Path
    phospho_path: Path
    pred_mat_path: Path | None = None
    phospho_encoding: str | None = None
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
        if not value.exists():
            msg = f"Path does not exist: {value}"
            raise ValueError(msg)
        if not value.is_file():
            msg = f"Path is not a file: {value}"
            raise ValueError(msg)
        return value

    @field_validator("comparisons")
    @classmethod
    def validate_comparisons(
        cls,
        value: tuple[ComparisonSpec, ...] | None,
    ) -> tuple[ComparisonSpec, ...] | None:
        if value is None:
            return None

        seen: set[ComparisonSpec] = set()
        for left_group, right_group in value:
            if left_group not in _VALID_COMPARISON_GROUPS:
                msg = f"Unknown comparison group: {left_group}"
                raise ValueError(msg)
            if right_group not in _VALID_COMPARISON_GROUPS:
                msg = f"Unknown comparison group: {right_group}"
                raise ValueError(msg)
            pair = (left_group, right_group)
            if pair in seen:
                msg = f"Duplicate comparison pair: {left_group!r}, {right_group!r}"
                raise ValueError(msg)
            seen.add(pair)
        return value

    @classmethod
    def validate_request(cls, **data: object) -> CorePipelineRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid core pipeline request",
                error=error,
            ) from error


class KinaseActivityRequest(PhospyRequestModel):
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


class KinaseWorkflowRequest(PhospyRequestModel):
    phospho_matrix: pd.DataFrame
    substrate_map: Mapping[str, Sequence[str]]
    site_sequences: Mapping[str, str] | Sequence[str] | pd.Series | None = None
    motif_sequences: Mapping[str, Sequence[str]] | None = None
    min_substrates: int = Field(default=1, ge=1)
    min_motif_size: int = Field(default=1, ge=1)
    allow_profile_only_fallback: bool = False
    ensemble_size: int = Field(default=10, ge=1)
    top: int = Field(default=50, ge=1)
    score_threshold: float = Field(default=0.8, ge=0.0, le=1.0)
    inclusion: int = Field(default=20, ge=1)
    n_iterations: int = Field(default=5, ge=1)
    random_state: int | None = None
    svm_mode: PredictionSvmMode | None = None
    degenerate_probability_policy: DegenerateProbabilityPolicy | None = None

    @field_validator("substrate_map")
    @classmethod
    def validate_substrate_map(
        cls,
        value: Mapping[str, Sequence[str]],
    ) -> Mapping[str, Sequence[str]]:
        if not value:
            msg = "substrate_map must not be empty"
            raise ValueError(msg)
        return value

    @field_validator("site_sequences")
    @classmethod
    def validate_site_sequences(
        cls,
        value: Mapping[str, str] | Sequence[str] | pd.Series | None,
    ) -> Mapping[str, str] | Sequence[str] | pd.Series | None:
        if value is None or isinstance(value, (pd.Series, Mapping)):
            return value
        msg = (
            "site_sequences must be provided as a mapping keyed by phosphosite ID "
            "or as a pandas Series with an explicit phosphosite index"
        )
        raise ValueError(msg)

    @field_validator("motif_sequences")
    @classmethod
    def validate_motif_sequences(
        cls,
        value: Mapping[str, Sequence[str]] | None,
    ) -> Mapping[str, Sequence[str]] | None:
        if value is not None and not value:
            msg = (
                "motif_sequences must not be empty; pass None and set "
                "allow_profile_only_fallback=True for profile-only prediction"
            )
            raise ValueError(msg)
        return value

    @model_validator(mode="after")
    def validate_cross_field_requirements(self) -> KinaseWorkflowRequest:
        if self.motif_sequences is None and not self.allow_profile_only_fallback:
            msg = (
                "motif_sequences are required for end-to-end prediction unless "
                "allow_profile_only_fallback=True"
            )
            raise ValueError(msg)

        if self.motif_sequences is not None and self.site_sequences is None:
            msg = "site_sequences are required when motif_sequences are provided"
            raise ValueError(msg)

        return self

    @classmethod
    def validate_request(cls, **data: object) -> KinaseWorkflowRequest:
        try:
            return cls.model_validate(data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid kinase workflow request",
                error=error,
            ) from error
