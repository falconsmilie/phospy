from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import pandas as pd
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    TypeAdapter,
    ValidationError,
    field_validator,
    model_validator,
)

from ...errors import RequestValidationError
from ...internal.defaults import DEFAULT_PREDICTION_DEBUG_TOP_N
from ...internal.types import (
    PREDICTION_TRACE_FORMAT_CSV,
    PREDICTION_TRACE_LEVEL_FULL,
    PREDICTION_TRACE_LEVEL_NONE,
    PREDICTION_TRACE_LEVEL_SUMMARY,
    PredictionSvmMode,
    PredictionTraceFormat,
    PredictionTraceLevel,
)
from ..schema.tables import PredictionScoreMatrixSchema

_PREDICTION_SVM_MODE_ADAPTER = TypeAdapter(PredictionSvmMode)
_PREDICTION_TRACE_FORMAT_ADAPTER = TypeAdapter(PredictionTraceFormat)
_PREDICTION_TRACE_LEVEL_ADAPTER = TypeAdapter(PredictionTraceLevel)


class PredictionRequest(BaseModel):
    """Raw boundary request for prediction execution."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    combined_scores: pd.DataFrame
    ensemble_size: int = Field(ge=1)
    top: int = Field(ge=1)
    score_threshold: float = Field(ge=0.0, le=1.0)
    inclusion: int = Field(ge=1)
    n_iterations: int = Field(ge=1)
    random_state: int | None = None
    debug_kinases: tuple[str, ...] | None = None
    debug_top_n: int = Field(default=DEFAULT_PREDICTION_DEBUG_TOP_N, ge=1)
    svm_mode: PredictionSvmMode
    sampling_trace: Any | None = None
    trace_level: PredictionTraceLevel = PREDICTION_TRACE_LEVEL_NONE
    trace_sink_format: PredictionTraceFormat = PREDICTION_TRACE_FORMAT_CSV
    trace_sink: Any | None = None

    @field_validator("combined_scores")
    @classmethod
    def validate_combined_scores(cls, value: pd.DataFrame) -> pd.DataFrame:
        return PredictionScoreMatrixSchema.validate(
            value,
            context="combined_scores",
        )

    @field_validator("debug_kinases", mode="before")
    @classmethod
    def normalize_debug_kinases(
        cls,
        value: Sequence[str] | None,
    ) -> tuple[str, ...] | None:
        if value is None:
            return None
        return tuple(value)

    @model_validator(mode="after")
    def validate_trace_sink_requirements(self) -> PredictionRequest:
        if (
            self.trace_level != PREDICTION_TRACE_LEVEL_FULL
            and self.trace_sink is not None
        ):
            msg = "trace_sink may only be provided when trace_level='full'"
            raise ValueError(msg)
        return self

    @classmethod
    def validate_request(
        cls,
        *,
        default_svm_mode: PredictionSvmMode,
        capture_debug_trace: bool = False,
        trace_sink_format: PredictionTraceFormat = PREDICTION_TRACE_FORMAT_CSV,
        **data: object,
    ) -> PredictionRequest:
        try:
            resolved_trace_level = _PREDICTION_TRACE_LEVEL_ADAPTER.validate_python(
                PREDICTION_TRACE_LEVEL_SUMMARY
                if data.get("trace_level") is None and capture_debug_trace
                else data.get("trace_level") or PREDICTION_TRACE_LEVEL_NONE
            )
        except ValidationError as error:
            details = error.errors(include_url=False)
            message = (
                str(details[0].get("msg", "Invalid value")) if details else str(error)
            )
            raise RequestValidationError(
                f"Invalid prediction request: trace_level: {message}"
            ) from error

        try:
            resolved_trace_format = _PREDICTION_TRACE_FORMAT_ADAPTER.validate_python(
                trace_sink_format
            )
        except ValidationError as error:
            details = error.errors(include_url=False)
            message = (
                str(details[0].get("msg", "Invalid value")) if details else str(error)
            )
            raise RequestValidationError(
                f"Invalid prediction request: trace_sink_format: {message}"
            ) from error

        try:
            resolved_svm_mode = _PREDICTION_SVM_MODE_ADAPTER.validate_python(
                default_svm_mode if data.get("svm_mode") is None else data["svm_mode"]
            )
        except ValidationError as error:
            details = error.errors(include_url=False)
            message = (
                str(details[0].get("msg", "Invalid value")) if details else str(error)
            )
            raise RequestValidationError(
                f"Invalid prediction request: svm_mode: {message}"
            ) from error

        request_data = dict(data)
        request_data["svm_mode"] = resolved_svm_mode
        request_data["trace_level"] = resolved_trace_level
        request_data["trace_sink_format"] = resolved_trace_format

        try:
            return cls.model_validate(request_data)
        except ValidationError as error:
            raise RequestValidationError.from_pydantic(
                context="Invalid prediction request",
                error=error,
            ) from error
