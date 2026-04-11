from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import pandas as pd
from pydantic import (
    Field,
    TypeAdapter,
    ValidationError,
    field_validator,
    model_validator,
)

from ...errors import RequestValidationError
from ...internal.types import (
    PredictionSvmMode,
    PredictionTraceFormat,
    PredictionTraceLevel,
)
from ..schema.tables import PredictionScoreMatrixSchema
from .shared import PhospyRequestModel, validate_adapter_value

_PREDICTION_SVM_MODE_ADAPTER = TypeAdapter(PredictionSvmMode)
_PREDICTION_TRACE_FORMAT_ADAPTER = TypeAdapter(PredictionTraceFormat)
_PREDICTION_TRACE_LEVEL_ADAPTER = TypeAdapter(PredictionTraceLevel)


class PredictionRequest(PhospyRequestModel):
    """Validated boundary request for prediction execution."""

    combined_scores: pd.DataFrame
    ensemble_size: int = Field(ge=1)
    top: int = Field(ge=1)
    score_threshold: float = Field(ge=0.0, le=1.0)
    inclusion: int = Field(ge=1)
    n_iterations: int = Field(ge=1)
    random_state: int | None = None
    debug_kinases: tuple[str, ...] | None = None
    debug_top_n: int = Field(default=10, ge=1)
    svm_mode: PredictionSvmMode
    sampling_trace: Any | None = None
    trace_level: PredictionTraceLevel = "none"
    trace_sink_format: PredictionTraceFormat = "csv"
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
        if self.trace_level != "full" and self.trace_sink is not None:
            msg = "trace_sink may only be provided when trace_level='full'"
            raise ValueError(msg)
        return self

    @classmethod
    def validate_request(
        cls,
        *,
        default_svm_mode: PredictionSvmMode,
        capture_debug_trace: bool = False,
        trace_sink_format: PredictionTraceFormat = "csv",
        **data: object,
    ) -> PredictionRequest:
        resolved_trace_level = validate_adapter_value(
            value=(
                "summary"
                if data.get("trace_level") is None and capture_debug_trace
                else data.get("trace_level") or "none"
            ),
            adapter=_PREDICTION_TRACE_LEVEL_ADAPTER,
            field_name="trace_level",
            context="Invalid prediction request",
        )
        resolved_trace_format = validate_adapter_value(
            value=trace_sink_format,
            adapter=_PREDICTION_TRACE_FORMAT_ADAPTER,
            field_name="trace_sink_format",
            context="Invalid prediction request",
        )
        resolved_svm_mode = validate_adapter_value(
            value=(
                default_svm_mode if data.get("svm_mode") is None else data["svm_mode"]
            ),
            adapter=_PREDICTION_SVM_MODE_ADAPTER,
            field_name="svm_mode",
            context="Invalid prediction request",
        )

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
