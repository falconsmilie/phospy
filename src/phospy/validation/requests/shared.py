from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
from pydantic import BaseModel, ConfigDict, TypeAdapter, ValidationError

from ...errors import RequestValidationError

if TYPE_CHECKING:
    from ...prediction.results import PredMatResult


class PhospyRequestModel(BaseModel):
    """Base request model with shared validation behaviour."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)


def validate_adapter_value(
    *,
    value: object,
    adapter: TypeAdapter[object],
    field_name: str,
    context: str,
) -> object:
    """Validate an adapter-backed field and raise a package-level request error."""

    try:
        return adapter.validate_python(value)
    except ValidationError as error:
        details = error.errors(include_url=False)
        message = str(details[0].get("msg", "Invalid value")) if details else str(error)
        raise RequestValidationError(f"{context}: {field_name}: {message}") from error


def normalize_pred_mat_input(
    pred_mat: pd.DataFrame | PredMatResult | None,
) -> pd.DataFrame | None:
    """Normalize public predMat inputs to the internal DataFrame contract."""

    from ...prediction.results import PredMatResult

    if isinstance(pred_mat, PredMatResult):
        return pred_mat.to_frame(copy=False)
    return pred_mat


__all__ = [
    "PhospyRequestModel",
    "normalize_pred_mat_input",
    "validate_adapter_value",
]
