from __future__ import annotations

from pydantic import BaseModel, ConfigDict, TypeAdapter, ValidationError

from .errors import RequestValidationError


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


__all__ = ["PhospyRequestModel", "validate_adapter_value"]
