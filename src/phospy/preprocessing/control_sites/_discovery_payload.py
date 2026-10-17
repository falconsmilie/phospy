"""Private scalar and payload validation for SPS discovery codecs."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.phosphosites.validation import require_site_key_index


def require_site_key(value: object, *, field_name: str) -> str:
    text = require_non_empty_string(value, field_name=field_name)
    if value != text:
        raise PhosPyInputError(
            f"{field_name} must be an exact canonical site_key without "
            "surrounding whitespace"
        )
    require_site_key_index(
        pd.Index([text], name="site_key"),
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    return text


def require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value.strip()


def require_positive_int(value: object, *, field_name: str) -> int:
    return require_int_at_least(value, minimum=1, field_name=field_name)


def require_non_negative_int(value: object, *, field_name: str) -> int:
    return require_int_at_least(value, minimum=0, field_name=field_name)


def require_int_at_least(value: object, *, minimum: int, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise PhosPyInputError(f"{field_name} must be an integer >= {minimum}")
    return value


def require_finite_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be a finite number")
    resolved = float(value)
    if not math.isfinite(resolved):
        raise PhosPyInputError(f"{field_name} must be a finite number")
    return resolved


def require_payload_mapping(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be an object")
    resolved: dict[str, object] = {}
    for key, item in cast(Mapping[object, object], value).items():
        if not isinstance(key, str):
            raise PhosPyInputError(f"{field_name} keys must be strings")
        resolved[key] = item
    return resolved


def require_payload_sequence(value: object, *, field_name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise PhosPyInputError(f"{field_name} must be an array")
    return tuple(cast(Sequence[object], value))


def require_payload_string(value: object, *, field_name: str) -> str:
    return require_non_empty_string(value, field_name=field_name)


def optional_payload_string(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return require_payload_string(value, field_name=field_name)


def require_payload_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an integer")
    return value


def optional_payload_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return require_payload_int(value, field_name=field_name)


def require_payload_float(value: object, *, field_name: str) -> float:
    return require_finite_float(value, field_name=field_name)


__all__: list[str] = []
