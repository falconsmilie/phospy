from __future__ import annotations

from collections.abc import Mapping

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import require_str


def _parse_optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be an int")
    if isinstance(value, int):
        return int(value)
    raise PhosPyInputError(f"{field_name} must be an int")


def _resolve_optional_alias_string(
    *,
    payload: Mapping[str, object],
    canonical_key: str,
    alias_key: str,
    scope: str,
) -> str | None:
    canonical_value = payload.get(canonical_key)
    alias_value = payload.get(alias_key)
    if canonical_value is None and alias_value is None:
        return None
    canonical = (
        None
        if canonical_value is None
        else require_str(
            canonical_value,
            field_name=f"{scope}.signalome_config.{canonical_key}",
        )
    )
    alias = (
        None
        if alias_value is None
        else require_str(
            alias_value,
            field_name=f"{scope}.signalome_config.{alias_key}",
        )
    )
    if canonical is not None and alias is not None and canonical != alias:
        raise PhosPyInputError(
            f"{scope}.signalome_config.{canonical_key} conflicts with "
            f"{scope}.signalome_config.{alias_key}; provide matching values."
        )
    return canonical if canonical is not None else alias


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be an int")
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float) and float(value).is_integer():
        return int(value)
    raise PhosPyInputError(f"{field_name} must be an int")


def _reject_unsupported_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    allowed_fields: frozenset[str],
) -> None:
    unknown_fields = sorted(
        str(key) for key in payload.keys() if str(key) not in allowed_fields
    )
    if unknown_fields:
        unknown = ", ".join(unknown_fields)
        raise PhosPyInputError(f"{field_name} contains unsupported field(s): {unknown}")


def _require_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    required_fields: frozenset[str],
) -> None:
    missing_fields = sorted(
        str(key) for key in required_fields if str(key) not in payload
    )
    if not missing_fields:
        return
    missing = ", ".join(missing_fields)
    raise PhosPyInputError(f"{field_name} is missing required field(s): {missing}")
