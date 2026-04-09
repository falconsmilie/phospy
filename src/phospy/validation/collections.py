from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import pandas as pd

from .errors import PhospyValidationError


def resolve_required_columns(
    columns: Iterable[str],
    *,
    argument_name: str,
    context: str,
) -> list[str]:
    """Resolve a required list of column names and reject empty collections."""

    resolved_columns = list(columns)
    if not resolved_columns:
        raise PhospyValidationError(
            f"{context} requires at least one column name in '{argument_name}'"
        )
    return resolved_columns


def normalize_string_sequence(
    value: Sequence[str],
    *,
    field_name: str,
    empty_message: str,
    invalid_message: str | None = None,
    deduplicate: bool = False,
) -> tuple[str, ...]:
    """Normalize a sequence of string-like values for validated request fields."""

    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        msg = (
            invalid_message or f"{field_name} must be provided as a sequence of values"
        )
        raise ValueError(msg)

    normalized_values = tuple(str(item) for item in value)
    if deduplicate:
        normalized_values = tuple(dict.fromkeys(normalized_values))
    if not normalized_values:
        raise ValueError(empty_message)
    return normalized_values


def normalize_sequence_mapping(
    value: Mapping[str, Sequence[str]],
    *,
    field_name: str,
    empty_message: str,
) -> dict[str, tuple[str, ...]]:
    """Normalize mapping-backed sequence collections for validated requests."""

    if not isinstance(value, Mapping):
        msg = f"{field_name} must be provided as a mapping"
        raise ValueError(msg)
    if not value:
        raise ValueError(empty_message)

    normalized: dict[str, tuple[str, ...]] = {}
    for key, raw_values in value.items():
        if isinstance(raw_values, (str, bytes)):
            msg = (
                f"{field_name}[{key!r}] must be a sequence of values, "
                "not a plain string"
            )
            raise ValueError(msg)
        normalized[str(key)] = tuple(str(item) for item in raw_values)
    return normalized


def normalize_site_sequence_series(
    value: Mapping[str, str] | pd.Series | None,
) -> pd.Series | None:
    """Normalize site-sequence inputs to a detached Series keyed by site ID."""

    if value is None:
        return None
    if isinstance(value, pd.Series):
        normalized = value.copy(deep=True).astype(object)
        normalized.index = normalized.index.map(str)
        normalized[:] = normalized.map(str)
        return normalized
    if isinstance(value, Mapping):
        normalized = {
            str(site_id): str(sequence) for site_id, sequence in value.items()
        }
        return pd.Series(normalized, dtype=object)
    msg = (
        "site_sequences must be provided as a mapping keyed by phosphosite ID "
        "or as a pandas Series with an explicit phosphosite index; plain "
        "sequences are not supported"
    )
    raise ValueError(msg)


def normalize_site_to_protein_mapping(
    value: object,
) -> dict[str, str] | None:
    """Normalize site-to-protein mappings for validated signalome requests."""

    if value is None:
        return None
    if not isinstance(value, dict):
        try:
            value = dict(value)  # type: ignore[arg-type]
        except (TypeError, ValueError) as error:
            msg = (
                "site_to_protein must be provided as a mapping of site IDs to "
                "protein IDs"
            )
            raise ValueError(msg) from error

    normalized: dict[str, str] = {}
    for raw_site_id, raw_protein_id in value.items():
        site_id = str(raw_site_id)
        protein_id = str(raw_protein_id).strip()
        if not site_id:
            msg = "site_to_protein keys must be non-empty site IDs"
            raise ValueError(msg)
        if not protein_id:
            msg = "site_to_protein values must be non-empty protein IDs"
            raise ValueError(msg)
        normalized[site_id] = protein_id
    if not normalized:
        msg = "site_to_protein must contain at least one site-to-protein mapping"
        raise ValueError(msg)
    return normalized


__all__ = [
    "normalize_sequence_mapping",
    "normalize_site_sequence_series",
    "normalize_site_to_protein_mapping",
    "normalize_string_sequence",
    "resolve_required_columns",
]
