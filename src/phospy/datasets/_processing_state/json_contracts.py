"""JSON contracts and validation primitives for processing-state payloads."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TypeAlias

from phospy.errors.input import PhosPyInputError

JsonPrimitive: TypeAlias = None | str | bool | int | float
JsonValue: TypeAlias = JsonPrimitive | list["JsonValue"] | dict[str, "JsonValue"]

MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1 = 1
TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1 = 1

V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS = frozenset(
    (
        "diagnostics_schema_version",
        "missing_data_policy",
        "imputation_method_id",
        "imputation_method_family",
        "input_missing_cell_count",
        "output_missing_cell_count",
        "imputed_cell_count",
        "affected_row_count",
        "affected_column_count",
        "affected_row_ids",
        "affected_column_ids",
        "imputed_row_ids",
        "imputed_column_ids",
        "dropped_row_ids",
        "random_seed",
        "method_parameters",
        "matrix_scale_requirement",
        "stage_order",
        "missingness_mask_hash",
        "left_censored_assumption",
        "rows_not_imputable",
        "row_medians_used",
        "per_column_distribution_parameters",
        "dropped_rows_above_max_missing_fraction",
        "neighbour_count",
        "distance_metric",
    )
)

V1_KNOWN_TOTAL_PROTEIN_DIAGNOSTICS_FIELDS = frozenset(
    (
        "diagnostics_schema_version",
        "policy",
        "requested_policy",
        "resolved_policy",
        "formula",
        "requires_log_scale",
        "input_scale",
        "output_scale",
        "quantitative_meaning",
        "matched_rows",
        "identity_mode",
        "identity_matching_policy",
        "phosphosite_key",
        "total_protein_key",
        "mapping_phosphosite_key",
        "mapping_total_protein_key",
        "mapping_table_fingerprint",
        "duplicate_policy",
        "unmatched_policy",
        "phosphosite_row_count",
        "total_protein_row_count",
        "corrected_row_count",
        "uncorrected_row_count",
        "unused_total_protein_row_count",
        "total_rows_used_by_multiple_phosphosites",
        "corrected_phosphosite_row_ids",
        "corrected_phosphosite_to_total_protein_row_id",
        "unmatched_phosphosite_row_ids",
        "uncorrected_phosphosite_row_reasons",
        "unused_total_protein_row_ids",
        "gene_symbol_matching_used",
        "gene_symbol_identity_warning",
        "total_table_hash",
        "input_phospho_hash",
        "output_phospho_hash",
    )
)


def require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be an object")
    return value


def require_string_keys(value: Mapping[str, object], *, field_name: str) -> None:
    for key in value:
        if not isinstance(key, str):
            raise PhosPyInputError(
                f"{field_name} must contain only string keys; got key "
                f"{key!r} ({type(key).__name__})"
            )


def require_optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def require_required_str(value: object, *, field_name: str) -> str:
    if value is None:
        raise PhosPyInputError(f"{field_name} is required")
    parsed = require_optional_str(value, field_name=field_name)
    if parsed is None:  # pragma: no cover - defensive guard
        raise PhosPyInputError(f"{field_name} is required")
    return parsed


def require_optional_bool(value: object, *, field_name: str) -> bool | None:
    if value is None:
        return None
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return value


def require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an int")
    return value


def require_optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return require_int(value, field_name=field_name)


def require_optional_non_negative_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    parsed = require_int(value, field_name=field_name)
    if parsed < 0:
        raise PhosPyInputError(f"{field_name} must be >= 0")
    return parsed


def require_required_non_negative_int(value: object, *, field_name: str) -> int:
    if value is None:
        raise PhosPyInputError(f"{field_name} is required")
    parsed = require_optional_non_negative_int(value, field_name=field_name)
    if parsed is None:  # pragma: no cover - defensive guard
        raise PhosPyInputError(f"{field_name} is required")
    return parsed


def require_optional_string_tuple(
    value: object, *, field_name: str
) -> tuple[str, ...] | None:
    if value is None:
        return None
    if not isinstance(value, (list, tuple)):
        raise PhosPyInputError(f"{field_name} must be an array of strings")
    parsed: list[str] = []
    for position, item in enumerate(value):
        parsed.append(
            require_required_str(
                item,
                field_name=f"{field_name}[{position}]",
            )
        )
    return tuple(parsed)


def require_required_string_tuple(value: object, *, field_name: str) -> tuple[str, ...]:
    if value is None:
        raise PhosPyInputError(f"{field_name} is required")
    parsed = require_optional_string_tuple(value, field_name=field_name)
    if parsed is None:  # pragma: no cover - defensive guard
        raise PhosPyInputError(f"{field_name} is required")
    return parsed


def require_optional_string_to_string_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, str] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be an object of string mappings")
    parsed: dict[str, str] = {}
    for raw_key, raw_value in value.items():
        key = require_required_str(
            raw_key,
            field_name=f"{field_name}.<key>",
        )
        parsed[key] = require_required_str(
            raw_value,
            field_name=f"{field_name}.{key}",
        )
    return parsed


def require_optional_string_to_float_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, float] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be an object of numeric mappings")
    parsed: dict[str, float] = {}
    for raw_key, raw_value in value.items():
        key = require_required_str(
            raw_key,
            field_name=f"{field_name}.<key>",
        )
        if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float)):
            raise PhosPyInputError(f"{field_name}.{key} must be a float")
        parsed[key] = float(raw_value)
    return parsed


def require_json_value(value: object, *, field_name: str) -> JsonValue:
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, (list, tuple)):
        return [
            require_json_value(item, field_name=f"{field_name}[]") for item in value
        ]
    if isinstance(value, Mapping):
        normalized: dict[str, JsonValue] = {}
        for raw_key, raw_value in value.items():
            key = require_required_str(raw_key, field_name=f"{field_name}.<key>")
            normalized[key] = require_json_value(
                raw_value, field_name=f"{field_name}.{key}"
            )
        return normalized
    raise PhosPyInputError(
        f"{field_name} must be JSON-compatible (null, bool, int, float, string, array, or object)"
    )


def require_json_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, JsonValue]:
    mapping = require_mapping(value, field_name=field_name)
    normalized: dict[str, JsonValue] = {}
    for raw_key, raw_value in mapping.items():
        key = require_required_str(raw_key, field_name=f"{field_name}.<key>")
        normalized[key] = require_json_value(
            raw_value, field_name=f"{field_name}.{key}"
        )
    return normalized


def require_optional_json_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, JsonValue] | None:
    if value is None:
        return None
    return require_json_mapping(value, field_name=field_name)


def set_optional_payload_value(
    payload: dict[str, JsonValue],
    key: str,
    value: JsonValue,
) -> None:
    if value is not None:
        payload[key] = value
