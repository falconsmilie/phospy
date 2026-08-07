"""JSON contracts and validation primitives for processing-state payloads."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TypeAlias, cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import (
    FrozenJsonMapping,
    FrozenJsonValue,
    freeze_json_mapping,
    freeze_json_value,
    thaw_json_mapping,
    thaw_json_value,
)

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
        "imputed_row_count",
        "imputed_column_count",
        "dropped_row_count",
        "random_seed",
        "method_parameters",
        "matrix_scale_requirement",
        "imputation_input_scale",
        "imputation_input_scale_source",
        "imputation_operation_order",
        "stage_order",
        "missingness_mask_hash",
        "imputation_mask_hash",
        "left_censored_assumption",
        "rows_not_imputable",
        "row_medians_used",
        "per_column_distribution_parameters",
        "dropped_rows_above_max_missing_fraction",
        "neighbour_count",
        "distance_metric",
        "knn_no_overlap_policy",
        "knn_no_overlap_policy_version",
        "knn_nearest_neighbour_imputed_cell_count",
        "knn_nearest_neighbour_imputed_row_ids",
        "knn_nearest_neighbour_imputed_column_ids",
        "knn_column_mean_fallback_imputed_cell_count",
        "knn_column_mean_fallback_row_ids",
        "knn_column_mean_fallback_column_ids",
        "knn_nearest_neighbour_imputation_mask_hash",
        "knn_column_mean_fallback_imputation_mask_hash",
        "knn_fully_column_mean_fallback_row_ids",
        "diagnostic_caveat_codes",
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
    seen: set[str] = set()
    for key in value:
        if not isinstance(key, str):
            raise PhosPyInputError(
                f"{field_name} must contain only string keys; got key "
                f"{key!r} ({type(key).__name__})"
            )
        if key in seen:
            raise PhosPyInputError(
                f"{field_name} contains duplicate JSON object key {key!r}"
            )
        seen.add(key)


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
    frozen = require_optional_frozen_string_to_string_mapping(
        value,
        field_name=field_name,
    )
    if frozen is None:
        return None
    return cast(dict[str, str], thaw_json_mapping(frozen, field_name=field_name))


def require_optional_frozen_string_to_string_mapping(
    value: object,
    *,
    field_name: str,
) -> FrozenJsonMapping | None:
    if value is None:
        return None
    mapping = require_mapping(value, field_name=field_name)
    frozen = freeze_json_mapping(mapping, field_name=field_name)
    parsed: list[tuple[str, str]] = []
    seen: set[str] = set()
    for raw_key, raw_value in frozen.items():
        key = require_required_str(
            raw_key,
            field_name=f"{field_name}.<key>",
        )
        if key in seen:
            raise PhosPyInputError(
                f"{field_name} contains duplicate JSON object key {key!r}"
            )
        seen.add(key)
        parsed.append(
            (
                key,
                require_required_str(
                    raw_value,
                    field_name=f"{field_name}.{key}",
                ),
            )
        )
    return FrozenJsonMapping(parsed, field_name=field_name)


def require_optional_string_to_float_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, float] | None:
    frozen = require_optional_frozen_string_to_float_mapping(
        value,
        field_name=field_name,
    )
    if frozen is None:
        return None
    return cast(dict[str, float], thaw_json_mapping(frozen, field_name=field_name))


def require_optional_frozen_string_to_float_mapping(
    value: object,
    *,
    field_name: str,
) -> FrozenJsonMapping | None:
    if value is None:
        return None
    mapping = require_mapping(value, field_name=field_name)
    frozen = freeze_json_mapping(mapping, field_name=field_name)
    parsed: list[tuple[str, float]] = []
    seen: set[str] = set()
    for raw_key, raw_value in frozen.items():
        key = require_required_str(
            raw_key,
            field_name=f"{field_name}.<key>",
        )
        if key in seen:
            raise PhosPyInputError(
                f"{field_name} contains duplicate JSON object key {key!r}"
            )
        seen.add(key)
        if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float)):
            raise PhosPyInputError(f"{field_name}.{key} must be a float")
        parsed.append((key, float(raw_value)))
    return FrozenJsonMapping(parsed, field_name=field_name)


def require_json_value(value: object, *, field_name: str) -> JsonValue:
    frozen = freeze_json_value(value, field_name=field_name)
    return cast(JsonValue, thaw_json_value(frozen, field_name=field_name))


def require_frozen_json_mapping(
    value: object,
    *,
    field_name: str,
) -> FrozenJsonMapping:
    mapping = require_mapping(value, field_name=field_name)
    if isinstance(mapping, FrozenJsonMapping):
        return mapping
    return freeze_json_mapping(mapping, field_name=field_name)


def require_optional_frozen_json_mapping(
    value: object,
    *,
    field_name: str,
) -> FrozenJsonMapping | None:
    if value is None:
        return None
    return require_frozen_json_mapping(value, field_name=field_name)


def thaw_frozen_json_mapping(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> dict[str, JsonValue]:
    return cast(
        dict[str, JsonValue],
        thaw_json_mapping(value, field_name=field_name),
    )


def thaw_frozen_json_value(
    value: FrozenJsonValue,
    *,
    field_name: str,
) -> JsonValue:
    return cast(
        JsonValue,
        thaw_json_value(value, field_name=field_name),
    )


def require_json_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, JsonValue]:
    return thaw_frozen_json_mapping(
        require_frozen_json_mapping(value, field_name=field_name),
        field_name=field_name,
    )


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
