"""Table-fingerprint provenance payload serialization."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.provenance.models import TableFingerprint
from phospy.provenance.serialization._payload import (
    optional_mapping,
    optional_str,
    raise_legacy_provenance_schema,
    reject_legacy_provenance_fields,
    require_int,
    require_mapping,
    require_sequence,
    require_str,
    to_json_safe,
    to_json_value,
)

_LEGACY_TABLE_FINGERPRINT_FIELDS = frozenset({"hash_algorithm", "hash_value"})


def table_fingerprint_to_payload(fingerprint: TableFingerprint) -> dict[str, object]:
    """Serialize a table fingerprint to a JSON-safe payload."""

    return {
        "name": fingerprint.name,
        "rows": int(fingerprint.rows),
        "columns": int(fingerprint.columns),
        "index_name": fingerprint.index_name,
        "column_names": list(fingerprint.column_names),
        "dtypes": list(fingerprint.dtypes),
        "exact_hash_algorithm": fingerprint.exact_hash_algorithm,
        "exact_hash_value": fingerprint.exact_hash_value,
        "tolerance_hash_algorithm": fingerprint.tolerance_hash_algorithm,
        "tolerance_hash_value": fingerprint.tolerance_hash_value,
        "index_structure": (
            None
            if fingerprint.index_structure is None
            else to_json_safe(fingerprint.index_structure)
        ),
        "column_index_structure": (
            None
            if fingerprint.column_index_structure is None
            else to_json_safe(fingerprint.column_index_structure)
        ),
    }


def table_fingerprint_from_payload(payload: Mapping[str, object]) -> TableFingerprint:
    """Deserialize a table fingerprint from a decoded payload."""

    payload = require_mapping(payload, field_name="table_fingerprint")
    reject_legacy_provenance_fields(
        payload,
        field_name="table_fingerprint",
        legacy_fields=_LEGACY_TABLE_FINGERPRINT_FIELDS,
    )
    index_structure = optional_mapping(
        payload.get("index_structure"),
        field_name="table_fingerprint.index_structure",
    )
    column_index_structure = optional_mapping(
        payload.get("column_index_structure"),
        field_name="table_fingerprint.column_index_structure",
    )
    exact_hash_algorithm = require_str(
        payload.get("exact_hash_algorithm"),
        field_name="table_fingerprint.exact_hash_algorithm",
    )
    exact_hash_value = require_str(
        payload.get("exact_hash_value"),
        field_name="table_fingerprint.exact_hash_value",
    )
    tolerance_hash_algorithm = require_str(
        payload.get("tolerance_hash_algorithm"),
        field_name="table_fingerprint.tolerance_hash_algorithm",
    )
    tolerance_hash_value = require_str(
        payload.get("tolerance_hash_value"),
        field_name="table_fingerprint.tolerance_hash_value",
    )
    return TableFingerprint(
        name=require_str(payload.get("name"), field_name="table_fingerprint.name"),
        rows=require_int(payload.get("rows"), field_name="table_fingerprint.rows"),
        columns=require_int(
            payload.get("columns"),
            field_name="table_fingerprint.columns",
        ),
        index_name=optional_str(
            payload.get("index_name"),
            field_name="table_fingerprint.index_name",
        ),
        column_names=tuple(
            require_str(item, field_name="table_fingerprint.column_names[]")
            for item in require_sequence(
                payload.get("column_names"),
                field_name="table_fingerprint.column_names",
            )
        ),
        dtypes=tuple(
            require_str(item, field_name="table_fingerprint.dtypes[]")
            for item in require_sequence(
                payload.get("dtypes"),
                field_name="table_fingerprint.dtypes",
            )
        ),
        exact_hash_algorithm=exact_hash_algorithm,
        exact_hash_value=exact_hash_value,
        tolerance_hash_algorithm=tolerance_hash_algorithm,
        tolerance_hash_value=tolerance_hash_value,
        index_structure=None
        if index_structure is None
        else {key: to_json_value(value) for key, value in index_structure.items()},
        column_index_structure=None
        if column_index_structure is None
        else {
            key: to_json_value(value) for key, value in column_index_structure.items()
        },
    )


def table_fingerprints_from_payload(
    value: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if value is None:
        raise_legacy_provenance_schema()
    payload = require_sequence(value, field_name=field_name)
    return tuple(
        table_fingerprint_from_payload(
            require_mapping(item, field_name=f"{field_name}[{position}]")
        )
        for position, item in enumerate(payload)
    )
