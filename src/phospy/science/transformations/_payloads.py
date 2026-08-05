"""Transformation-domain provenance payload parsing helpers."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Final, cast

from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import InvalidTransformationStateError
from phospy.provenance.immutability import freeze_json_mapping_with_error_type
from phospy.provenance.models import TableFingerprint
from phospy.provenance.serialization import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)
from phospy.science.transformations.policy import QuantitativeMeaningEvidenceMode

_QUANTITATIVE_MEANING_OPERATION_PATTERN: Final[re.Pattern[str]] = re.compile(
    r"^[a-z0-9][a-z0-9_]*(?:\.[a-z0-9][a-z0-9_]*)*$"
)


def _normalize_required_provenance_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise InvalidTransformationStateError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise InvalidTransformationStateError(
            f"{field_name} must be a non-empty string"
        )
    return normalized


def _normalize_quantitative_meaning_evidence_mode(
    evidence_mode: QuantitativeMeaningEvidenceMode | str,
) -> QuantitativeMeaningEvidenceMode:
    raw_evidence_mode = cast(object, evidence_mode)
    if isinstance(raw_evidence_mode, QuantitativeMeaningEvidenceMode):
        return raw_evidence_mode
    try:
        return QuantitativeMeaningEvidenceMode(str(raw_evidence_mode).strip())
    except ValueError as exc:
        supported = ", ".join(
            member.value for member in QuantitativeMeaningEvidenceMode
        )
        raise InvalidTransformationStateError(
            "unsupported quantitative meaning evidence mode "
            f"{raw_evidence_mode!r}; supported: {supported}"
        ) from exc


def _normalize_quantitative_meaning_operation_id(value: object) -> str:
    raw = _normalize_required_provenance_text(
        value,
        field_name="quantitative_meaning_provenance.operation_id",
    )
    normalized = raw.strip().lower().replace("-", "_")
    normalized = re.sub(r"[\s_]+", "_", normalized)
    normalized = re.sub(r"\.+", ".", normalized)
    if not _QUANTITATIVE_MEANING_OPERATION_PATTERN.fullmatch(normalized):
        raise InvalidTransformationStateError(
            "quantitative_meaning_provenance.operation_id must be a stable "
            "dot-separated machine identifier using lowercase ASCII letters, "
            "digits, and underscores"
        )
    return normalized


def _normalize_diagnostic_caveat_codes(values: object) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise InvalidTransformationStateError(
            "quantitative_meaning_provenance.diagnostic_caveat_codes must be a "
            "sequence of strings"
        )
    if not isinstance(values, Sequence):
        raise InvalidTransformationStateError(
            "quantitative_meaning_provenance.diagnostic_caveat_codes must be a "
            "sequence of strings"
        )
    raw_values = tuple(cast(Sequence[object], values))
    codes: list[str] = []
    for code in raw_values:
        normalized = str(code).strip().lower().replace("-", "_").replace(" ", "_")
        normalized = re.sub(r"_+", "_", normalized)
        if not normalized:
            continue
        if not re.fullmatch(r"[a-z0-9][a-z0-9_]*", normalized):
            raise InvalidTransformationStateError(
                "quantitative_meaning_provenance.diagnostic_caveat_codes must "
                "contain stable lowercase ASCII caveat codes"
            )
        if normalized not in codes:
            codes.append(normalized)
    return tuple(codes)


def _canonical_table_fingerprint_payload_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[Mapping[str, object], ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise InvalidTransformationStateError(f"{field_name} must be a sequence")
    if not isinstance(values, Sequence):
        raise InvalidTransformationStateError(f"{field_name} must be a sequence")
    raw_values = tuple(cast(Sequence[object], values))
    return tuple(
        _canonical_table_fingerprint_payload(
            value,
            field_name=f"{field_name}[{position}]",
        )
        for position, value in enumerate(raw_values)
    )


def _canonical_table_fingerprint_payload(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    try:
        if isinstance(value, TableFingerprint):
            payload = table_fingerprint_to_payload(value)
        elif isinstance(value, Mapping):
            fingerprint = table_fingerprint_from_payload(
                cast(Mapping[str, object], value)
            )
            payload = table_fingerprint_to_payload(fingerprint)
        else:
            raise InvalidTransformationStateError(
                f"{field_name} must be a TableFingerprint or table fingerprint payload"
            )
    except PhosPyInputError as exc:
        raise InvalidTransformationStateError(
            f"{field_name} must be a valid table fingerprint payload: {exc}"
        ) from exc
    return cast(
        Mapping[str, object],
        freeze_json_mapping_with_error_type(
            payload,
            field_name=field_name,
            error_type=InvalidTransformationStateError,
        ),
    )


def _require_payload_mapping(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        raw_mapping = cast(Mapping[object, object], value)
        result: dict[str, object] = {}
        for key, item in raw_mapping.items():
            if not isinstance(key, str):
                raise InvalidTransformationStateError(
                    f"{field_name} JSON object keys must be strings; "
                    f"got {type(key).__name__}"
                )
            if key in result:
                raise InvalidTransformationStateError(
                    f"{field_name} contains duplicate JSON object key {key!r}"
                )
            result[key] = item
        return result
    raise InvalidTransformationStateError(f"{field_name} must be an object")


def _require_payload_sequence(value: object, *, field_name: str) -> list[object]:
    if isinstance(value, (str, bytes, bytearray)):
        raise InvalidTransformationStateError(f"{field_name} must be an array")
    if isinstance(value, Sequence):
        return list(cast(Sequence[object], value))
    raise InvalidTransformationStateError(f"{field_name} must be an array")


def _require_payload_str(value: object, *, field_name: str) -> str:
    return _normalize_required_provenance_text(value, field_name=field_name)


def _optional_payload_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_payload_str(value, field_name=field_name)


def _require_payload_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidTransformationStateError(f"{field_name} must be an int")
    return int(value)


__all__ = [
    "_canonical_table_fingerprint_payload",
    "_canonical_table_fingerprint_payload_tuple",
    "_normalize_diagnostic_caveat_codes",
    "_normalize_quantitative_meaning_evidence_mode",
    "_normalize_quantitative_meaning_operation_id",
    "_normalize_required_provenance_text",
    "_optional_payload_str",
    "_require_payload_int",
    "_require_payload_mapping",
    "_require_payload_sequence",
    "_require_payload_str",
]
