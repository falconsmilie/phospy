"""Protein-scoped phosphosite key model and reversible encoding helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar
from urllib.parse import quote, unquote

ErrorType = TypeVar("ErrorType", bound=Exception)

_ENCODING_VERSION_PREFIX = "phospy:v1"
_REQUIRED_FIELDS = (
    "organism",
    "protein_namespace",
    "protein_identifier",
    "residue",
    "position",
)
_OPTIONAL_FIELDS = ("isoform_id",)
_ENCODING_FIELD_ORDER = (*_REQUIRED_FIELDS, *_OPTIONAL_FIELDS)
_VALID_RESIDUES = {"S", "T", "Y"}


@dataclass(frozen=True, slots=True)
class ProteinScopedPhosphositeKey:
    organism: str
    protein_namespace: str
    protein_identifier: str
    residue: str
    position: int
    isoform_id: str | None = None


def build_protein_scoped_site_key(
    *,
    organism: object,
    protein_namespace: object,
    protein_identifier: object,
    residue: object,
    position: object,
    isoform_id: object = None,
    field_name: str,
    error_type: type[ErrorType],
) -> ProteinScopedPhosphositeKey:
    return ProteinScopedPhosphositeKey(
        organism=_required_text(
            organism,
            field_name=f"{field_name}.organism",
            error_type=error_type,
        ),
        protein_namespace=_required_text(
            protein_namespace,
            field_name=f"{field_name}.protein_namespace",
            error_type=error_type,
        ),
        protein_identifier=_required_text(
            protein_identifier,
            field_name=f"{field_name}.protein_identifier",
            error_type=error_type,
        ),
        residue=_residue_token(
            residue,
            field_name=f"{field_name}.residue",
            error_type=error_type,
        ),
        position=_positive_integer(
            position,
            field_name=f"{field_name}.position",
            error_type=error_type,
        ),
        isoform_id=_optional_text(
            isoform_id,
            field_name=f"{field_name}.isoform_id",
            error_type=error_type,
        ),
    )


def encode_site_key(key: ProteinScopedPhosphositeKey) -> str:
    canonical = build_protein_scoped_site_key(
        organism=key.organism,
        protein_namespace=key.protein_namespace,
        protein_identifier=key.protein_identifier,
        residue=key.residue,
        position=key.position,
        isoform_id=key.isoform_id,
        field_name="site_key",
        error_type=ValueError,
    )
    payload = {
        "organism": canonical.organism,
        "protein_namespace": canonical.protein_namespace,
        "protein_identifier": canonical.protein_identifier,
        "residue": canonical.residue,
        "position": str(canonical.position),
    }
    if canonical.isoform_id is not None:
        payload["isoform_id"] = canonical.isoform_id

    encoded_parts = [_ENCODING_VERSION_PREFIX]
    for field_name in _ENCODING_FIELD_ORDER:
        if field_name not in payload:
            continue
        raw_value = payload[field_name]
        encoded_value = quote(str(raw_value), safe="")
        encoded_parts.append(f"{field_name}={encoded_value}")
    return "|".join(encoded_parts)


def decode_site_key(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> ProteinScopedPhosphositeKey:
    if not isinstance(value, str):
        raise error_type(f"{field_name} must be a string encoded site key")
    token = value.strip()
    if token == "":
        raise error_type(f"{field_name} must be a non-empty encoded site key")

    parts = token.split("|")
    if not parts or parts[0] != _ENCODING_VERSION_PREFIX:
        raise error_type(f"{field_name} must start with '{_ENCODING_VERSION_PREFIX}'")

    values: dict[str, str] = {}
    for part in parts[1:]:
        if "=" not in part:
            raise error_type(
                f"{field_name} contains malformed key-value segment {part!r}"
            )
        key_name, encoded_value = part.split("=", 1)
        if key_name == "":
            raise error_type(f"{field_name} contains an empty key name")
        if key_name not in _ENCODING_FIELD_ORDER:
            raise error_type(
                f"{field_name} contains unsupported key field {key_name!r}"
            )
        if key_name in values:
            raise error_type(f"{field_name} contains duplicate key field {key_name!r}")
        values[key_name] = unquote(encoded_value)

    missing = [name for name in _REQUIRED_FIELDS if name not in values]
    if missing:
        missing_csv = ", ".join(missing)
        raise error_type(
            f"{field_name} is missing required encoded fields: {missing_csv}"
        )

    try:
        parsed_position = int(values["position"])
    except (TypeError, ValueError):
        raise error_type(
            f"{field_name}.position must be an integer in encoded site key"
        ) from None

    return build_protein_scoped_site_key(
        organism=values["organism"],
        protein_namespace=values["protein_namespace"],
        protein_identifier=values["protein_identifier"],
        residue=values["residue"],
        position=parsed_position,
        isoform_id=values.get("isoform_id"),
        field_name=field_name,
        error_type=error_type,
    )


def _required_text(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    if not isinstance(value, str):
        raise error_type(f"{field_name} must be a string")
    token = value.strip()
    if token == "":
        raise error_type(f"{field_name} must be a non-empty string")
    return token


def _optional_text(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise error_type(f"{field_name} must be a string when provided")
    token = value.strip()
    if token == "":
        raise error_type(f"{field_name} must be a non-empty string when provided")
    return token


def _residue_token(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    token = _required_text(value, field_name=field_name, error_type=error_type).upper()
    if token not in _VALID_RESIDUES:
        raise error_type(f"{field_name} must be one of 'S', 'T', or 'Y'")
    return token


def _positive_integer(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise error_type(f"{field_name} must be an integer")
    if value <= 0:
        raise error_type(f"{field_name} must be a positive integer")
    return int(value)


__all__ = [
    "ProteinScopedPhosphositeKey",
    "build_protein_scoped_site_key",
    "encode_site_key",
    "decode_site_key",
]
