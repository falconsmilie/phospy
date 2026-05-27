"""Dataset-boundary validation for protein-scoped phosphosite identity."""

from __future__ import annotations

import numbers
from typing import Any, TypeVar, cast

import pandas as pd

from phospy.science.sites.identifiers import (
    ParsedSiteToken,
    canonicalize_site_series,
    try_parse_site_token,
)
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    decode_site_key,
    encode_site_key,
)
from phospy.science.sites.validation import require_no_mixed_site_key_isoform_scope

ErrorType = TypeVar("ErrorType", bound=Exception)
_SITE_POSITION_CANDIDATE_COLUMNS = ("position", "site_position")


def enforce_site_key_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "site_key",
) -> pd.Series:
    """Require one present site_key column with decodable encoded key values."""

    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    values = pd.Series(site_metadata[column_name], dtype="object")
    normalized: list[str] = []
    for row_id, raw_value in values.items():
        key = decode_site_key(
            raw_value,
            field_name=f"{field_name}.{column_name}[{row_id!r}]",
            error_type=error_type,
        )
        normalized.append(encode_site_key(key))
    return pd.Series(
        normalized,
        index=pd.Index(site_metadata.index),
        name=column_name,
        dtype="object",
    )


def enforce_display_id_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "display_id",
) -> pd.Series:
    """Require one present display_id column with canonical site identifiers."""

    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    return canonicalize_site_series(
        pd.Series(site_metadata[column_name], dtype="object"),
        field_name=f"{field_name}.{column_name}",
        error_type=error_type,
    )


def enforce_unique_site_key_identity(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str = "site_key",
    preview_limit: int = 5,
) -> pd.Series:
    """Require unique site_key identity and strict mixed-isoform consistency."""

    site_keys = enforce_site_key_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=site_key_column,
    )
    duplicate_mask = site_keys.duplicated(keep=False)
    if bool(duplicate_mask.any()):
        duplicate_values = list(
            dict.fromkeys(site_keys.loc[duplicate_mask].astype(str).tolist())
        )
        preview = ", ".join(repr(value) for value in duplicate_values[:preview_limit])
        suffix = "" if len(duplicate_values) <= preview_limit else " ..."
        raise error_type(
            f"{field_name}.{site_key_column} must be unique; "
            f"duplicate_values=[{preview}{suffix}]"
        )

    require_no_mixed_site_key_isoform_scope(
        site_keys=site_keys,
        field_name=f"{field_name}.{site_key_column}",
        error_type=error_type,
        preview_limit=preview_limit,
    )
    return site_keys


def enforce_site_key_matches_metadata(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str = "site_key",
    isoform_column: str = "isoform_id",
    preview_limit: int = 5,
) -> pd.Series:
    """Require encoded site_key values to match metadata-derived protein keys."""

    site_keys = enforce_site_key_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=site_key_column,
    )

    mismatches: list[str] = []
    for row_id, encoded_site_key in site_keys.items():
        decoded_key = decode_site_key(
            encoded_site_key,
            field_name=f"{field_name}.{site_key_column}[{row_id!r}]",
            error_type=error_type,
        )
        metadata_key = build_protein_scoped_site_key(
            organism=_required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="organism",
                field_name=field_name,
                error_type=error_type,
            ),
            protein_namespace=_required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="protein_namespace",
                field_name=field_name,
                error_type=error_type,
            ),
            protein_identifier=_required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="protein_identifier",
                field_name=field_name,
                error_type=error_type,
            ),
            residue=_resolve_row_residue(
                site_metadata=site_metadata,
                row_id=row_id,
                field_name=field_name,
                error_type=error_type,
            ),
            position=_resolve_row_position(
                site_metadata=site_metadata,
                row_id=row_id,
                field_name=field_name,
                error_type=error_type,
            ),
            isoform_id=_optional_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name=isoform_column,
                field_name=field_name,
                error_type=error_type,
            ),
            field_name=f"{field_name}[{row_id!r}]",
            error_type=error_type,
        )
        if decoded_key != metadata_key:
            mismatches.append(
                f"{row_id!r}:observed={encoded_site_key!r}:"
                f"expected={encode_site_key(metadata_key)!r}"
            )

    if mismatches:
        preview = ", ".join(mismatches[:preview_limit])
        suffix = "" if len(mismatches) <= preview_limit else " ..."
        raise error_type(
            f"{field_name}.{site_key_column} must match metadata-derived "
            f"ProteinScopedPhosphositeKey values; mismatches=[{preview}{suffix}]"
        )
    return site_keys


def enforce_site_key_index(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str = "site_key",
    preview_limit: int = 5,
) -> None:
    """Require site_metadata index labels to match site_key values."""

    site_keys = enforce_site_key_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=site_key_column,
    )
    mismatches: list[str] = []
    for row_id, encoded_site_key in site_keys.items():
        index_value = row_id
        if not isinstance(index_value, str):
            mismatches.append(f"{row_id!r}:index={index_value!r}")
            continue
        if index_value != encoded_site_key:
            mismatches.append(
                f"{row_id!r}:index={index_value!r}:site_key={encoded_site_key!r}"
            )
    if not mismatches:
        return
    preview = ", ".join(mismatches[:preview_limit])
    suffix = "" if len(mismatches) <= preview_limit else " ..."
    raise error_type(
        f"{field_name}.index must match {field_name}.{site_key_column} when "
        f"enforced; mismatches=[{preview}{suffix}]"
    )


def _required_text_value(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    column_name: str,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    value = site_metadata.at[row_id, column_name]
    if not isinstance(value, str):
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a non-empty string"
        )
    token = value.strip()
    if token == "":
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a non-empty string"
        )
    return token


def _optional_text_value(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    column_name: str,
    field_name: str,
    error_type: type[ErrorType],
) -> str | None:
    if column_name not in site_metadata.columns:
        return None
    value = site_metadata.at[row_id, column_name]
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a string when provided"
        )
    token = value.strip()
    if token == "":
        return None
    return token


def _resolve_row_residue(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    explicit_residue = (
        _optional_text_value(
            site_metadata=site_metadata,
            row_id=row_id,
            column_name="residue",
            field_name=field_name,
            error_type=error_type,
        )
        if "residue" in site_metadata.columns
        else None
    )
    parsed_site = _parse_row_site_token(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    if explicit_residue is None:
        if parsed_site is None:
            raise error_type(
                f"{field_name}[{row_id!r}] requires residue metadata or strict site "
                "token parsing to derive ProteinScopedPhosphositeKey"
            )
        return parsed_site.residue
    token = explicit_residue.upper()
    if len(token) != 1 or token not in {"S", "T", "Y"}:
        raise error_type(
            f"{field_name}[{row_id!r}].residue must be one of 'S', 'T', or 'Y'"
        )
    if parsed_site is not None and parsed_site.residue != token:
        raise error_type(
            f"{field_name}[{row_id!r}] has inconsistent residue metadata; "
            f"site_token={parsed_site.residue!r}, residue_column={token!r}"
        )
    return token


def _resolve_row_position(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> int:
    explicit_position = _resolve_explicit_position(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    parsed_site = _parse_row_site_token(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    if explicit_position is None:
        if parsed_site is None:
            raise error_type(
                f"{field_name}[{row_id!r}] requires position metadata or strict site "
                "token parsing to derive ProteinScopedPhosphositeKey"
            )
        return int(parsed_site.position)
    if parsed_site is not None and explicit_position != int(parsed_site.position):
        raise error_type(
            f"{field_name}[{row_id!r}] has inconsistent position metadata; "
            f"site_token={parsed_site.position}, position_column={explicit_position}"
        )
    return explicit_position


def _resolve_explicit_position(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> int | None:
    for column_name in _SITE_POSITION_CANDIDATE_COLUMNS:
        if column_name not in site_metadata.columns:
            continue
        raw_value = site_metadata.at[row_id, column_name]
        if _is_missing(raw_value):
            continue
        if isinstance(raw_value, bool):
            raise error_type(
                f"{field_name}[{row_id!r}].{column_name} must be an integer >= 1"
            )
        if isinstance(raw_value, numbers.Integral):
            integer_value = int(raw_value)
            if integer_value < 1:
                raise error_type(
                    f"{field_name}[{row_id!r}].{column_name} must be an integer >= 1"
                )
            return integer_value
        if isinstance(raw_value, numbers.Real):
            numeric_value = float(cast(Any, raw_value))
            if not numeric_value.is_integer() or numeric_value < 1:
                raise error_type(
                    f"{field_name}[{row_id!r}].{column_name} must be an integer >= 1"
                )
            return int(numeric_value)
        if isinstance(raw_value, str):
            stripped = raw_value.strip()
            if stripped == "":
                continue
            try:
                parsed = int(stripped)
            except ValueError as exc:
                raise error_type(
                    f"{field_name}[{row_id!r}].{column_name} must be an integer >= 1"
                ) from exc
            if parsed < 1:
                raise error_type(
                    f"{field_name}[{row_id!r}].{column_name} must be an integer >= 1"
                )
            return int(parsed)
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be an integer >= 1"
        )
    return None


def _parse_row_site_token(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> ParsedSiteToken | None:
    if "site" not in site_metadata.columns:
        return None
    parsed = try_parse_site_token(site_metadata.at[row_id, "site"])
    if parsed is not None:
        return parsed
    raw_value = site_metadata.at[row_id, "site"]
    if _is_missing(raw_value):
        return None
    if isinstance(raw_value, str) and raw_value.strip() == "":
        return None
    raise error_type(
        f"{field_name}[{row_id!r}].site must use strict 'S/T/Y<position>' tokens"
    )


def _is_missing(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


__all__ = [
    "enforce_display_id_column",
    "enforce_site_key_column",
    "enforce_site_key_index",
    "enforce_site_key_matches_metadata",
    "enforce_unique_site_key_identity",
]
