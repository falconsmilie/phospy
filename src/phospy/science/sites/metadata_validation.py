"""Science-owned phosphosite metadata validation helpers."""

from __future__ import annotations

from typing import TypeVar

import pandas as pd

from phospy.science.sites.identifiers import ParsedSiteToken, try_parse_site_token
from phospy.science.sites.identity import (
    build_phosphosite_identity,
    validate_identity_optional_columns,
)
from phospy.science.sites.site_keys import require_positive_integer_position

ErrorType = TypeVar("ErrorType", bound=Exception)
_EXAMPLE_LIMIT = 5
_PHOSPHORYLATABLE_RESIDUES = frozenset({"S", "T", "Y"})
_CANONICAL_AMINO_ACID_RESIDUES = frozenset(
    {
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    }
)
_SUPPORTED_UNKNOWN_SEQUENCE_CHARACTERS = frozenset({"X"})
_SUPPORTED_BASE_SEQUENCE_SPECIAL_CHARACTERS = frozenset({"_", "-"})
_SUPPORTED_BASE_SEQUENCE_CHARACTERS = frozenset(
    _CANONICAL_AMINO_ACID_RESIDUES
    | _SUPPORTED_UNKNOWN_SEQUENCE_CHARACTERS
    | _SUPPORTED_BASE_SEQUENCE_SPECIAL_CHARACTERS
)
_SITE_POSITION_CANDIDATE_COLUMNS = ("site_position", "position")


def validate_site_identity_metadata(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_column: str = "site",
    site_sequence_column: str = "site_sequence",
    residue_column: str = "residue",
    allow_opaque_site_values: bool = False,
) -> None:
    """Validate row-level phosphosite identity coherence metadata."""

    validate_identity_optional_columns(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
    )
    malformed_site_values: list[str] = []
    inconsistent_residue_rows: list[str] = []
    inconsistent_position_rows: list[str] = []
    non_phospho_centre_rows: list[str] = []
    centre_mismatch_rows: list[str] = []
    invalid_explicit_residue_rows: list[str] = []
    invalid_explicit_position_rows: list[str] = []

    site_positions = _resolve_site_position_series(site_metadata)

    for site_id in site_metadata.index.tolist():
        site_value = site_metadata.at[site_id, site_column]
        parsed_site = try_parse_site_token(site_value)
        if parsed_site is None and not allow_opaque_site_values:
            malformed_site_values.append(f"{site_id!r}:{site_value!r}")

        explicit_residue = _resolve_optional_residue(
            site_metadata.at[site_id, residue_column]
            if residue_column in site_metadata.columns
            else None
        )
        if residue_column in site_metadata.columns and explicit_residue is None:
            raw_residue = site_metadata.at[site_id, residue_column]
            if not _is_missing(raw_residue):
                invalid_explicit_residue_rows.append(f"{site_id!r}:{raw_residue!r}")

        explicit_position: int | None = None
        if site_positions.name is not None:
            raw_position = site_positions.at[site_id]
            try:
                explicit_position = require_positive_integer_position(
                    raw_position,
                    field_name=f"{field_name}.{site_positions.name}",
                    error_type=error_type,
                )
            except error_type:
                invalid_explicit_position_rows.append(f"{site_id!r}:{raw_position!r}")

        if parsed_site is not None and explicit_residue is not None:
            if explicit_residue != parsed_site.residue:
                inconsistent_residue_rows.append(
                    f"{site_id!r}: site={parsed_site.residue!r}, "
                    f"residue_column={explicit_residue!r}"
                )
        if parsed_site is not None and explicit_position is not None:
            if explicit_position != parsed_site.position:
                inconsistent_position_rows.append(
                    f"{site_id!r}: site={parsed_site.position}, "
                    f"site_position_column={explicit_position}"
                )

        if site_sequence_column not in site_metadata.columns:
            continue
        parsed_sequence = _resolve_optional_sequence(
            site_metadata.at[site_id, site_sequence_column]
        )
        if parsed_sequence is None:
            continue
        if not _sequence_supports_central_residue_check(parsed_sequence):
            continue
        centre_residue = _resolve_central_residue(parsed_sequence)
        if centre_residue is None:
            continue
        if centre_residue not in _PHOSPHORYLATABLE_RESIDUES:
            non_phospho_centre_rows.append(
                f"{site_id!r}: centre={centre_residue!r}, sequence={parsed_sequence!r}"
            )
            continue
        expected_residue = _resolve_expected_residue(parsed_site, explicit_residue)
        if expected_residue is not None and centre_residue != expected_residue:
            centre_mismatch_rows.append(
                f"{site_id!r}: expected={expected_residue!r}, "
                f"observed={centre_residue!r}"
            )

    validate_site_sequence_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=site_sequence_column,
    )

    details: list[str] = []
    if malformed_site_values:
        details.append(
            "site values must use strict 'S/T/Y<position>' tokens (example: "
            "'S123') unless opaque-site mode is explicitly enabled; "
            + _summarise_examples(malformed_site_values)
        )
    if invalid_explicit_residue_rows:
        details.append(
            f"{field_name}.residue must be one residue letter when provided; "
            + _summarise_examples(invalid_explicit_residue_rows)
        )
    if invalid_explicit_position_rows:
        details.append(
            f"{field_name}.{site_positions.name} must contain positive integer "
            "values when the column is present; "
            + _summarise_examples(invalid_explicit_position_rows)
        )
    if inconsistent_residue_rows:
        details.append(
            "residue column must match parsed site residue when both are present; "
            + _summarise_examples(inconsistent_residue_rows)
        )
    if inconsistent_position_rows:
        details.append(
            "site position column must match parsed site position when both are "
            "present; " + _summarise_examples(inconsistent_position_rows)
        )
    if non_phospho_centre_rows:
        details.append(
            f"{field_name}.site_sequence must contain a centred phosphorylatable "
            "residue (S/T/Y); " + _summarise_examples(non_phospho_centre_rows)
        )
    if centre_mismatch_rows:
        details.append(
            "site_sequence central residue must agree with site/residue metadata; "
            + _summarise_examples(centre_mismatch_rows)
        )

    if details:
        raise error_type(
            f"{field_name} phosphosite identity metadata validation failed; "
            + "; ".join(details)
        )


def validate_site_sequence_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "site_sequence",
) -> None:
    """Validate site-sequence strings as plausible amino-acid contexts."""

    if column_name not in site_metadata.columns:
        return
    invalid_rows: list[str] = []
    values = site_metadata[column_name]
    for site_id, raw_value in values.items():
        if not isinstance(raw_value, str):
            invalid_rows.append(f"{site_id!r}:{raw_value!r}")
            continue
        sequence = raw_value.strip().upper()
        if sequence == "":
            invalid_rows.append(f"{site_id!r}:{raw_value!r}:blank_sequence")
            continue
        if len(sequence) < 3:
            invalid_rows.append(f"{site_id!r}:{raw_value!r}:sequence_too_short")
            continue
        unsupported_characters = sorted(
            character
            for character in sequence
            if character not in _SUPPORTED_BASE_SEQUENCE_CHARACTERS
        )
        if unsupported_characters:
            invalid_rows.append(
                f"{site_id!r}:{raw_value!r}:"
                f"unsupported_characters={''.join(unsupported_characters)!r}"
            )
            continue
        if not any(
            character in _CANONICAL_AMINO_ACID_RESIDUES for character in sequence
        ):
            invalid_rows.append(f"{site_id!r}:{raw_value!r}:no_residue_letters")
            continue
    if invalid_rows:
        raise error_type(
            f"{field_name}.{column_name} must be plausible amino-acid context "
            "strings (allowed residues: ACDEFGHIKLMNPQRSTVWY; allowed unknown: X; "
            "allowed gap placeholders: '_' and '-'); "
            + _summarise_examples(invalid_rows)
        )


def enforce_site_identity_rows(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    allow_opaque_site_values: bool = False,
) -> None:
    """Enforce row-level phosphosite identity parsing for table boundaries."""

    for row_position, site_id in enumerate(site_metadata.index.tolist()):
        row = site_metadata.iloc[row_position]  # pyright: ignore[reportUnknownMemberType]
        _ = build_phosphosite_identity(
            display_id=site_id,
            gene_symbol=row["gene_symbol"],
            site=row["site"],
            allow_opaque_site_values=allow_opaque_site_values,
            protein_id=(
                None if "protein_id" not in site_metadata.columns else row["protein_id"]
            ),
            protein_accession=(
                None
                if "protein_accession" not in site_metadata.columns
                else row["protein_accession"]
            ),
            field_name=f"{field_name}[{row_position}:{site_id!r}]",
            error_type=error_type,
        )


def _resolve_site_position_series(site_metadata: pd.DataFrame) -> pd.Series:
    for column_name in _SITE_POSITION_CANDIDATE_COLUMNS:
        if column_name in site_metadata.columns:
            return pd.Series(site_metadata[column_name], index=site_metadata.index)
    return pd.Series(index=site_metadata.index, dtype="object")


def _resolve_optional_residue(value: object) -> str | None:
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        return None
    token = value.strip().upper()
    if len(token) != 1:
        return None
    if not token.isalpha():
        return None
    return token


def _resolve_optional_sequence(value: object) -> str | None:
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        return None
    stripped = value.strip().upper()
    if stripped == "":
        return None
    return stripped


def _resolve_central_residue(site_sequence: str) -> str | None:
    sequence_length = len(site_sequence)
    if sequence_length == 0 or sequence_length % 2 == 0:
        return None
    return site_sequence[sequence_length // 2]


def _sequence_supports_central_residue_check(site_sequence: str) -> bool:
    if not site_sequence.isalpha():
        return False
    return len(site_sequence) % 2 == 1


def _resolve_expected_residue(
    parsed_site: ParsedSiteToken | None,
    explicit_residue: str | None,
) -> str | None:
    if explicit_residue is not None:
        return explicit_residue
    if parsed_site is not None:
        return parsed_site.residue
    return None


def _is_missing(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _summarise_examples(values: list[str], *, limit: int = _EXAMPLE_LIMIT) -> str:
    suffix = "" if len(values) <= limit else " ..."
    return "[" + ", ".join(values[:limit]) + suffix + "]"


__all__ = [
    "enforce_site_identity_rows",
    "validate_site_identity_metadata",
    "validate_site_sequence_column",
]
