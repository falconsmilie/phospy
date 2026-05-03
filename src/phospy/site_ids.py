"""Shared site-identifier canonicalisation helpers."""

from __future__ import annotations

import re
from typing import TypeVar

import pandas as pd

ErrorType = TypeVar("ErrorType", bound=Exception)
_SITE_IDENTIFIER_PATTERN = re.compile(
    r"^\s*(?P<gene_symbol>[^;]+?)\s*;\s*(?P<site>[^;]+?)\s*;?\s*$"
)
_GENE_TOKEN_PATTERN = re.compile(r"^[^;\s]+$")
_SITE_TOKEN_PATTERN = re.compile(r"^[^;\s]+$")
_SITE_IDENTIFIER_EXPECTATION = (
    "site identifiers must use 'GENE;SITE;' format (example: 'MAPK1;S123;')"
)


def canonicalize_site_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[ErrorType],
    require_unique: bool = True,
    index_name: str | None = None,
) -> pd.Index:
    """Canonicalise one site-identifier index to stripped string labels."""

    canonical = pd.Index(
        [
            canonicalize_site_identifier(
                value,
                field_name=field_name,
                error_type=error_type,
            )
            for value in index.tolist()
        ],
        name=index.name if index_name is None else index_name,
    )
    if require_unique and not canonical.is_unique:
        duplicates = list(dict.fromkeys(canonical[canonical.duplicated(keep=False)]))
        preview = ", ".join(duplicates[:5])
        suffix = "" if len(duplicates) <= 5 else " ..."
        raise error_type(
            f"{field_name} contains duplicate site identifiers after "
            f"canonicalization: {preview}{suffix}"
        )
    return canonical


def canonicalize_site_series(
    series: pd.Series,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> pd.Series:
    """Canonicalise one site-identifier series to stripped string labels."""

    return pd.Series(
        [
            canonicalize_site_identifier(
                value,
                field_name=field_name,
                error_type=error_type,
            )
            for value in series.tolist()
        ],
        index=series.index.copy(),
        name=series.name,
        dtype="object",
    )


def canonicalize_site_identifier(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    """Canonicalise one supported phosphosite identifier to ``GENE;SITE;``."""

    if _is_missing(value):
        raise error_type(f"{field_name} must not contain missing site identifiers")
    raw_identifier = str(value)
    if raw_identifier.strip() == "":
        raise error_type(f"{field_name} must contain non-empty site identifiers")
    match = _SITE_IDENTIFIER_PATTERN.fullmatch(raw_identifier)
    if match is None:
        raise _malformed_site_identifier_error(
            field_name=field_name,
            raw_value=raw_identifier,
            error_type=error_type,
        )
    gene_symbol = match.group("gene_symbol").strip().upper()
    site = match.group("site").strip().upper()
    if not _GENE_TOKEN_PATTERN.fullmatch(gene_symbol):
        raise _malformed_site_identifier_error(
            field_name=field_name,
            raw_value=raw_identifier,
            error_type=error_type,
        )
    if not _SITE_TOKEN_PATTERN.fullmatch(site):
        raise _malformed_site_identifier_error(
            field_name=field_name,
            raw_value=raw_identifier,
            error_type=error_type,
        )
    return f"{gene_symbol};{site};"


def parse_canonical_site_identifier(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> tuple[str, str]:
    """Return canonical ``(gene_symbol, site)`` components for one site ID."""

    canonical = canonicalize_site_identifier(
        value,
        field_name=field_name,
        error_type=error_type,
    )
    gene_symbol, site, _ = canonical.split(";")
    return gene_symbol, site


def canonicalize_site_components(
    gene_symbol: object,
    site: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    """Canonicalise one ``gene_symbol``/``site`` pair to ``GENE;SITE;``."""

    if _is_missing(gene_symbol) or _is_missing(site):
        raise error_type(
            f"{field_name} requires non-missing gene_symbol/site tokens to build "
            "site identifiers"
        )
    gene_symbol_token = str(gene_symbol).strip()
    site_token = str(site).strip()
    if gene_symbol_token == "" or site_token == "":
        raise error_type(
            f"{field_name} requires non-empty gene_symbol/site tokens to build "
            "site identifiers"
        )
    return canonicalize_site_identifier(
        f"{gene_symbol_token};{site_token};",
        field_name=field_name,
        error_type=error_type,
    )


def canonicalize_site_components_series(
    *,
    gene_symbol: pd.Series,
    site: pd.Series,
    field_name: str,
    error_type: type[ErrorType],
    output_name: str = "site_id",
) -> pd.Series:
    """Canonicalise paired ``gene_symbol``/``site`` columns to one site-ID series."""

    if len(gene_symbol.index) != len(site.index):
        raise error_type(
            f"{field_name} requires gene_symbol/site series with equal length"
        )
    canonical = [
        canonicalize_site_components(
            gene_value,
            site_value,
            field_name=field_name,
            error_type=error_type,
        )
        for gene_value, site_value in zip(
            gene_symbol.tolist(),
            site.tolist(),
            strict=False,
        )
    ]
    return pd.Series(
        canonical,
        index=gene_symbol.index.copy(),
        name=output_name,
        dtype="object",
    )


def _malformed_site_identifier_error(
    *,
    field_name: str,
    raw_value: str,
    error_type: type[ErrorType],
    detail: str | None = None,
) -> ErrorType:
    extras = "" if detail is None else f"; {detail}"
    return error_type(
        f"{field_name} {_SITE_IDENTIFIER_EXPECTATION}; got '{raw_value}'{extras}"
    )


def _is_missing(value: object) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False
