from __future__ import annotations

import re
from collections.abc import Iterable

import pandas as pd

from ...errors import TableSchemaError
from ...internal.constants import GENE_P_SITE_COLUMN

_SITE_TOKEN_PATTERN = re.compile(r"^[A-Za-z]+\d+$")
_CANONICAL_SITE_ID_PATTERN = re.compile(
    r"^(?P<entity>[^;\s]+);(?P<site>[A-Za-z]+\d+);$"
)


def normalize_identifier_series(series: pd.Series) -> pd.Series:
    """Normalize identifier values for case/whitespace-insensitive joins."""

    return series.astype("string").str.strip().str.upper()


def require_splitable_gene_p_site(
    series: pd.Series,
    *,
    context: str,
    column_name: str = GENE_P_SITE_COLUMN,
) -> None:
    """Validate gene-site identifiers that must split on a single underscore."""

    normalized = series.astype("string")
    split_columns = normalized.str.split("_", n=1, expand=True)
    underscore_count = normalized.str.count("_")
    if split_columns.shape[1] < 2:
        invalid_mask = pd.Series(True, index=series.index)
    else:
        site_tokens = split_columns[1].astype("string").str.strip()
        invalid_mask = (
            normalized.isna()
            | (underscore_count != 1)
            | split_columns[0].isna()
            | split_columns[1].isna()
            | (split_columns[0].str.strip().str.len() == 0)
            | (split_columns[1].str.strip().str.len() == 0)
            | (~site_tokens.str.fullmatch(_SITE_TOKEN_PATTERN.pattern).fillna(False))
        )
    if invalid_mask.any():
        sample_values = series.loc[invalid_mask].astype(str).unique()[:3]
        sample_preview = ", ".join(str(value) for value in sample_values)
        msg = (
            f"{context} contains malformed {column_name} values that cannot be split "
            "into non-empty gene and site parts using a single underscore "
            "with a site token like 'S123': "
            f"{sample_preview}"
        )
        raise TableSchemaError(msg)


def parse_canonical_site_id(site_id: object) -> tuple[str, str] | None:
    """Parse a canonical phosphosite ID in ``ENTITY;SITE;`` format."""

    match = _CANONICAL_SITE_ID_PATTERN.fullmatch(str(site_id).strip())
    if match is None:
        return None
    entity = match.group("entity").strip().upper()
    site_token = match.group("site").strip().upper()
    return entity, site_token


def build_canonical_site_id(
    *,
    entity: object,
    site_token: object,
    context: str,
) -> str:
    """Build a canonical phosphosite ID in ``ENTITY;SITE;`` format."""

    normalized_entity = str(entity).strip().upper()
    normalized_site_token = str(site_token).strip().upper()
    candidate = f"{normalized_entity};{normalized_site_token};"
    if parse_canonical_site_id(candidate) is None:
        msg = (
            f"{context} must produce canonical site IDs in 'ENTITY;SITE;' format "
            f"with a site token like 'S123'. Invalid value: {candidate}"
        )
        raise TableSchemaError(msg)
    return candidate


def require_canonical_site_ids(
    site_ids: Iterable[object],
    *,
    context: str,
    label: str = "site IDs",
) -> None:
    """Validate canonical phosphosite identifiers in ``ENTITY;SITE;`` format."""

    invalid_values = [
        str(site_id) for site_id in site_ids if parse_canonical_site_id(site_id) is None
    ]
    if invalid_values:
        preview = ", ".join(invalid_values[:3])
        msg = (
            f"{context} must contain canonical {label} in 'ENTITY;SITE;' format "
            f"with a site token like 'S123'. Invalid values: {preview}"
        )
        if len(invalid_values) > 3:
            msg += ", ..."
        raise TableSchemaError(msg)


__all__ = [
    "build_canonical_site_id",
    "normalize_identifier_series",
    "parse_canonical_site_id",
    "require_canonical_site_ids",
    "require_splitable_gene_p_site",
]
