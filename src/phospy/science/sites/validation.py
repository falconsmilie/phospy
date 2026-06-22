"""Phosphosite-specific validation helpers."""

from __future__ import annotations

from typing import TypeVar, cast

import pandas as pd

from phospy.science.sites.identifiers import (
    canonicalize_site_index,
    canonicalize_site_series,
    parse_canonical_site_identifier,
)
from phospy.science.sites.site_keys import (
    ProteinScopedPhosphositeKey,
    decode_site_key,
    encode_site_key,
)

ErrorType = TypeVar("ErrorType", bound=Exception)


def require_canonical_site_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[ErrorType],
    strict_supported_format: bool = True,
) -> pd.Index:
    """Require one site index to already use the expected strict format."""

    if not strict_supported_format:
        _require_stripped_site_identifiers(
            index.tolist(),
            field_name=field_name,
            error_type=error_type,
        )
        return index
    canonical = canonicalize_site_index(
        index,
        field_name=field_name,
        error_type=error_type,
        require_unique=True,
        index_name=(str(index.name) if index.name is not None else None),
    )
    if canonical.tolist() != index.tolist():
        raise error_type(
            f"{field_name} must use the recommended site identifier format 'GENE;SITE;'"
        )
    if not index.is_unique:
        duplicate_count, duplicate_labels = _resolve_duplicate_labels(index)
        preview = ", ".join(repr(label) for label in duplicate_labels[:5])
        suffix = "" if len(duplicate_labels) <= 5 else " ..."
        raise error_type(
            f"{field_name} must be unique; duplicate_count={duplicate_count}, "
            f"duplicate_labels={preview}{suffix}"
        )
    return index


def require_canonical_site_series(
    series: pd.Series,
    *,
    field_name: str,
    error_type: type[ErrorType],
    strict_supported_format: bool = True,
) -> pd.Series:
    """Require one site-id series to already use the expected strict format."""

    if not strict_supported_format:
        _require_stripped_site_identifiers(
            series.tolist(),
            field_name=field_name,
            error_type=error_type,
        )
        return series
    canonical = canonicalize_site_series(
        series,
        field_name=field_name,
        error_type=error_type,
    )
    if canonical.tolist() != series.tolist():
        raise error_type(
            f"{field_name} must use the recommended site identifier format 'GENE;SITE;'"
        )
    return series


def require_site_key_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[ErrorType],
    require_unique: bool = True,
) -> pd.Index:
    """Require one index of encoded site_key values."""

    normalised: list[str] = []
    for position, value in enumerate(index.tolist()):
        key = decode_site_key(
            value,
            field_name=f"{field_name}[{position}]",
            error_type=error_type,
        )
        normalised.append(encode_site_key(key))
    if normalised != [str(value) for value in index.tolist()]:
        raise error_type(
            f"{field_name} must contain analysis-ready encoded site_key values"
        )
    if require_unique and not index.is_unique:
        duplicate_count, duplicate_labels = _resolve_duplicate_labels(index)
        preview = ", ".join(repr(label) for label in duplicate_labels[:5])
        suffix = "" if len(duplicate_labels) <= 5 else " ..."
        raise error_type(
            f"{field_name} must be unique; duplicate_count={duplicate_count}, "
            f"duplicate_labels={preview}{suffix}"
        )
    return index


def require_site_key_series(
    series: pd.Series,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> pd.Series:
    """Require one series of analysis-ready encoded site_key values."""

    normalised: list[str] = []
    for row_id, value in series.items():
        key = decode_site_key(
            value,
            field_name=f"{field_name}[{row_id!r}]",
            error_type=error_type,
        )
        normalised.append(encode_site_key(key))
    if normalised != [str(value) for value in series.tolist()]:
        raise error_type(
            f"{field_name} must contain analysis-ready encoded site_key values"
        )
    return series


def require_site_identity_coherence(
    *,
    site_index: pd.Index,
    site_metadata: pd.DataFrame,
    site_index_field_name: str,
    site_metadata_field_name: str,
    gene_symbol_column: str = "gene_symbol",
    site_column: str = "site",
    error_type: type[ErrorType],
    error_preview_limit: int = 5,
) -> None:
    """Require display-site IDs to agree with metadata gene/site columns."""

    unparseable_site_ids: list[str] = []
    mismatched_rows: list[str] = []

    for site_id in site_index:
        parsed = _parse_site_identity(
            site_id,
            field_name=site_index_field_name,
            error_type=error_type,
        )
        if parsed is None:
            unparseable_site_ids.append(str(site_id))
            continue

        expected_gene_symbol, expected_site = parsed
        observed_gene_symbol = site_metadata.at[site_id, gene_symbol_column]
        observed_site = site_metadata.at[site_id, site_column]
        if (
            observed_gene_symbol != expected_gene_symbol
            or observed_site != expected_site
        ):
            mismatched_rows.append(
                f"{site_id} expected(gene_symbol={expected_gene_symbol!r}, "
                f"site={expected_site!r}) observed(gene_symbol={observed_gene_symbol!r}, "
                f"site={observed_site!r})"
            )

    if not unparseable_site_ids and not mismatched_rows:
        return

    details: list[str] = []
    if unparseable_site_ids:
        preview = ", ".join(
            repr(site_id) for site_id in unparseable_site_ids[:error_preview_limit]
        )
        suffix = "" if len(unparseable_site_ids) <= error_preview_limit else " ..."
        details.append(
            f"unparseable site IDs for '<gene_symbol>;<site>;': {preview}{suffix}"
        )
    if mismatched_rows:
        preview = "; ".join(mismatched_rows[:error_preview_limit])
        suffix = "" if len(mismatched_rows) <= error_preview_limit else " ..."
        details.append(f"mismatched rows: {preview}{suffix}")

    joined_details = "; ".join(details)
    raise error_type(
        "dataset site-identity coherence failed: "
        f"{site_index_field_name} display-site IDs must agree with "
        f"{site_metadata_field_name}.{gene_symbol_column} and "
        f"{site_metadata_field_name}.{site_column}; {joined_details}"
    )


def require_no_mixed_site_key_isoform_scope(
    site_keys: pd.Series,
    *,
    field_name: str,
    error_type: type[ErrorType],
    preview_limit: int = 5,
) -> None:
    """Reject mixed missing/specified isoform keys for the same protein-scoped site."""

    decoded_by_row: list[tuple[object, ProteinScopedPhosphositeKey]] = [
        (
            row_id,
            decode_site_key(
                value,
                field_name=f"{field_name}[{row_id!r}]",
                error_type=error_type,
            ),
        )
        for row_id, value in site_keys.items()
    ]

    isoform_scope: dict[tuple[object, ...], set[bool]] = {}
    rows_by_scope: dict[tuple[object, ...], list[object]] = {}
    for row_id, key in decoded_by_row:
        scope = (
            key.organism,
            key.protein_namespace,
            key.protein_identifier,
            key.residue,
            int(key.position),
        )
        isoform_scope.setdefault(scope, set()).add(key.isoform_id is None)
        rows_by_scope.setdefault(scope, []).append(row_id)

    mixed_scopes = [
        (scope, rows_by_scope[scope])
        for scope, flags in isoform_scope.items()
        if len(flags) > 1
    ]
    if not mixed_scopes:
        return

    previews: list[str] = []
    for scope, row_ids in mixed_scopes[:preview_limit]:
        row_preview = ", ".join(repr(str(row_id)) for row_id in row_ids[:preview_limit])
        row_suffix = "" if len(row_ids) <= preview_limit else " ..."
        previews.append(
            "("
            f"organism={scope[0]!r}, "
            f"protein_namespace={scope[1]!r}, "
            f"protein_identifier={scope[2]!r}, "
            f"residue={scope[3]!r}, "
            f"position={scope[4]!r}"
            ") rows=["
            f"{row_preview}{row_suffix}]"
        )
    suffix = "" if len(mixed_scopes) <= preview_limit else " ..."
    raise error_type(
        f"{field_name} contains mixed isoform scope for protein-scoped site "
        "identities; rows for the same "
        "(organism, protein_namespace, protein_identifier, residue, position) "
        "must either all define isoform_id or all omit it; "
        f"conflicts={'; '.join(previews)}{suffix}"
    )


def _parse_site_identity(
    site_id: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> tuple[str, str] | None:
    if not isinstance(site_id, str):
        return None
    try:
        return parse_canonical_site_identifier(
            site_id,
            field_name=field_name,
            error_type=error_type,
        )
    except error_type:
        return None


def _require_stripped_site_identifiers(
    values: list[object],
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> None:
    if any(_is_missing_site_identifier(value) for value in values):
        raise error_type(f"{field_name} must not contain missing site identifiers")
    if not all(isinstance(value, str) for value in values):
        raise error_type(
            f"{field_name} must contain non-empty stripped site identifiers"
        )
    raw_values = [value for value in values if isinstance(value, str)]
    stripped_values = [value.strip() for value in raw_values]
    if any(value == "" for value in stripped_values):
        raise error_type(f"{field_name} must contain non-empty site identifiers")

    collisions: dict[str, set[str]] = {}
    for raw_value, stripped_value in zip(raw_values, stripped_values, strict=False):
        collisions.setdefault(stripped_value, set()).add(raw_value)
    colliding = [value for value, raw_set in collisions.items() if len(raw_set) > 1]
    if colliding:
        preview = ", ".join(colliding[:5])
        suffix = "" if len(colliding) <= 5 else " ..."
        raise error_type(
            f"{field_name} contains colliding site identifiers when stripped: "
            f"{preview}{suffix}"
        )
    if any(
        raw_value != stripped_value
        for raw_value, stripped_value in zip(raw_values, stripped_values, strict=False)
    ):
        raise error_type(
            f"{field_name} must contain non-empty stripped site identifiers"
        )


def _is_missing_site_identifier(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _resolve_duplicate_labels(index: pd.Index) -> tuple[int, list[object]]:
    values = cast(list[object], index.tolist())
    duplicate_count = 0
    duplicate_labels: list[object] = []

    for value in values:
        occurrences = 0
        for candidate in values:
            if candidate == value:
                occurrences += 1
        if occurrences <= 1:
            continue
        duplicate_count += 1
        if any(existing == value for existing in duplicate_labels):
            continue
        duplicate_labels.append(value)

    return duplicate_count, duplicate_labels


__all__ = [
    "require_no_mixed_site_key_isoform_scope",
    "require_canonical_site_index",
    "require_canonical_site_series",
    "require_site_key_index",
    "require_site_key_series",
    "require_site_identity_coherence",
]
