"""Phosphosite identity domain model and collision validation helpers."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import TypeVar

import pandas as pd

from phospy.sites.identifiers import (
    canonicalize_site_components,
    canonicalize_site_identifier,
    parse_canonical_site_identifier,
    try_parse_site_token,
)

ErrorType = TypeVar("ErrorType", bound=Exception)

PHOSPHOSITE_IDENTITY_OPTIONAL_COLUMNS = (
    "organism",
    "protein_id",
    "protein_accession",
    "isoform_id",
    "source_namespace",
    "source_site_id",
)
PHOSPHOSITE_PROTEIN_CONTEXT_COLUMNS = ("protein_id", "protein_accession")


@dataclass(frozen=True, slots=True)
class PhosphositeIdentity:
    """Scientific phosphosite identity with display-ID compatibility fields."""

    display_id: str
    gene_symbol: str
    residue: str | None
    position: int | None
    organism: str | None = None
    protein_id: str | None = None
    protein_accession: str | None = None
    isoform_id: str | None = None
    source_namespace: str | None = None
    source_site_id: str | None = None

    def scientific_context_key(self) -> tuple[object, ...]:
        return (
            self.gene_symbol,
            self.residue,
            self.position,
            self.organism,
            self.protein_id,
            self.protein_accession,
            self.isoform_id,
            self.source_namespace,
            self.source_site_id,
        )

    def has_protein_context(self) -> bool:
        return bool(_has_text(self.protein_id) or _has_text(self.protein_accession))


def build_phosphosite_identity(
    *,
    display_id: object,
    gene_symbol: object,
    site: object,
    organism: object = None,
    protein_id: object = None,
    protein_accession: object = None,
    isoform_id: object = None,
    source_namespace: object = None,
    source_site_id: object = None,
    allow_opaque_site_values: bool = False,
    field_name: str,
    error_type: type[ErrorType],
) -> PhosphositeIdentity:
    """Build one structured phosphosite identity from row-level metadata."""

    canonical_display_id = canonicalize_site_identifier(
        display_id,
        field_name=f"{field_name}.display_id",
        error_type=error_type,
    )
    canonical_components_id = canonicalize_site_components(
        gene_symbol,
        site,
        field_name=f"{field_name}.gene_symbol/site",
        error_type=error_type,
    )
    if canonical_components_id != canonical_display_id:
        raise error_type(
            f"{field_name} has inconsistent display and metadata identity fields; "
            f"display_id={canonical_display_id!r}, "
            f"gene_symbol/site={canonical_components_id!r}"
        )
    canonical_gene_symbol, canonical_site = parse_canonical_site_identifier(
        canonical_display_id,
        field_name=f"{field_name}.display_id",
        error_type=error_type,
    )
    parsed_site = try_parse_site_token(canonical_site)
    if parsed_site is None and not allow_opaque_site_values:
        raise error_type(
            f"{field_name}.site must use '<residue><position>' tokens to build a "
            f"scientific phosphosite identity; got {canonical_site!r}"
        )

    return PhosphositeIdentity(
        display_id=canonical_display_id,
        gene_symbol=canonical_gene_symbol,
        residue=(None if parsed_site is None else parsed_site.residue),
        position=(None if parsed_site is None else int(parsed_site.position)),
        organism=_optional_text(
            organism,
            field_name=f"{field_name}.organism",
            error_type=error_type,
        ),
        protein_id=_optional_text(
            protein_id,
            field_name=f"{field_name}.protein_id",
            error_type=error_type,
        ),
        protein_accession=_optional_text(
            protein_accession,
            field_name=f"{field_name}.protein_accession",
            error_type=error_type,
        ),
        isoform_id=_optional_text(
            isoform_id,
            field_name=f"{field_name}.isoform_id",
            error_type=error_type,
        ),
        source_namespace=_optional_text(
            source_namespace,
            field_name=f"{field_name}.source_namespace",
            error_type=error_type,
        ),
        source_site_id=_optional_text(
            source_site_id,
            field_name=f"{field_name}.source_site_id",
            error_type=error_type,
        ),
    )


def validate_identity_optional_columns(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    optional_columns: Iterable[str] = PHOSPHOSITE_IDENTITY_OPTIONAL_COLUMNS,
) -> None:
    """Validate optional identity columns as missing-or-canonical-string fields."""

    invalid_rows: list[str] = []
    for column_name in optional_columns:
        if column_name not in site_metadata.columns:
            continue
        series = site_metadata.loc[:, column_name]
        for site_id, raw_value in series.items():
            try:
                _optional_text(
                    raw_value,
                    field_name=f"{field_name}.{column_name}",
                    error_type=error_type,
                )
            except error_type:
                invalid_rows.append(f"{site_id!r}:{column_name!r}:{raw_value!r}")
    if invalid_rows:
        preview = ", ".join(invalid_rows[:5])
        suffix = "" if len(invalid_rows) <= 5 else " ..."
        raise error_type(
            f"{field_name} optional identity columns must contain missing values or "
            f"canonical non-empty strings; invalid_rows={preview}{suffix}"
        )


def validate_no_conflicting_identity_collisions(
    *,
    site_metadata: pd.DataFrame,
    display_ids: pd.Series,
    field_name: str,
    error_type: type[ErrorType],
    preview_limit: int = 5,
) -> None:
    """Reject ambiguous same-display-ID rows with conflicting scientific identity."""

    if len(site_metadata.index) != len(display_ids.index):
        raise error_type(
            f"{field_name} requires aligned site_metadata/display_ids lengths; "
            f"site_metadata_rows={int(len(site_metadata.index))}, "
            f"display_id_rows={int(len(display_ids.index))}"
        )

    if site_metadata.empty:
        return

    display_key_to_contexts: dict[str, set[tuple[object, ...]]] = {}
    display_key_to_rows: dict[str, list[str]] = {}

    for position, row_id in enumerate(site_metadata.index.tolist()):
        display_id = display_ids.iloc[position]
        row = site_metadata.iloc[position]
        identity = build_phosphosite_identity(
            display_id=display_id,
            gene_symbol=row["gene_symbol"],
            site=row["site"],
            organism=_series_value_or_none(
                site_metadata,
                row_position=position,
                column="organism",
            ),
            protein_id=_series_value_or_none(
                site_metadata, row_position=position, column="protein_id"
            ),
            protein_accession=_series_value_or_none(
                site_metadata, row_position=position, column="protein_accession"
            ),
            isoform_id=_series_value_or_none(
                site_metadata, row_position=position, column="isoform_id"
            ),
            source_namespace=_series_value_or_none(
                site_metadata, row_position=position, column="source_namespace"
            ),
            source_site_id=_series_value_or_none(
                site_metadata, row_position=position, column="source_site_id"
            ),
            field_name=f"{field_name}[{row_id!r}]",
            error_type=error_type,
        )
        display_key_to_contexts.setdefault(identity.display_id, set()).add(
            identity.scientific_context_key()
        )
        display_key_to_rows.setdefault(identity.display_id, []).append(str(row_id))

    conflicts: list[str] = []
    for display_id, contexts in display_key_to_contexts.items():
        if len(contexts) <= 1:
            continue
        rows_for_display = display_key_to_rows.get(display_id, ())
        if not _requires_strict_collision_rejection(
            site_metadata=site_metadata,
            row_ids=rows_for_display,
        ):
            continue
        row_ids = display_key_to_rows.get(display_id, ())
        row_preview = ", ".join(repr(value) for value in tuple(row_ids)[:preview_limit])
        row_suffix = "" if len(row_ids) <= preview_limit else " ..."
        conflicts.append(
            f"{display_id!r} has {int(len(contexts))} conflicting scientific "
            f"identities across source rows [{row_preview}{row_suffix}]"
        )

    if conflicts:
        preview = "; ".join(conflicts[:preview_limit])
        suffix = "" if len(conflicts) <= preview_limit else " ..."
        raise error_type(
            f"{field_name} contains conflicting scientific identities for duplicate "
            f"display site IDs; {preview}{suffix}. Provide explicit disambiguation "
            "upstream or remove ambiguous duplicates."
        )


def _series_value_or_none(
    frame: pd.DataFrame,
    *,
    row_position: int,
    column: str,
) -> object:
    if column not in frame.columns:
        return None
    return frame.iloc[row_position][column]


def _optional_text(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> str | None:
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        raise error_type(f"{field_name} must be a string when provided")
    token = value.strip()
    if token == "":
        return None
    return token


def _has_text(value: str | None) -> bool:
    return value is not None and value.strip() != ""


def _requires_strict_collision_rejection(
    *,
    site_metadata: pd.DataFrame,
    row_ids: Iterable[str],
) -> bool:
    strict_columns = (
        "protein_accession",
        "isoform_id",
        "source_namespace",
        "source_site_id",
    )
    row_index = set(str(row_id) for row_id in row_ids)
    if not row_index:
        return False
    for column in strict_columns:
        if column not in site_metadata.columns:
            continue
        for row_id in site_metadata.index.tolist():
            if str(row_id) not in row_index:
                continue
            value = site_metadata.at[row_id, column]
            if _has_strict_identity_value(value):
                return True
    return False


def _is_missing(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _has_strict_identity_value(value: object) -> bool:
    if _is_missing(value):
        return False
    if isinstance(value, str):
        return value.strip() != ""
    return True


__all__ = [
    "PHOSPHOSITE_IDENTITY_OPTIONAL_COLUMNS",
    "PHOSPHOSITE_PROTEIN_CONTEXT_COLUMNS",
    "PhosphositeIdentity",
    "build_phosphosite_identity",
    "validate_identity_optional_columns",
    "validate_no_conflicting_identity_collisions",
]
