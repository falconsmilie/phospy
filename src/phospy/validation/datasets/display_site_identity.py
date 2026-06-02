"""Legacy validation for phosphosite display-label duplicate guards."""

from __future__ import annotations

from typing import TypeVar

import pandas as pd

from phospy.science.sites.identifiers import canonicalize_site_series

ErrorType = TypeVar("ErrorType", bound=Exception)

DISPLAY_SITE_CONTEXT_COLUMNS = (
    "protein_id",
    "protein_accession",
    "protein_group",
    "isoform_id",
    "source",
    "source_namespace",
    "site_id",
    "source_site_id",
)


def enforce_unique_display_site_identity_rows(
    *,
    site_metadata: pd.DataFrame,
    display_site_ids: pd.Series,
    field_name: str,
    error_type: type[ErrorType],
    context_columns: tuple[str, ...] = DISPLAY_SITE_CONTEXT_COLUMNS,
    preview_limit: int = 5,
) -> pd.Series:
    """Reject duplicate display labels for callers using the legacy guard.

    Current analysis-ready row identity is ``site_key``. This helper is retained
    for compatibility with legacy duplicate checks and must not be treated as
    the current dataset row-identity boundary.
    """

    if len(site_metadata.index) != len(display_site_ids.index):
        raise error_type(
            f"{field_name} requires aligned site_metadata/display_site_ids lengths; "
            f"site_metadata_rows={int(len(site_metadata.index))}, "
            f"display_site_id_rows={int(len(display_site_ids.index))}"
        )
    if site_metadata.empty:
        return display_site_ids.copy()

    normalised_display_ids = canonicalize_site_series(
        display_site_ids,
        field_name=f"{field_name}.display_site_ids",
        error_type=error_type,
    )
    duplicate_mask = normalised_display_ids.duplicated(keep=False)
    if not bool(duplicate_mask.any()):
        return normalised_display_ids

    duplicate_ids = list(
        dict.fromkeys(normalised_display_ids.loc[duplicate_mask].astype(str).tolist())
    )
    conflicting_messages: list[str] = []
    plain_messages: list[str] = []

    normalized_labels = normalised_display_ids.astype(str).tolist()
    for duplicate_id in duplicate_ids:
        duplicate_positions = tuple(
            position
            for position, normalized_label in enumerate(normalized_labels)
            if normalized_label == duplicate_id
        )
        conflicts = _describe_context_conflicts(
            site_metadata=site_metadata,
            duplicate_positions=duplicate_positions,
            context_columns=context_columns,
            preview_limit=preview_limit,
        )
        if conflicts:
            conflicting_messages.append(
                "Duplicate phosphosite display identifier "
                f"{duplicate_id!r} maps to multiple protein or isoform contexts"
                f" ({conflicts}). This legacy display-label duplicate guard "
                "cannot derive a safe site_key for those rows. Resolve or "
                "disambiguate the protein context before using this legacy path."
            )
            continue
        plain_messages.append(
            "Duplicate phosphosite display identifier "
            f"{duplicate_id!r} appears more than once in this legacy display-label "
            "duplicate guard. Aggregate, remove, or derive distinct site_key rows "
            "before using this legacy path."
        )

    duplicate_preview = ", ".join(
        repr(site_id) for site_id in duplicate_ids[:preview_limit]
    )
    duplicate_suffix = "" if len(duplicate_ids) <= preview_limit else " ..."
    summary = (
        f"{field_name} contains duplicate normalised phosphosite display "
        f"identifiers: {duplicate_preview}{duplicate_suffix}"
    )
    details = conflicting_messages + plain_messages
    detail_text = " ".join(details[: max(preview_limit, len(details))])
    raise error_type(f"{summary}. {detail_text}")


def _describe_context_conflicts(
    *,
    site_metadata: pd.DataFrame,
    duplicate_positions: tuple[int, ...],
    context_columns: tuple[str, ...],
    preview_limit: int,
) -> str:
    conflicts: list[str] = []
    for column_name in context_columns:
        if column_name not in site_metadata.columns:
            continue
        column = site_metadata[column_name]
        column_values = column.tolist()
        values = {
            _normalise_context_value(column_values[position])
            for position in duplicate_positions
        }
        if len(values) <= 1:
            continue
        ordered_values = sorted(values, key=repr)
        value_preview = ", ".join(
            repr(value) for value in ordered_values[:preview_limit]
        )
        value_suffix = "" if len(ordered_values) <= preview_limit else " ..."
        conflicts.append(f"{column_name}=[{value_preview}{value_suffix}]")
    return ", ".join(conflicts)


def _normalise_context_value(value: object) -> object:
    missing = bool(pd.Series((value,), dtype="object").isna().iat[0])
    if missing:
        return None
    if isinstance(value, str):
        token = value.strip()
        return None if token == "" else token
    return value


__all__ = [
    "DISPLAY_SITE_CONTEXT_COLUMNS",
    "enforce_unique_display_site_identity_rows",
]
