"""Single-organism coherence helpers for analysis-ready datasets."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

import pandas as pd

from phospy.provenance.models import RunProvenance
from phospy.science.references.models import Organism
from phospy.science.sites.organisms import normalize_organism
from phospy.science.sites.site_keys import decode_site_key, encode_site_key

ErrorType = TypeVar("ErrorType", bound=Exception)


@dataclass(frozen=True, slots=True)
class NormalizedDatasetOrganismState:
    """Dataset tables after organism normalization."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame


def normalize_dataset_organism_state(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    error_type: type[ErrorType],
) -> NormalizedDatasetOrganismState:
    """Normalize organism-bearing site-key and row metadata values."""

    normalized_phospho = _normalize_site_key_indexed_frame(
        phospho,
        field_name="dataset.phospho.index",
        error_type=error_type,
    )
    normalized_site_metadata = _normalize_site_key_indexed_frame(
        site_metadata,
        field_name="dataset.site_metadata.index",
        error_type=error_type,
    )
    normalized_site_metadata = _normalize_site_key_column(
        normalized_site_metadata,
        field_name="dataset.site_metadata",
        error_type=error_type,
    )
    normalized_site_metadata = _normalize_row_organism_column(
        normalized_site_metadata,
        field_name="dataset.site_metadata",
        error_type=error_type,
    )
    return NormalizedDatasetOrganismState(
        phospho=normalized_phospho,
        site_metadata=normalized_site_metadata,
    )


def resolve_single_dataset_organism(
    *,
    site_metadata: pd.DataFrame,
    organism: Organism | None,
    error_type: type[ErrorType],
    preview_limit: int = 5,
) -> Organism:
    """Require all rows to describe one organism and resolve dataset organism."""

    row_organisms: list[tuple[object, Organism]] = []
    for row_id, raw_value in site_metadata.loc[:, "organism"].items():
        row_organisms.append(
            (
                row_id,
                normalize_organism(
                    raw_value,
                    field_name=f"dataset.site_metadata[{row_id!r}].organism",
                    error_type=error_type,
                ),
            )
        )
    unique_organisms = tuple(dict.fromkeys(item[1] for item in row_organisms))
    if len(unique_organisms) != 1:
        examples = _format_row_organism_examples(
            row_organisms,
            preview_limit=preview_limit,
        )
        raise error_type(
            "dataset.site_metadata organism values must resolve to one "
            "Organism; mixed-organism AnalysisReadyPhosphoDataset rows are not "
            f"supported; row_examples=[{examples}]"
        )
    resolved = unique_organisms[0]
    if organism is None:
        return resolved
    if organism is resolved:
        return organism

    mismatched_rows = [
        (row_id, row_organism)
        for row_id, row_organism in row_organisms
        if row_organism is not organism
    ]
    examples = _format_row_organism_examples(
        mismatched_rows,
        preview_limit=preview_limit,
    )
    raise error_type(
        "dataset.organism must match every dataset.site_metadata organism row; "
        f"dataset.organism={organism.value!r}; row_examples=[{examples}]"
    )


def require_dataset_provenance_organism_coherence(
    *,
    organism: Organism,
    provenance: RunProvenance | None,
    error_type: type[ErrorType],
) -> None:
    """Require supplied run/reference provenance to agree with dataset organism."""

    if provenance is None:
        return
    values: list[tuple[str, object]] = [("dataset.organism", organism)]
    if provenance.reference_context is not None:
        values.append(
            (
                "dataset.provenance.reference_context.organism",
                provenance.reference_context.organism,
            )
        )
    if provenance.reference is not None:
        values.append(
            (
                "dataset.provenance.reference.organism",
                provenance.reference.organism,
            )
        )
        if provenance.reference.reference_context is not None:
            values.append(
                (
                    "dataset.provenance.reference.reference_context.organism",
                    provenance.reference.reference_context.organism,
                )
            )
    _require_same_dataset_organism(
        values=values,
        error_type=error_type,
    )


def _normalize_site_key_indexed_frame(
    frame: pd.DataFrame,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> pd.DataFrame:
    if not _all_values_are_strict_encoded_site_keys(frame.index.tolist()):
        return frame
    normalized_index = _normalize_site_key_index(
        frame.index,
        field_name=field_name,
        error_type=error_type,
    )
    if normalized_index.equals(frame.index):
        return frame
    normalized = frame.copy(deep=False)
    normalized.index = normalized_index
    return normalized


def _normalize_site_key_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> pd.Index:
    values: list[str] = []
    for position, raw_value in enumerate(index.tolist()):
        key = decode_site_key(
            raw_value,
            field_name=f"{field_name}[{position}]",
            error_type=error_type,
        )
        values.append(encode_site_key(key))
    return pd.Index(values, name=index.name)


def _normalize_site_key_column(
    site_metadata: pd.DataFrame,
    *,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "site_key",
) -> pd.DataFrame:
    if column_name not in site_metadata.columns:
        return site_metadata
    if not _all_values_are_strict_encoded_site_keys(
        site_metadata.loc[:, column_name].tolist()
    ):
        return site_metadata
    values: list[str] = []
    for row_id, raw_value in site_metadata.loc[:, column_name].items():
        key = decode_site_key(
            raw_value,
            field_name=f"{field_name}.{column_name}[{row_id!r}]",
            error_type=error_type,
        )
        values.append(encode_site_key(key))
    if values == site_metadata.loc[:, column_name].astype(str).tolist():
        return site_metadata
    normalized = site_metadata.copy(deep=False)
    normalized.loc[:, column_name] = values
    return normalized


def _all_values_are_strict_encoded_site_keys(values: list[object]) -> bool:
    if not values:
        return False
    for value in values:
        if not isinstance(value, str):
            return False
        if value != value.strip():
            return False
        if not value.startswith("phospy:v1|"):
            return False
    return True


def _require_same_dataset_organism(
    *,
    values: list[tuple[str, object]],
    error_type: type[ErrorType],
) -> None:
    normalized = [
        (
            field_name,
            normalize_organism(
                value,
                field_name=field_name,
                error_type=error_type,
            ),
            value,
        )
        for field_name, value in values
    ]
    expected_field, expected_organism, _ = normalized[0]
    conflicts = [
        (field_name, organism, raw_value)
        for field_name, organism, raw_value in normalized[1:]
        if organism is not expected_organism
    ]
    if not conflicts:
        return
    conflict_text = "; ".join(
        f"{field_name}={_format_organism_value(raw_value)!r}"
        f" resolved_to={organism.value!r}"
        for field_name, organism, raw_value in conflicts
    )
    raise error_type(
        "dataset organism identity conflict; "
        f"{expected_field}={expected_organism.value!r}; {conflict_text}"
    )


def _format_organism_value(value: object) -> str:
    if isinstance(value, Organism):
        return value.value
    return str(value)


def _normalize_row_organism_column(
    site_metadata: pd.DataFrame,
    *,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "organism",
) -> pd.DataFrame:
    if column_name not in site_metadata.columns:
        return site_metadata
    values: list[str] = []
    for row_id, raw_value in site_metadata.loc[:, column_name].items():
        values.append(
            normalize_organism(
                raw_value,
                field_name=f"{field_name}[{row_id!r}].{column_name}",
                error_type=error_type,
            ).value
        )
    if values == site_metadata.loc[:, column_name].astype(str).tolist():
        return site_metadata
    normalized = site_metadata.copy(deep=False)
    normalized.loc[:, column_name] = values
    return normalized


def _format_row_organism_examples(
    rows: list[tuple[object, Organism]],
    *,
    preview_limit: int,
) -> str:
    preview = ", ".join(
        f"{row_id!r}:organism={organism.value!r}"
        for row_id, organism in rows[:preview_limit]
    )
    suffix = "" if len(rows) <= preview_limit else " ..."
    return f"{preview}{suffix}"


__all__ = [
    "NormalizedDatasetOrganismState",
    "normalize_dataset_organism_state",
    "require_dataset_provenance_organism_coherence",
    "resolve_single_dataset_organism",
]
