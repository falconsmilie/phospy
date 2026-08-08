"""Signalome bundle table normalization helpers."""

from __future__ import annotations

import ast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.signalomes.constants import (
    LEGACY_PROTEIN_GROUP_ID_COLUMN,
    PROTEIN_GROUP_ID_COLUMN,
)


def normalize_module_assignments_table(table: pd.DataFrame) -> pd.DataFrame:
    """Normalize tuple/list/dict-serialized signalome assignment fields."""

    normalized = table.copy(deep=True)
    normalized = _normalize_current_site_key_column(normalized)
    normalized = migrate_signalome_protein_group_id_column(
        normalized,
        field_name="bundle signalome module_assignments",
    )
    normalized = _normalize_string_columns(normalized)
    candidate_columns = [
        str(column)
        for column in normalized.columns
        if str(column).endswith("_candidates")
    ]
    for candidates_column in candidate_columns:
        candidates_index = _unique_column_position(
            normalized,
            column_name=candidates_column,
            field_name="bundle signalome module_assignments",
        )
        candidates = (
            normalized.loc[:, candidates_column]
            .map(_parse_kinase_candidates)
            .astype(object)
        )
        normalized = normalized.drop(columns=[candidates_column])
        normalized.insert(candidates_index, candidates_column, candidates)
    weight_columns = [
        str(column) for column in normalized.columns if str(column).endswith("_weights")
    ]
    for weight_column in weight_columns:
        weight_index = _unique_column_position(
            normalized,
            column_name=weight_column,
            field_name="bundle signalome module_assignments",
        )
        weights = (
            normalized.loc[:, weight_column].map(_parse_kinase_weights).astype(object)
        )
        normalized = normalized.drop(columns=[weight_column])
        normalized.insert(weight_index, weight_column, weights)
    return normalized


def migrate_signalome_protein_group_id_column(
    table: pd.DataFrame,
    *,
    field_name: str,
) -> pd.DataFrame:
    """Migrate legacy Signalome grouping column name to protein_group_id."""

    has_current = PROTEIN_GROUP_ID_COLUMN in table.columns
    has_legacy = LEGACY_PROTEIN_GROUP_ID_COLUMN in table.columns
    if has_current and has_legacy:
        current = (
            table.loc[:, PROTEIN_GROUP_ID_COLUMN].fillna("").astype(str).str.strip()
        )
        legacy = (
            table.loc[:, LEGACY_PROTEIN_GROUP_ID_COLUMN]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        mismatch_mask = current.ne(legacy)
        if bool(mismatch_mask.any()):
            mismatch_rows = [
                str(row_id)
                for row_id in table.index[mismatch_mask.to_numpy()].astype(str).tolist()
            ]
            preview = ", ".join(mismatch_rows[:5])
            suffix = "" if len(mismatch_rows) <= 5 else " ..."
            raise PhosPyInputError(
                f"{field_name} has conflicting Signalome grouping columns "
                "protein_group_id and legacy protein_id; "
                f"mismatch_rows=[{preview}{suffix}]"
            )
        return table.drop(columns=[LEGACY_PROTEIN_GROUP_ID_COLUMN])
    if has_current or not has_legacy:
        return table
    return table.rename(
        columns={LEGACY_PROTEIN_GROUP_ID_COLUMN: PROTEIN_GROUP_ID_COLUMN}
    )


def _normalize_current_site_key_column(table: pd.DataFrame) -> pd.DataFrame:
    if "site_key" in table.columns or "site_key.1" not in table.columns:
        return table
    normalized = table.copy(deep=True)
    return normalized.rename(columns={"site_key.1": "site_key"})


def _normalize_string_columns(table: pd.DataFrame) -> pd.DataFrame:
    normalized = table.copy(deep=True)
    for column_name in (
        "site_key",
        "display_id",
        "gene_symbol",
        "site",
        "protein_group_id",
        "protein_accession",
        "isoform_id",
        "top_kinase",
        "top_kinase_selection_policy",
        "module_top_kinase",
        "module_top_kinase_selection_policy",
    ):
        if column_name in normalized.columns:
            column_index = _unique_column_position(
                normalized,
                column_name=column_name,
                field_name="bundle signalome module_assignments",
            )
            series = (
                normalized.loc[:, column_name].astype(object).fillna("").astype(str)
            )
            normalized = normalized.drop(columns=[column_name])
            normalized.insert(column_index, column_name, series)
    return normalized


def _unique_column_position(
    table: pd.DataFrame,
    *,
    column_name: str,
    field_name: str,
) -> int:
    position = table.columns.get_loc(column_name)
    if isinstance(position, int):
        return position
    raise PhosPyInputError(f"{field_name}.{column_name} must be a unique column")


def _parse_kinase_candidates(value: object) -> tuple[str, ...]:
    if isinstance(value, tuple):
        return tuple(str(item) for item in value)
    if isinstance(value, list):
        return tuple(str(item) for item in value)
    raw = str(value).strip()
    if raw == "" or raw.lower() == "nan":
        return ()
    try:
        parsed = ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return (raw,)
    if isinstance(parsed, tuple):
        return tuple(str(item) for item in parsed)
    if isinstance(parsed, list):
        return tuple(str(item) for item in parsed)
    return (str(parsed),)


def _parse_kinase_weights(value: object) -> tuple[tuple[str, float], ...]:
    if isinstance(value, dict):
        return tuple((str(key), float(weight)) for key, weight in value.items())
    if isinstance(value, (tuple, list)):
        return _normalize_kinase_weight_pairs(value)
    raw = str(value).strip()
    if raw == "" or raw.lower() == "nan":
        return ()
    try:
        parsed = ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return ()
    if isinstance(parsed, dict):
        return tuple((str(key), float(weight)) for key, weight in parsed.items())
    if isinstance(parsed, (tuple, list)):
        return _normalize_kinase_weight_pairs(parsed)
    return ()


def _normalize_kinase_weight_pairs(
    values: tuple[object, ...] | list[object],
) -> tuple[tuple[str, float], ...]:
    normalized_pairs: list[tuple[str, float]] = []
    for value in values:
        if not isinstance(value, (tuple, list)) or len(value) != 2:
            continue
        kinase, weight = value
        try:
            normalized_pairs.append((str(kinase), float(weight)))
        except (TypeError, ValueError):
            continue
    return tuple(normalized_pairs)
