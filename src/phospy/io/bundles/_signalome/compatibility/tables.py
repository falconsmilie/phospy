"""Signalome bundle table normalization helpers."""

from __future__ import annotations

import ast


def normalize_module_assignments_table(table):
    """Normalize tuple/list/dict-serialized signalome assignment fields."""

    normalized = table.copy(deep=True)
    normalized = _normalize_identity_columns(normalized)
    normalized = _normalize_string_columns(normalized)
    candidate_columns = [
        str(column)
        for column in normalized.columns
        if str(column).endswith("_candidates")
    ]
    for candidates_column in candidate_columns:
        candidates_index = normalized.columns.get_loc(candidates_column)
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
        weight_index = normalized.columns.get_loc(weight_column)
        weights = (
            normalized.loc[:, weight_column].map(_parse_kinase_weights).astype(object)
        )
        normalized = normalized.drop(columns=[weight_column])
        normalized.insert(weight_index, weight_column, weights)
    return normalized


def _normalize_identity_columns(table):
    normalized = table.copy(deep=True)
    for column_name in ("site_key", "display_id", "site_id"):
        duplicate_name = f"{column_name}.1"
        if (
            column_name not in normalized.columns
            and duplicate_name in normalized.columns
        ):
            normalized = normalized.rename(columns={duplicate_name: column_name})
    index_name = str(normalized.index.name) if normalized.index.name is not None else ""
    if index_name:
        duplicate_index_name = f"{index_name}.1"
        if (
            index_name not in normalized.columns
            and duplicate_index_name in normalized.columns
        ):
            normalized = normalized.rename(columns={duplicate_index_name: index_name})
    return normalized


def _normalize_string_columns(table):
    normalized = table.copy(deep=True)
    for column_name in (
        "site_key",
        "display_id",
        "site_id",
        "gene_symbol",
        "site",
        "protein_id",
        "protein_accession",
        "isoform_id",
        "top_kinase",
        "top_kinase_selection_policy",
        "module_top_kinase",
        "module_top_kinase_selection_policy",
    ):
        if column_name in normalized.columns:
            column_index = normalized.columns.get_loc(column_name)
            series = (
                normalized.loc[:, column_name].astype(object).fillna("").astype(str)
            )
            normalized = normalized.drop(columns=[column_name])
            normalized.insert(column_index, column_name, series)
    return normalized


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
