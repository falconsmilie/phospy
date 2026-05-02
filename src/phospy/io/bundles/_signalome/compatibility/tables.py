"""Signalome bundle table normalization helpers."""

from __future__ import annotations

import ast


def normalize_module_assignments_table(table):
    """Normalize tuple/list/dict-serialized signalome assignment fields."""

    normalized = table.copy(deep=True)
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
