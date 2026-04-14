from __future__ import annotations

import json
from collections.abc import Mapping, Sequence

import pandas as pd

from .constants import TOP_KINASE_CANDIDATES_COLUMN, TOP_KINASE_WEIGHTS_COLUMN

__all__ = [
    "normalize_top_kinase_weights",
    "serialize_site_assignments_for_export",
    "serialize_top_kinase_candidates",
    "serialize_top_kinase_weights",
]


def serialize_top_kinase_candidates(value: object) -> str:
    """Serialize candidate-kinase tuples/lists into a deterministic JSON array."""

    if isinstance(value, str):
        return value
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return json.dumps([str(kinase) for kinase in value])

    msg = f"{TOP_KINASE_CANDIDATES_COLUMN} must be a sequence of kinase names"
    raise TypeError(msg)


def serialize_top_kinase_weights(value: object) -> str:
    """Serialize weighted tie assignments into a deterministic JSON object."""

    if isinstance(value, str):
        return value

    if isinstance(value, Mapping):
        items = [(str(kinase), float(weight)) for kinase, weight in value.items()]
        return json.dumps(dict(items))

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        items: list[tuple[str, float]] = []
        for pair in value:
            if not (isinstance(pair, Sequence) and not isinstance(pair, (str, bytes))):
                msg = (
                    f"{TOP_KINASE_WEIGHTS_COLUMN} sequence entries must be "
                    "(kinase, weight)"
                )
                raise TypeError(msg)
            pair_values = tuple(pair)
            if len(pair_values) != 2:
                msg = (
                    f"{TOP_KINASE_WEIGHTS_COLUMN} sequence entries must be "
                    "(kinase, weight)"
                )
                raise TypeError(msg)
            kinase, weight = pair_values
            items.append((str(kinase), float(weight)))
        return json.dumps(dict(items))

    msg = (
        f"{TOP_KINASE_WEIGHTS_COLUMN} must be a mapping or sequence of (kinase, weight)"
    )
    raise TypeError(msg)


def normalize_top_kinase_weights(value: object) -> tuple[tuple[str, float], ...]:
    """Normalize top-kinase weights into an ordered tuple-of-pairs payload."""

    if isinstance(value, str):
        decoded = json.loads(value)
        return normalize_top_kinase_weights(decoded)

    if isinstance(value, Mapping):
        return tuple((str(kinase), float(weight)) for kinase, weight in value.items())

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        items: list[tuple[str, float]] = []
        for pair in value:
            if not (isinstance(pair, Sequence) and not isinstance(pair, (str, bytes))):
                msg = (
                    f"{TOP_KINASE_WEIGHTS_COLUMN} sequence entries must be "
                    "(kinase, weight)"
                )
                raise TypeError(msg)
            pair_values = tuple(pair)
            if len(pair_values) != 2:
                msg = (
                    f"{TOP_KINASE_WEIGHTS_COLUMN} sequence entries must be "
                    "(kinase, weight)"
                )
                raise TypeError(msg)
            kinase, weight = pair_values
            items.append((str(kinase), float(weight)))
        return tuple(items)

    msg = f"{TOP_KINASE_WEIGHTS_COLUMN} must be a mapping, JSON object, or sequence"
    raise TypeError(msg)


def serialize_site_assignments_for_export(
    site_assignments: pd.DataFrame,
) -> pd.DataFrame:
    """Return a CSV-ready copy with tie fields encoded as JSON strings."""

    exported = site_assignments.copy(deep=True)
    if TOP_KINASE_CANDIDATES_COLUMN in exported.columns:
        exported[TOP_KINASE_CANDIDATES_COLUMN] = exported[
            TOP_KINASE_CANDIDATES_COLUMN
        ].map(serialize_top_kinase_candidates)
    if TOP_KINASE_WEIGHTS_COLUMN in exported.columns:
        exported[TOP_KINASE_WEIGHTS_COLUMN] = exported[TOP_KINASE_WEIGHTS_COLUMN].map(
            serialize_top_kinase_weights
        )
    return exported
