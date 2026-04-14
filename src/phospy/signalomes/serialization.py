from __future__ import annotations

import json
from collections.abc import Mapping, Sequence

import pandas as pd

__all__ = [
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

    msg = "top_kinase_candidates must be a sequence of kinase names"
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
                msg = "top_kinase_weights sequence entries must be (kinase, weight)"
                raise TypeError(msg)
            pair_values = tuple(pair)
            if len(pair_values) != 2:
                msg = "top_kinase_weights sequence entries must be (kinase, weight)"
                raise TypeError(msg)
            kinase, weight = pair_values
            items.append((str(kinase), float(weight)))
        return json.dumps(dict(items))

    msg = "top_kinase_weights must be a mapping or sequence of (kinase, weight)"
    raise TypeError(msg)


def serialize_site_assignments_for_export(
    site_assignments: pd.DataFrame,
) -> pd.DataFrame:
    """Return a CSV-ready copy with tie fields encoded as JSON strings."""

    exported = site_assignments.copy(deep=True)
    if "top_kinase_candidates" in exported.columns:
        exported["top_kinase_candidates"] = exported["top_kinase_candidates"].map(
            serialize_top_kinase_candidates
        )
    if "top_kinase_weights" in exported.columns:
        exported["top_kinase_weights"] = exported["top_kinase_weights"].map(
            serialize_top_kinase_weights
        )
    return exported
