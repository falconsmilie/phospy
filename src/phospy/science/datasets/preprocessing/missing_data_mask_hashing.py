"""Stable fingerprints for missing-data mechanism masks."""

from __future__ import annotations

import pandas as pd

from phospy.provenance.hashing import hash_table_tolerance


def hash_imputation_mask(mask: pd.DataFrame) -> str:
    """Return a stable fingerprint for a binary imputation mask."""

    return hash_table_tolerance(
        mask.astype("int8"),
        name="missing_data.imputation_mask",
    )


def hash_group_aware_imputation_mask(mask: pd.DataFrame) -> str:
    """Return a bundle-stable fingerprint for a group-aware union mask."""

    return hash_table_tolerance(
        _canonical_mask(mask),
        name="missing_data.imputation_mask",
    )


def hash_group_aware_mechanism_mask(
    mask: pd.DataFrame,
    *,
    mechanism: str,
    mask_kind: str,
) -> str:
    """Return a stable, mechanism-labelled group-aware mask fingerprint."""

    if mechanism not in {"knn", "minprob"}:
        raise ValueError("group-aware mechanism must be 'knn' or 'minprob'")
    if mask_kind not in {"target", "imputation"}:
        raise ValueError("group-aware mask kind must be 'target' or 'imputation'")
    return hash_table_tolerance(
        _canonical_mask(mask),
        name=f"missing_data.group_aware.{mechanism}_{mask_kind}_mask",
    )


def _canonical_mask(mask: pd.DataFrame) -> pd.DataFrame:
    canonical = mask.astype("int8")
    canonical.index = pd.Index(canonical.index.tolist(), dtype="object", name=None)
    canonical.columns = pd.Index(canonical.columns.tolist(), dtype="object", name=None)
    return canonical


__all__ = [
    "hash_group_aware_imputation_mask",
    "hash_group_aware_mechanism_mask",
    "hash_imputation_mask",
]
