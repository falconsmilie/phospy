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


__all__ = ["hash_imputation_mask"]
