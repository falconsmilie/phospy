"""Protein-module derivation from site-level module assignments."""

from __future__ import annotations

import pandas as pd

from phospy.signalomes.constants import (
    MODULE_ID_COLUMN,
    PROTEIN_COLUMN,
    SITE_ID_COLUMN,
)


def derive_protein_modules(
    *,
    site_clusters: pd.Series,
    site_to_protein: pd.Series,
) -> pd.Series:
    """Collapse site-level clusters into protein-level module assignments."""

    aligned_site_to_protein = site_to_protein.copy()
    aligned_site_to_protein.index = pd.Index(
        aligned_site_to_protein.index.astype(str),
        name=SITE_ID_COLUMN,
    )
    cluster_index = pd.Index(site_clusters.index.astype(str), name=SITE_ID_COLUMN)
    missing_sites = [
        site_id
        for site_id in cluster_index
        if site_id not in aligned_site_to_protein.index
    ]
    if missing_sites:
        preview = ", ".join(missing_sites[:3])
        suffix = "..." if len(missing_sites) > 3 else ""
        raise ValueError(
            f"site_to_protein is missing clustered site mappings: {preview}{suffix}"
        )
    aligned_site_to_protein = aligned_site_to_protein.loc[cluster_index].astype(str)

    proteins = pd.Index(aligned_site_to_protein.tolist(), dtype=object)
    membership = pd.crosstab(site_clusters, proteins)
    membership = (membership > 0).astype(int)

    assignments: dict[str, int] = {}
    pattern_to_module: dict[tuple[int, ...], int] = {}
    next_module_id = 1
    for protein in membership.columns:
        pattern = tuple(int(value) for value in membership.loc[:, protein].tolist())
        if pattern not in pattern_to_module:
            pattern_to_module[pattern] = next_module_id
            next_module_id += 1
        assignments[str(protein)] = pattern_to_module[pattern]

    protein_modules = pd.Series(assignments, dtype="int64", name=MODULE_ID_COLUMN)
    protein_modules.index.name = PROTEIN_COLUMN
    return protein_modules


__all__ = ["derive_protein_modules"]
