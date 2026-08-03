"""Protein-module derivation from site-level module assignments."""

from __future__ import annotations

import pandas as pd

from phospy.science.signalomes.constants import (
    MODULE_ID_COLUMN,
    PROTEIN_COLUMN,
    SITE_ID_COLUMN,
)


def derive_protein_modules(
    *,
    site_clusters: pd.Series,
    site_to_protein_group_id: pd.Series | None = None,
    site_to_protein: pd.Series | None = None,
) -> pd.Series:
    """Collapse site-level clusters into protein-group-level module assignments."""

    if site_to_protein_group_id is None:
        if site_to_protein is None:
            raise TypeError("derive_protein_modules requires site_to_protein_group_id")
        site_to_protein_group_id = site_to_protein
    elif site_to_protein is not None and not site_to_protein.equals(
        site_to_protein_group_id
    ):
        raise ValueError(
            "derive_protein_modules received conflicting site_to_protein_group_id "
            "and legacy site_to_protein mappings"
        )
    aligned_site_to_protein = site_to_protein_group_id.copy()
    aligned_site_to_protein.index = pd.Index(
        aligned_site_to_protein.index.astype(str),
        name=SITE_ID_COLUMN,
    )
    cluster_index = pd.Index(site_clusters.index.astype(str), name=SITE_ID_COLUMN)
    missing_sites: list[str] = [
        str(site_id)
        for site_id in cluster_index
        if site_id not in aligned_site_to_protein.index
    ]
    if missing_sites:
        preview = ", ".join(missing_sites[:3])
        suffix = "..." if len(missing_sites) > 3 else ""
        raise ValueError(
            "site_to_protein_group_id is missing clustered site mappings: "
            f"{preview}{suffix}"
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
