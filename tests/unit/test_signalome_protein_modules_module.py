from __future__ import annotations

import pandas as pd
import pytest

from phospy.science.signalomes.clustering.protein_modules import derive_protein_modules


def test_derive_protein_modules_groups_by_site_membership_patterns() -> None:
    site_clusters = pd.Series(
        [1, 2, 1, 2, 3],
        index=pd.Index(["S1", "S2", "S3", "S4", "S5"], name="site_id"),
        dtype="int64",
    )
    site_to_protein = pd.Series(
        ["P1", "P1", "P2", "P2", "P3"],
        index=site_clusters.index.copy(),
        dtype=str,
    )

    modules = derive_protein_modules(
        site_clusters=site_clusters,
        site_to_protein=site_to_protein,
    )

    assert modules.at["P1"] == modules.at["P2"]
    assert modules.at["P3"] != modules.at["P1"]


def test_derive_protein_modules_fails_when_mapping_is_missing_sites() -> None:
    site_clusters = pd.Series(
        [1, 2],
        index=pd.Index(["S1", "S2"], name="site_id"),
        dtype="int64",
    )
    site_to_protein = pd.Series(
        ["P1"],
        index=pd.Index(["S1"], name="site_id"),
        dtype=str,
    )

    with pytest.raises(ValueError, match="missing clustered site mappings"):
        derive_protein_modules(
            site_clusters=site_clusters,
            site_to_protein=site_to_protein,
        )
