from __future__ import annotations

import pandas as pd

from phospy.io.bundles._signalome.reconstruction import (
    _normalize_site_metadata_for_dataset_contract,
)


def test_signalome_bundle_site_metadata_restores_site_key_column_from_site_key_index() -> (
    None
):
    table = pd.DataFrame(
        {"display_id": ["MAPK14;Y182;"]},
        index=pd.Index(["test.site_key.row_a"], name="site_key"),
    )

    normalized = _normalize_site_metadata_for_dataset_contract(table)

    assert normalized.loc[:, "site_key"].tolist() == ["test.site_key.row_a"]


def test_signalome_bundle_site_metadata_does_not_repair_display_index_identity() -> (
    None
):
    table = pd.DataFrame(
        {"display_id": ["MAPK14;Y182;"]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    normalized = _normalize_site_metadata_for_dataset_contract(table)

    assert "site_key" not in normalized.columns
