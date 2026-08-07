from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    EnrichmentWorkflowResult,
    PtmSetCollection,
)

pytestmark = pytest.mark.contract


def test_public_consumer_runs_site_key_enrichment_smoke() -> None:
    site_keys = (
        "rat|MAPK14|Y182",
        "rat|GSK3B|S9",
        "rat|AKT1|T308",
    )
    collection = PtmSetCollection(
        sets={
            "MAPK_GSK3_AXIS": site_keys[:2],
            "AKT_ONLY": (site_keys[2],),
        },
        identifier_kind="site_key",
        term_names={
            "MAPK_GSK3_AXIS": "MAPK/GSK3 axis",
            "AKT_ONLY": "AKT-only control",
        },
        source_name="contract in-memory PTM sets",
        source_version="2026-08",
    )

    result = EnrichmentWorkflow().run(
        EnrichmentWorkflowRequest(
            identifier_column="site_key",
            identifier_kind="site_key",
            set_collection=collection,
            input_table=pd.DataFrame({"site_key": [site_keys[0], f" {site_keys[1]} "]}),
            background_universe=site_keys,
            config=EnrichmentConfig(min_set_size=1, max_set_size=3),
        )
    )

    assert isinstance(result, EnrichmentWorkflowResult)
    assert result.identifier_kind == "site_key"
    assert result.set_collection_summary["collection_kind"] == "ptm_set"
    assert result.background_summary["selected_identifier_source"] == "input_table"

    table = result.table
    assert {
        "term_id",
        "identifier_kind",
        "input_overlap_count",
        "overlap_identifiers",
        "p_value",
        "adjusted_p_value",
    } <= set(table.columns)

    rows = table.set_index("term_id")
    matched = rows.loc["MAPK_GSK3_AXIS"]
    assert matched["identifier_kind"] == "site_key"
    assert int(matched["input_overlap_count"]) == 2
    assert set(matched["overlap_identifiers"]) == set(site_keys[:2])
