from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd
import pandas.testing as pdt

from phospy.workflows.signalome.science import (
    build_expanded_signalome_table,
    build_signalome_module_table,
)


def _legacy_build_signalome_module_table(
    *,
    module_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    kinase_order: Sequence[str],
) -> pd.DataFrame:
    module_index = pd.Index(
        sorted(
            {
                int(value)
                for value in module_assignments.loc[:, "module_id"]
                if int(value) > 0
            }
        ),
        name="module_id",
    )
    kinase_index = pd.Index([str(kinase) for kinase in kinase_order], name="kinase")
    module_table = pd.DataFrame(
        0.0,
        index=module_index.copy(),
        columns=kinase_index.copy(),
    )

    protein_to_module = (
        module_assignments.loc[:, ["protein_id", "module_id"]]
        .drop_duplicates(subset=["protein_id"])
        .set_index("protein_id")
        .loc[:, "module_id"]
        .astype("int64")
    )
    protein_to_module = protein_to_module.loc[protein_to_module > 0]
    site_to_protein = module_assignments.loc[:, "protein_id"].astype(str)
    site_to_protein.index = pd.Index(
        site_to_protein.index.astype(str),
        name="site_id",
    )

    for kinase in kinase_index:
        substrate_sites = pd.Index(
            [str(site_id) for site_id in kinase_substrates.get(str(kinase), ())],
            name="site_id",
        )
        if substrate_sites.empty:
            continue
        substrate_proteins = (
            site_to_protein.reindex(substrate_sites).dropna().astype(str)
        )
        if substrate_proteins.empty:
            continue
        unique_proteins = pd.Index(sorted(set(substrate_proteins.tolist())))
        module_hits = (
            protein_to_module.reindex(unique_proteins).dropna().astype("int64")
        )
        if module_hits.empty:
            continue
        counts = module_hits.value_counts().astype(float)
        module_table.loc[counts.index.astype(int), kinase] = counts.to_numpy(
            dtype=float, copy=False
        )

    row_totals = module_table.sum(axis=1)
    non_zero_rows = row_totals > 0.0
    if non_zero_rows.any():
        module_table.loc[non_zero_rows] = (
            module_table.loc[non_zero_rows].div(row_totals.loc[non_zero_rows], axis=0)
            * 100.0
        )
    return module_table.astype(float).round(3)


def test_build_signalome_module_table_matches_legacy_semantics() -> None:
    module_assignments = pd.DataFrame(
        {
            "protein_id": ["P1", "P1", "P2", "P3", "P4", "P5", "P6", "P3"],
            "module_id": [1, 1, 2, 2, 0, 3, 3, 3],
        },
        index=pd.Index(
            ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"],
            name="site_id",
        ),
    )
    kinase_substrates = {
        "K1": ("S1", "S2", "S4", "S8"),
        "K2": ("S3", "S5", "S6", "MISSING"),
        "K3": ("S7",),
    }
    kinase_order = ["K2", "K1", "K3", "K1"]

    expected = _legacy_build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates=kinase_substrates,
        kinase_order=kinase_order,
    )
    observed = build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates=kinase_substrates,
        kinase_order=kinase_order,
    )

    pdt.assert_frame_equal(observed, expected, check_dtype=False)


def test_build_signalome_module_table_weighted_top_propagates_fractional_support() -> (
    None
):
    module_assignments = pd.DataFrame(
        {
            "protein_id": ["P1", "P1", "P2", "P2"],
            "module_id": [1, 1, 2, 2],
            "top_kinase_weights": [
                (("K1", 0.5), ("K2", 0.5)),
                (("K1", 1.0),),
                (("K2", 1.0),),
                (("K2", 1.0),),
            ],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )

    weighted = build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates={"K1": (), "K2": ()},
        kinase_order=["K1", "K2"],
        assignment_policy="weighted_top",
    )
    binary = build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates={"K1": (), "K2": ()},
        kinase_order=["K1", "K2"],
        assignment_policy="cutoff_binary",
    )

    assert weighted.loc[1, "K1"] == 66.667
    assert weighted.loc[1, "K2"] == 33.333
    assert weighted.loc[2, "K2"] == 100.0
    assert binary.loc[1, "K1"] == 0.0
    assert binary.loc[1, "K2"] == 0.0


def test_build_expanded_signalome_table_tracks_membership_and_site_order() -> None:
    module_assignments = pd.DataFrame(
        {
            "protein_id": ["P1", "P2", "P3", "P4"],
            "module_id": [2, 1, 2, 3],
            "top_kinase": ["K2", "K1", "K2", "K3"],
            "top_score": [0.91, 0.93, 0.92, 0.88],
            "top_kinase_weights": [
                (("K2", 1.0),),
                (("K1", 1.0),),
                (("K2", 1.0),),
                (("K3", 1.0),),
            ],
        },
        index=pd.Index(["S3", "S1", "S4", "S2"], name="site_id"),
    )
    signalome_modules = pd.DataFrame(
        {
            "K1": [5.0, 85.0, 0.0],
            "K2": [90.0, 10.0, 0.0],
            "K3": [0.0, 0.0, 100.0],
        },
        index=pd.Index([1, 2, 3], name="module_id"),
    )
    kinase_network_edges = pd.DataFrame(
        {
            "source_kinase": ["K1"],
            "target_kinase": ["K2"],
            "correlation": [0.95],
        }
    )
    kinase_substrates = {
        "K1": ("S1", "S2"),
        "K2": ("S3", "S4"),
    }

    expanded = build_expanded_signalome_table(
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates=kinase_substrates,
        assignment_policy="cutoff_binary",
    )

    k1_sites = expanded.loc[
        (expanded.loc[:, "kinase"] == "K1")
        & (expanded.loc[:, "row_kind"] == "site")
        & (expanded.loc[:, "site_id"] != ""),
        :,
    ]
    assert k1_sites.loc[:, "site_id"].tolist() == ["S3", "S1", "S4"]
    assert k1_sites.loc[:, "site_order"].tolist() == [0, 1, 2]
    assert k1_sites.loc[:, "linked_kinases"].iloc[0] == '["K1","K2"]'
    assert k1_sites.loc[:, "regulated_module_ids"].iloc[0] == "[1,2]"


def test_build_expanded_signalome_table_weighted_top_uses_fractional_site_support() -> (
    None
):
    module_assignments = pd.DataFrame(
        {
            "protein_id": ["P1", "P1", "P2"],
            "module_id": [1, 1, 2],
            "top_kinase": ["K1", "K1", "K2"],
            "top_score": [0.95, 0.96, 0.92],
            "top_kinase_weights": [
                (("K1", 0.5), ("K2", 0.5)),
                (("K1", 1.0),),
                (("K2", 1.0),),
            ],
        },
        index=pd.Index(["S1", "S2", "S3"], name="site_id"),
    )
    signalome_modules = pd.DataFrame(
        {"K1": [90.0, 0.0], "K2": [10.0, 100.0]},
        index=pd.Index([1, 2], name="module_id"),
    )
    kinase_network_edges = pd.DataFrame(
        {
            "source_kinase": ["K1"],
            "target_kinase": ["K2"],
            "correlation": [0.91],
        }
    )

    expanded = build_expanded_signalome_table(
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates={"K1": (), "K2": ()},
        assignment_policy="weighted_top",
    )

    site_s1 = expanded.loc[
        (expanded.loc[:, "kinase"] == "K1") & (expanded.loc[:, "site_id"] == "S1"),
        :,
    ].iloc[0]
    assert site_s1["support_kinases"] == '["K1","K2"]'
    assert float(site_s1["support_weight"]) == 1.0
