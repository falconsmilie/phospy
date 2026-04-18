from __future__ import annotations

import pandas as pd

from phospy.workflows.signalome.science import (
    LEXICOGRAPHIC_TIE_BREAK_POLICY,
    build_module_assignments,
)


def test_build_module_assignments_surfaces_site_level_score_ties() -> None:
    prediction_matrix = pd.DataFrame(
        {
            "K2": [0.8, 0.1],
            "K1": [0.8, 0.9],
        },
        index=pd.Index(["P1;S1;", "P2;S2;"], name="site_id"),
        dtype=float,
    )
    site_to_protein = pd.Series(
        ["P1", "P2"],
        index=prediction_matrix.index.copy(),
        name="protein_id",
        dtype=str,
    )

    assignments = build_module_assignments(
        prediction_matrix=prediction_matrix,
        site_to_protein=site_to_protein,
    )

    tied_site = assignments.loc["P1;S1;"]
    assert tied_site["top_kinase"] == "K1"
    assert tied_site["top_kinase_candidates"] == ("K1", "K2")
    assert int(tied_site["top_kinase_tie_count"]) == 2
    assert bool(tied_site["top_kinase_is_ambiguous"])
    assert tied_site["top_kinase_selection_policy"] == LEXICOGRAPHIC_TIE_BREAK_POLICY


def test_build_module_assignments_surfaces_protein_level_equal_ranking() -> None:
    prediction_matrix = pd.DataFrame(
        {
            "K1": [0.9, 0.1],
            "K2": [0.1, 0.9],
        },
        index=pd.Index(["P1;S1;", "P1;S2;"], name="site_id"),
        dtype=float,
    )
    site_to_protein = pd.Series(
        ["P1", "P1"],
        index=prediction_matrix.index.copy(),
        name="protein_id",
        dtype=str,
    )

    assignments = build_module_assignments(
        prediction_matrix=prediction_matrix,
        site_to_protein=site_to_protein,
    )

    assert (
        assignments.loc["P1;S1;", "module_id"] == assignments.loc["P1;S2;", "module_id"]
    )
    assert assignments.loc["P1;S1;", "module_top_kinase"] == "K1"
    assert assignments.loc["P1;S2;", "module_top_kinase"] == "K1"
    assert assignments.loc["P1;S1;", "module_top_kinase_candidates"] == ("K1", "K2")
    assert assignments.loc["P1;S2;", "module_top_kinase_candidates"] == ("K1", "K2")
    assert int(assignments.loc["P1;S1;", "module_top_kinase_tie_count"]) == 2
    assert int(assignments.loc["P1;S2;", "module_top_kinase_tie_count"]) == 2
    assert bool(assignments.loc["P1;S1;", "module_top_kinase_is_ambiguous"])
    assert bool(assignments.loc["P1;S2;", "module_top_kinase_is_ambiguous"])
    assert (
        assignments.loc["P1;S1;", "module_top_kinase_selection_policy"]
        == LEXICOGRAPHIC_TIE_BREAK_POLICY
    )
    assert (
        assignments.loc["P1;S2;", "module_top_kinase_selection_policy"]
        == LEXICOGRAPHIC_TIE_BREAK_POLICY
    )
