from __future__ import annotations

import pandas as pd

from phospy.science.signalomes.science import (
    LEXICOGRAPHIC_TIE_BREAK_POLICY,
    NO_SUPPORT_SELECTION_POLICY,
    UNSUPPORTED_KINASE,
    build_module_assignments,
)
from tests.support.site_keys import site_key_index_from_display_ids


def _site_metadata(display_ids: list[str]) -> tuple[pd.Index, pd.DataFrame]:
    site_index = site_key_index_from_display_ids(display_ids)
    genes = [display_id.split(";")[0] for display_id in display_ids]
    sites = [display_id.split(";")[1] for display_id in display_ids]
    return site_index, pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": genes,
            "site": sites,
        },
        index=site_index.copy(),
    )


def test_build_module_assignments_surfaces_site_level_score_ties() -> None:
    site_index, site_metadata = _site_metadata(["P1;S1;", "P2;S2;"])
    prediction_matrix = pd.DataFrame(
        {
            "K2": [0.8, 0.1],
            "K1": [0.8, 0.9],
        },
        index=site_index.copy(),
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
        site_to_protein_group_id=site_to_protein,
        site_metadata=site_metadata,
    )

    tied_site = assignments.loc[site_index[0]]
    assert tied_site["top_kinase"] == "K1"
    assert tied_site["top_kinase_candidates"] == ("K1", "K2")
    assert tied_site["top_kinase_weights"] == (("K1", 0.5), ("K2", 0.5))
    assert int(tied_site["top_kinase_tie_count"]) == 2
    assert bool(tied_site["top_kinase_is_ambiguous"])
    assert tied_site["top_kinase_selection_policy"] == LEXICOGRAPHIC_TIE_BREAK_POLICY


def test_build_module_assignments_surfaces_protein_level_equal_ranking() -> None:
    site_index, site_metadata = _site_metadata(["P1;S1;", "P1;S2;"])
    prediction_matrix = pd.DataFrame(
        {
            "K1": [0.9, 0.1],
            "K2": [0.1, 0.9],
        },
        index=site_index.copy(),
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
        site_to_protein_group_id=site_to_protein,
        site_metadata=site_metadata,
    )
    p1s1, p1s2 = site_index.astype(str).tolist()

    assert assignments.loc[p1s1, "module_id"] == assignments.loc[p1s2, "module_id"]
    assert assignments.loc[p1s1, "module_top_kinase"] == "K1"
    assert assignments.loc[p1s2, "module_top_kinase"] == "K1"
    assert assignments.loc[p1s1, "module_top_kinase_candidates"] == ("K1", "K2")
    assert assignments.loc[p1s2, "module_top_kinase_candidates"] == ("K1", "K2")
    assert int(assignments.loc[p1s1, "module_top_kinase_tie_count"]) == 2
    assert int(assignments.loc[p1s2, "module_top_kinase_tie_count"]) == 2
    assert bool(assignments.loc[p1s1, "module_top_kinase_is_ambiguous"])
    assert bool(assignments.loc[p1s2, "module_top_kinase_is_ambiguous"])
    assert (
        assignments.loc[p1s1, "module_top_kinase_selection_policy"]
        == LEXICOGRAPHIC_TIE_BREAK_POLICY
    )
    assert (
        assignments.loc[p1s2, "module_top_kinase_selection_policy"]
        == LEXICOGRAPHIC_TIE_BREAK_POLICY
    )


def test_build_module_assignments_marks_zero_evidence_rows_without_false_winner() -> (
    None
):
    site_index, site_metadata = _site_metadata(["P1;S1;", "P2;S2;"])
    prediction_matrix = pd.DataFrame(
        {
            "K1": [float("nan"), 0.7],
            "K2": [float("nan"), 0.2],
        },
        index=site_index.copy(),
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
        site_to_protein_group_id=site_to_protein,
        site_metadata=site_metadata,
    )

    unsupported = assignments.loc[site_index[0]]
    assert unsupported["top_kinase"] == UNSUPPORTED_KINASE
    assert unsupported["top_kinase_candidates"] == ()
    assert unsupported["top_kinase_weights"] == ()
    assert int(unsupported["top_kinase_tie_count"]) == 0
    assert not bool(unsupported["top_kinase_is_ambiguous"])
    assert unsupported["top_kinase_selection_policy"] == NO_SUPPORT_SELECTION_POLICY
    assert pd.isna(unsupported["top_score"])
    assert int(unsupported["module_id"]) == 0
    assert unsupported["module_top_kinase"] == UNSUPPORTED_KINASE
    assert (
        unsupported["module_top_kinase_selection_policy"] == NO_SUPPORT_SELECTION_POLICY
    )
