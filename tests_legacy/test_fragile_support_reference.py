from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.profiles import build_kinase_substrate_profiles
from phospy.prediction.scoring import combine_profile_and_motif_scores

from phospy.prediction import KinaseScorer

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = ROOT / "tests" / "fixtures" / "fragile_support_reference"


def _read_indexed_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(FIXTURE_DIR / name, index_col=0)


def _read_grouped_mapping(name: str) -> dict[str, list[str]]:
    frame = pd.read_csv(FIXTURE_DIR / name)
    return {
        str(kinase): group.loc[:, "site_id"].astype(str).tolist()
        for kinase, group in frame.groupby("kinase", sort=False)
    }


def test_fragile_support_reference_has_mixed_support_and_candidate_states() -> None:
    summary = pd.read_csv(FIXTURE_DIR / "screening_summary.csv")

    assert set(summary.loc[:, "candidate_status"]) == {
        "dropped",
        "just_above_inclusion",
        "just_below_inclusion",
        "robust",
    }
    assert int(summary.loc[:, "substrate_count"].min()) == 1
    assert int(summary.loc[:, "substrate_count"].max()) >= 8
    assert int(summary.loc[:, "candidate_count"].min()) <= 2
    assert int(summary.loc[:, "candidate_count"].max()) >= 30
    assert summary.loc[
        summary.loc[:, "candidate_status"] == "dropped", "kinase"
    ].tolist() == ["LCK"]
    assert summary.loc[
        summary.loc[:, "candidate_status"] == "just_below_inclusion", "kinase"
    ].tolist() == ["PRKAA1"]
    assert summary.loc[
        summary.loc[:, "candidate_status"] == "just_above_inclusion", "kinase"
    ].tolist() == ["PRKAA2", "AKT1"]
    assert summary.loc[
        summary.loc[:, "candidate_status"] == "robust", "kinase"
    ].tolist() == [
        "MAPK1",
        "IRAK1",
    ]


def test_fragile_support_reference_recomputes_deterministic_seam_outputs() -> None:
    phospho_matrix = _read_indexed_csv("phospho_matrix.csv")
    substrate_map = _read_grouped_mapping("substrate_map.csv")
    expected_profile_matrix = _read_indexed_csv("profile_matrix.csv")
    expected_profile_scores = _read_indexed_csv("profile_scores.csv")
    motif_scores = _read_indexed_csv("motif_scores.csv")
    motif_sizes = pd.read_csv(FIXTURE_DIR / "motif_sizes.csv").set_index("kinase")[
        "motif_size"
    ]
    expected_combined_scores = _read_indexed_csv("combined_scores.csv")
    expected_combined_weights = (
        pd.read_csv(FIXTURE_DIR / "combined_weights.csv")
        .set_index("kinase")
        .sort_index()
    )
    expected_candidate_substrates = _read_grouped_mapping("candidate_substrates.csv")

    profile_result = build_kinase_substrate_profiles(
        substrate_map=substrate_map,
        phospho_matrix=phospho_matrix,
        min_substrates=1,
    )
    actual_profile_matrix = profile_result.profile_matrix.sort_index().sort_index(
        axis=1
    )
    actual_profile_scores = (
        KinaseScorer(actual_profile_matrix)
        .score_phosphosite_profiles(phospho_matrix)
        .sort_index()
        .sort_index(axis=1)
    )
    actual_combined_scores, actual_combined_weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=actual_profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_result.substrate_counts.astype(float),
    )
    actual_combined_scores = actual_combined_scores.sort_index().sort_index(axis=1)
    actual_combined_weights = actual_combined_weights.sort_index()
    actual_candidate_substrates = build_candidate_substrate_list(
        combined_scores=actual_combined_scores,
        top=50,
        score_threshold=0.8,
        inclusion=20,
    )

    pdt.assert_frame_equal(
        actual_profile_matrix,
        expected_profile_matrix.sort_index().sort_index(axis=1),
        check_dtype=False,
        check_names=False,
        atol=1e-8,
        rtol=1e-6,
    )
    pdt.assert_frame_equal(
        actual_profile_scores,
        expected_profile_scores.sort_index().sort_index(axis=1),
        check_dtype=False,
        check_names=False,
        atol=1e-8,
        rtol=1e-6,
    )
    pdt.assert_frame_equal(
        actual_combined_scores,
        expected_combined_scores.sort_index().sort_index(axis=1),
        check_dtype=False,
        check_names=False,
        atol=1e-8,
        rtol=1e-6,
    )
    pdt.assert_frame_equal(
        actual_combined_weights,
        expected_combined_weights,
        check_dtype=False,
        check_names=False,
        atol=1e-8,
        rtol=1e-6,
    )
    assert actual_candidate_substrates == expected_candidate_substrates
