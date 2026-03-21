from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phosrpy.scoring import (
    KinaseScorer,
    KinaseSubstrateScoreResult,
    combine_profile_and_motif_scores,
    kinase_substrate_score,
)


def make_kinase_profiles() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [1.0, 3.0],
            "sample_2": [2.0, 2.0],
            "sample_3": [3.0, 1.0],
        },
        index=["KINASE_A", "KINASE_B"],
    )


def test_score_phosphosite_profiles_rescales_correlations() -> None:
    scorer = KinaseScorer(make_kinase_profiles())
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 3.0, 1.0],
            "sample_2": [2.0, 2.0, 1.0],
            "sample_3": [3.0, 1.0, 1.0],
        },
        index=["SITE_A", "SITE_B", "SITE_CONST"],
    )

    result = scorer.score_phosphosite_profiles(phospho_matrix)

    assert float(result.loc["SITE_A", "KINASE_A"]) == pytest.approx(1.0)
    assert float(result.loc["SITE_A", "KINASE_B"]) == pytest.approx(0.0)
    assert float(result.loc["SITE_B", "KINASE_A"]) == pytest.approx(0.0)
    assert float(result.loc["SITE_B", "KINASE_B"]) == pytest.approx(1.0)
    assert np.isnan(result.loc["SITE_CONST", "KINASE_A"])
    assert np.isnan(result.loc["SITE_CONST", "KINASE_B"])


def test_score_phosphosite_profiles_aligns_columns_by_name() -> None:
    scorer = KinaseScorer(make_kinase_profiles())
    phospho_matrix = pd.DataFrame(
        {
            "sample_3": [3.0],
            "sample_1": [1.0],
            "sample_2": [2.0],
        },
        index=["SITE_A"],
    )

    result = scorer.score_phosphosite_profiles(phospho_matrix)
    assert float(result.loc["SITE_A", "KINASE_A"]) == pytest.approx(1.0)


def test_score_phosphosite_profiles_requires_matching_columns() -> None:
    scorer = KinaseScorer(make_kinase_profiles())
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0],
            "sample_2": [2.0],
            "sample_x": [3.0],
        },
        index=["SITE_A"],
    )

    with pytest.raises(ValueError, match="must match kinase profile columns"):
        scorer.score_phosphosite_profiles(phospho_matrix)


def test_combine_profile_and_motif_scores_uses_rank_weights() -> None:
    motif_scores = pd.DataFrame(
        {
            "KINASE_A": [0.9, 0.1],
            "KINASE_B": [0.2, 0.8],
        },
        index=["SITE_1", "SITE_2"],
    )
    profile_scores = pd.DataFrame(
        {
            "KINASE_B": [0.6, 0.4],
            "KINASE_A": [0.3, 0.7],
            "KINASE_C": [0.5, 0.5],
        },
        index=["SITE_1", "SITE_2"],
    )
    motif_sizes = pd.Series({"KINASE_A": 10, "KINASE_B": 20})
    profile_sizes = pd.Series({"KINASE_A": 30, "KINASE_B": 5, "KINASE_C": 15})

    combined, weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
    )

    expected_a = (np.log(2.0) * 0.9 + np.log(3.0) * 0.3) / (np.log(2.0) + np.log(3.0))
    expected_b = (np.log(3.0) * 0.2 + np.log(2.0) * 0.6) / (np.log(3.0) + np.log(2.0))

    assert list(combined.columns) == ["KINASE_A", "KINASE_B"]
    assert float(combined.loc["SITE_1", "KINASE_A"]) == pytest.approx(expected_a)
    assert float(combined.loc["SITE_1", "KINASE_B"]) == pytest.approx(expected_b)
    assert float(weights.loc["KINASE_A", "motif_weight"]) == pytest.approx(
        np.log(2.0) / (np.log(2.0) + np.log(3.0))
    )
    assert float(weights.loc["KINASE_B", "profile_weight"]) == pytest.approx(
        np.log(2.0) / (np.log(2.0) + np.log(3.0))
    )


def test_combine_profile_and_motif_scores_can_fall_back_to_profile_only() -> None:
    motif_scores = pd.DataFrame(
        {"KINASE_X": [0.9]},
        index=["SITE_1"],
    )
    profile_scores = pd.DataFrame(
        {"KINASE_A": [0.4], "KINASE_B": [0.6]},
        index=["SITE_1"],
    )
    motif_sizes = pd.Series({"KINASE_X": 10})
    profile_sizes = pd.Series({"KINASE_A": 2, "KINASE_B": 3})

    combined, weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
        allow_profile_only_fallback=True,
    )

    pd.testing.assert_frame_equal(combined, profile_scores)
    assert (weights["motif_weight"] == 0.0).all()
    assert (weights["profile_weight"] == 1.0).all()


def test_kinase_substrate_score_builds_profile_only_result() -> None:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 3.0],
            "sample_2": [2.0, 2.0],
            "sample_3": [3.0, 1.0],
        },
        index=["SITE_A", "SITE_B"],
    )

    result = kinase_substrate_score(
        substrate_map={"KINASE_A": ["SITE_A"], "KINASE_B": ["SITE_B"]},
        phospho_matrix=phospho_matrix,
        allow_profile_only_fallback=True,
    )

    assert isinstance(result, KinaseSubstrateScoreResult)
    assert result.motif_scores is None
    pd.testing.assert_frame_equal(result.profile_scores, result.combined_scores)
    assert list(result.ks_activity_matrix.index) == ["KINASE_A", "KINASE_B"]
    assert float(result.profile_scores.loc["SITE_A", "KINASE_A"]) == pytest.approx(1.0)


def test_kinase_substrate_score_combines_filtered_motif_and_profile_inputs() -> None:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 3.0],
            "sample_2": [2.0, 2.0],
            "sample_3": [3.0, 1.0],
        },
        index=["SITE_A", "SITE_B"],
    )
    motif_scores = pd.DataFrame(
        {
            "KINASE_A": [0.8, 0.3],
            "KINASE_B": [0.4, 0.9],
            "KINASE_X": [0.7, 0.6],
        },
        index=phospho_matrix.index.copy(),
    )
    motif_sizes = pd.Series({"KINASE_A": 4, "KINASE_B": 2, "KINASE_X": 1})

    result = kinase_substrate_score(
        substrate_map={"KINASE_A": ["SITE_A"], "KINASE_B": ["SITE_B"]},
        phospho_matrix=phospho_matrix,
        motif_scores=motif_scores,
        motif_sizes=motif_sizes,
        min_motif_size=2,
    )

    assert list(result.motif_scores.columns) == ["KINASE_A", "KINASE_B"]
    assert list(result.combined_scores.columns) == ["KINASE_A", "KINASE_B"]
    assert list(result.ks_activity_matrix.index) == ["KINASE_A", "KINASE_B"]
    expected_profile_b = np.log(2.5) / (np.log(2.0) + np.log(2.5))
    assert float(result.weights.loc["KINASE_A", "motif_weight"]) == pytest.approx(
        np.log(3.0) / (np.log(3.0) + np.log(2.5))
    )
    assert float(result.weights.loc["KINASE_B", "profile_weight"]) == pytest.approx(
        expected_profile_b
    )


def test_kinase_substrate_score_can_fall_back_when_no_motif_overlap_exists() -> None:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 3.0],
            "sample_2": [2.0, 2.0],
            "sample_3": [3.0, 1.0],
        },
        index=["SITE_A", "SITE_B"],
    )
    motif_scores = pd.DataFrame(
        {"KINASE_X": [0.8, 0.2]},
        index=phospho_matrix.index.copy(),
    )
    motif_sizes = pd.Series({"KINASE_X": 5})

    result = kinase_substrate_score(
        substrate_map={"KINASE_A": ["SITE_A"], "KINASE_B": ["SITE_B"]},
        phospho_matrix=phospho_matrix,
        motif_scores=motif_scores,
        motif_sizes=motif_sizes,
        allow_profile_only_fallback=True,
    )

    pd.testing.assert_frame_equal(result.combined_scores, result.profile_scores)
    assert list(result.ks_activity_matrix.index) == ["KINASE_A", "KINASE_B"]
