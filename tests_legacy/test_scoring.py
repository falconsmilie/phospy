from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from phospy.prediction.scoring import combine_profile_and_motif_scores

from phospy.prediction import KinaseScorer


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


def test_score_phosphosite_profiles_matches_unbatched_result() -> None:
    scorer = KinaseScorer(make_kinase_profiles(), correlation_batch_size=1)
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 3.0, 2.0, 4.0],
            "sample_2": [2.0, 2.0, 1.0, 4.0],
            "sample_3": [3.0, 1.0, 0.0, 4.0],
        },
        index=["SITE_A", "SITE_B", "SITE_C", "SITE_CONST"],
    )

    chunked = scorer.score_phosphosite_profiles(phospho_matrix)
    unbatched = scorer.score_phosphosite_profiles(
        phospho_matrix,
        correlation_batch_size=None,
    )

    pd.testing.assert_frame_equal(chunked, unbatched)


def test_score_phosphosite_profiles_rejects_invalid_batch_size() -> None:
    scorer = KinaseScorer(make_kinase_profiles())
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0],
            "sample_2": [2.0],
            "sample_3": [3.0],
        },
        index=["SITE_A"],
    )

    with pytest.raises(ValueError, match="correlation_batch_size must be at least 1"):
        scorer.score_phosphosite_profiles(
            phospho_matrix,
            correlation_batch_size=0,
        )


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


def test_combine_profile_and_motif_scores_preserves_profile_only_kinases_on_partial_overlap() -> (
    None
):
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
            "KINASE_C": [0.5, 0.9],
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
        allow_profile_only_fallback=True,
    )

    expected_overlap = pd.DataFrame(
        {
            "KINASE_A": [
                (np.log(2.0) * 0.9 + np.log(3.0) * 0.3) / (np.log(2.0) + np.log(3.0)),
                (np.log(2.0) * 0.1 + np.log(3.0) * 0.7) / (np.log(2.0) + np.log(3.0)),
            ],
            "KINASE_B": [
                (np.log(3.0) * 0.2 + np.log(2.0) * 0.6) / (np.log(3.0) + np.log(2.0)),
                (np.log(3.0) * 0.8 + np.log(2.0) * 0.4) / (np.log(3.0) + np.log(2.0)),
            ],
            "KINASE_C": [0.5, 0.9],
        },
        index=["SITE_1", "SITE_2"],
    )

    pd.testing.assert_frame_equal(combined, expected_overlap)
    assert float(weights.loc["KINASE_C", "motif_weight"]) == pytest.approx(0.0)
    assert float(weights.loc["KINASE_C", "profile_weight"]) == pytest.approx(1.0)
    assert float(weights.loc["KINASE_C", "motif_rank_weight"]) == pytest.approx(0.0)
    assert float(weights.loc["KINASE_C", "profile_rank_weight"]) == pytest.approx(1.0)


def test_kinase_scoring_result_tables_are_detached_from_input_matrix() -> None:
    scorer = KinaseScorer(make_kinase_profiles())
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
            "KINASE_A": [0.9, 0.1],
            "KINASE_B": [0.2, 0.8],
        },
        index=phospho_matrix.index.copy(),
    )
    motif_sizes = pd.Series({"KINASE_A": 10, "KINASE_B": 20}, dtype=float)
    profile_sizes = pd.Series({"KINASE_A": 30, "KINASE_B": 5}, dtype=float)
    original = phospho_matrix.copy(deep=True)

    result = scorer.score(
        phospho_matrix=phospho_matrix,
        motif_scores=motif_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
    )

    result.profile_scores.loc["SITE_A", "KINASE_A"] = -999.0
    assert result.combined_scores is not None
    result.combined_scores.loc["SITE_A", "KINASE_A"] = -999.0
    assert result.weights is not None
    result.weights.loc["KINASE_A", "motif_weight"] = -999.0

    pd.testing.assert_frame_equal(phospho_matrix, original)
