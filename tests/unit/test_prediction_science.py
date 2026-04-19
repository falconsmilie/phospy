from __future__ import annotations

import pandas as pd
import pytest

from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.scoring import (
    DOWNSTREAM_SCORE_SOURCE_COMBINED,
    DOWNSTREAM_SCORE_SOURCE_PROFILE,
    combine_profile_and_motif_scores,
    select_downstream_score_matrix,
)


def test_combine_profile_and_motif_scores_falls_back_when_motif_column_is_all_nan() -> (
    None
):
    profile_scores = pd.DataFrame(
        {"K1": [0.8, 0.3], "K2": [0.2, 0.7]},
        index=["S1", "S2"],
    )
    motif_scores = pd.DataFrame(
        {"K1": [float("nan"), float("nan")], "K2": [0.1, 0.9]},
        index=profile_scores.index.copy(),
    )
    profile_sizes = pd.Series({"K1": 20.0, "K2": 25.0})
    motif_sizes = pd.Series({"K1": 12.0, "K2": 14.0})

    combined, _weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
    )

    assert combined.loc[:, "K1"].tolist() == pytest.approx(
        profile_scores.loc[:, "K1"].tolist()
    )


def test_build_candidate_substrate_list_can_restrict_sites_per_kinase() -> None:
    scores = pd.DataFrame(
        {"K1": [0.9, 0.8, 0.7], "K2": [0.2, 0.95, 0.85]},
        index=["S1", "S2", "S3"],
    )

    candidates = build_candidate_substrate_list(
        scores=scores,
        top=3,
        score_threshold=0.0,
        inclusion=1,
        allowed_sites_by_kinase={"K1": ["S1", "S3"], "K2": ["S2"]},
    )

    assert candidates == {"K1": ["S1", "S3"], "K2": ["S2"]}


def test_select_downstream_score_matrix_prefers_combined_scores() -> None:
    profile_scores = pd.DataFrame({"K1": [0.1, 0.2]}, index=["S1", "S2"])
    combined_scores = pd.DataFrame({"K1": [0.7, 0.6]}, index=["S1", "S2"])

    selected, source = select_downstream_score_matrix(
        profile_scores=profile_scores,
        combined_scores=combined_scores,
    )

    assert selected is combined_scores
    assert source == DOWNSTREAM_SCORE_SOURCE_COMBINED


def test_select_downstream_score_matrix_falls_back_to_profile_scores() -> None:
    profile_scores = pd.DataFrame({"K1": [0.1, 0.2]}, index=["S1", "S2"])

    selected, source = select_downstream_score_matrix(
        profile_scores=profile_scores,
        combined_scores=None,
    )

    assert selected is profile_scores
    assert source == DOWNSTREAM_SCORE_SOURCE_PROFILE


def test_combine_profile_and_motif_scores_can_skip_weight_table() -> None:
    profile_scores = pd.DataFrame({"K1": [0.2, 0.8]}, index=["S1", "S2"])
    motif_scores = pd.DataFrame({"K1": [0.7, 0.3]}, index=["S1", "S2"])
    profile_sizes = pd.Series({"K1": 4.0})
    motif_sizes = pd.Series({"K1": 4.0})

    combined_scores, weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
        emit_weights=False,
    )

    assert weights is None
    assert list(combined_scores.columns) == ["K1"]
