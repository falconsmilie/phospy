from __future__ import annotations

import pandas as pd
import pytest

from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.scoring import (
    DOWNSTREAM_SCORE_SOURCE_PROFILE,
    DOWNSTREAM_SCORE_SOURCE_RANK_WEIGHTED_FUSION,
    MotifProfileRankFusionPolicy,
    fuse_profile_and_motif_scores_by_rank_weight,
    select_downstream_score_matrix,
)
from phospy.scientific_policies import ScientificPolicyId


def test_fuse_profile_and_motif_scores_by_rank_weight_falls_back_when_motif_column_is_all_nan() -> (
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

    rank_weighted_fusion_scores, _score_fusion_weights = (
        fuse_profile_and_motif_scores_by_rank_weight(
            motif_scores=motif_scores,
            profile_scores=profile_scores,
            motif_sizes=motif_sizes,
            profile_sizes=profile_sizes,
        )
    )

    assert rank_weighted_fusion_scores.loc[:, "K1"].tolist() == pytest.approx(
        profile_scores.loc[:, "K1"].tolist()
    )


def test_fuse_profile_and_motif_scores_by_rank_weight_preserves_profile_cell_when_motif_is_missing() -> (
    None
):
    profile_scores = pd.DataFrame({"K1": [0.8, 0.6]}, index=["S1", "S2"])
    motif_scores = pd.DataFrame({"K1": [float("nan"), 0.2]}, index=["S1", "S2"])
    profile_sizes = pd.Series({"K1": 10.0})
    motif_sizes = pd.Series({"K1": 10.0})

    rank_weighted_fusion_scores, _score_fusion_weights = (
        fuse_profile_and_motif_scores_by_rank_weight(
            motif_scores=motif_scores,
            profile_scores=profile_scores,
            motif_sizes=motif_sizes,
            profile_sizes=profile_sizes,
        )
    )

    assert rank_weighted_fusion_scores.at["S1", "K1"] == pytest.approx(
        profile_scores.at["S1", "K1"]
    )
    assert rank_weighted_fusion_scores.at["S2", "K1"] == pytest.approx(0.4)


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


def test_select_downstream_score_matrix_prefers_rank_weighted_fusion_scores() -> None:
    profile_scores = pd.DataFrame({"K1": [0.1, 0.2]}, index=["S1", "S2"])
    rank_weighted_fusion_scores = pd.DataFrame({"K1": [0.7, 0.6]}, index=["S1", "S2"])

    selected, source = select_downstream_score_matrix(
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted_fusion_scores,
    )

    assert selected is rank_weighted_fusion_scores
    assert source == DOWNSTREAM_SCORE_SOURCE_RANK_WEIGHTED_FUSION


def test_select_downstream_score_matrix_falls_back_to_profile_scores() -> None:
    profile_scores = pd.DataFrame({"K1": [0.1, 0.2]}, index=["S1", "S2"])

    selected, source = select_downstream_score_matrix(
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=None,
    )

    assert selected is profile_scores
    assert source == DOWNSTREAM_SCORE_SOURCE_PROFILE


def test_fuse_profile_and_motif_scores_by_rank_weight_can_skip_weight_table() -> None:
    profile_scores = pd.DataFrame({"K1": [0.2, 0.8]}, index=["S1", "S2"])
    motif_scores = pd.DataFrame({"K1": [0.7, 0.3]}, index=["S1", "S2"])
    profile_sizes = pd.Series({"K1": 4.0})
    motif_sizes = pd.Series({"K1": 4.0})

    rank_weighted_fusion_scores, score_fusion_weights = (
        fuse_profile_and_motif_scores_by_rank_weight(
            motif_scores=motif_scores,
            profile_scores=profile_scores,
            motif_sizes=motif_sizes,
            profile_sizes=profile_sizes,
            emit_weights=False,
        )
    )

    assert score_fusion_weights is None
    assert list(rank_weighted_fusion_scores.columns) == ["K1"]


def test_fuse_profile_and_motif_scores_policy_metadata_and_default_behavior() -> None:
    profile_scores = pd.DataFrame({"K1": [0.3, 0.7]}, index=["S1", "S2"])
    motif_scores = pd.DataFrame({"K1": [0.8, 0.2]}, index=["S1", "S2"])
    profile_sizes = pd.Series({"K1": 12.0})
    motif_sizes = pd.Series({"K1": 8.0})

    wrapper_scores, wrapper_weights = fuse_profile_and_motif_scores_by_rank_weight(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
    )
    policy = MotifProfileRankFusionPolicy()
    policy_scores, policy_weights = policy.fuse(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
    )

    pd.testing.assert_frame_equal(wrapper_scores, policy_scores)
    assert wrapper_weights is not None
    assert policy_weights is not None
    pd.testing.assert_frame_equal(wrapper_weights, policy_weights)
    assert policy.record.id == ScientificPolicyId.MOTIF_PROFILE_RANK_FUSION
