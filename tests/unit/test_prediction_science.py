from __future__ import annotations

import pandas as pd
import pytest

from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.scoring import (
    DOWNSTREAM_SCORE_SOURCE_PROFILE,
    DOWNSTREAM_SCORE_SOURCE_RANK_WEIGHTED_FUSION,
    KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE,
    KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT,
    KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP,
    KINASE_SCORE_SOURCE_SUMMARY_COLUMNS,
    KINASE_SCORE_SOURCE_UNAVAILABLE_NO_SCORE,
    MotifProfileRankFusionPolicy,
    build_kinase_score_source_diagnostics,
    fuse_profile_and_motif_scores_by_rank_weight,
    select_downstream_score_matrix,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyId


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


def test_build_kinase_score_source_diagnostics_tracks_fused_and_profile_fallback() -> (
    None
):
    profile_scores = pd.DataFrame(
        {"K1": [0.9, 0.6, float("nan")]},
        index=["S1", "S2", "S3"],
    )
    motif_scores = pd.DataFrame(
        {"K1": [0.2, float("nan"), 0.5]},
        index=profile_scores.index.copy(),
    )
    rank_weighted_fusion_scores = pd.DataFrame(
        {"K1": [0.55, 0.6, float("nan")]},
        index=profile_scores.index.copy(),
    )

    source_matrix, source_summary = build_kinase_score_source_diagnostics(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted_fusion_scores,
    )

    assert (
        source_matrix.at["S1", "K1"] == KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE
    )
    assert source_matrix.at["S2", "K1"] == (
        KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT
    )
    assert source_matrix.at["S3", "K1"] == KINASE_SCORE_SOURCE_UNAVAILABLE_NO_SCORE
    assert (
        tuple(source_summary.columns.astype(str)) == KINASE_SCORE_SOURCE_SUMMARY_COLUMNS
    )
    assert source_summary.at["K1", "fused_motif_profile_evidence_count"] == 1
    assert source_summary.at["K1", "profile_only_motif_missing_or_constant_count"] == 1
    assert source_summary.at["K1", "profile_only_no_motif_overlap_count"] == 0
    assert source_summary.at["K1", "unavailable_no_score_count"] == 1
    assert source_summary.at["K1", "sites_with_score_count"] == 2
    assert source_summary.at["K1", "total_sites_count"] == 3


def test_build_kinase_score_source_diagnostics_tracks_profile_only_no_overlap() -> None:
    profile_scores = pd.DataFrame(
        {"K_PROFILE_ONLY": [0.3, float("nan")]},
        index=["S1", "S2"],
    )
    motif_scores = pd.DataFrame(index=profile_scores.index.copy())
    rank_weighted_fusion_scores = profile_scores.copy(deep=True)

    source_matrix, source_summary = build_kinase_score_source_diagnostics(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted_fusion_scores,
    )

    assert source_matrix.at["S1", "K_PROFILE_ONLY"] == (
        KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP
    )
    assert source_matrix.at["S2", "K_PROFILE_ONLY"] == (
        KINASE_SCORE_SOURCE_UNAVAILABLE_NO_SCORE
    )
    assert (
        source_summary.at["K_PROFILE_ONLY", "profile_only_no_motif_overlap_count"] == 1
    )
    assert source_summary.at["K_PROFILE_ONLY", "unavailable_no_score_count"] == 1


def test_build_kinase_score_source_diagnostics_marks_all_nan_motif_column_as_profile_only() -> (
    None
):
    profile_scores = pd.DataFrame({"K1": [0.4, 0.7]}, index=["S1", "S2"])
    motif_scores = pd.DataFrame(
        {"K1": [float("nan"), float("nan")]},
        index=profile_scores.index.copy(),
    )
    rank_weighted_fusion_scores = profile_scores.copy(deep=True)

    source_matrix, source_summary = build_kinase_score_source_diagnostics(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted_fusion_scores,
    )

    assert (
        source_matrix.loc[:, "K1"]
        == KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT
    ).all()
    assert source_summary.at["K1", "fused_motif_profile_evidence_count"] == 0
    assert source_summary.at["K1", "profile_only_motif_missing_or_constant_count"] == 2
