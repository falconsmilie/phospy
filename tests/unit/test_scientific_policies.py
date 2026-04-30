from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.scoring import MotifProfileRankFusionPolicy
from phospy.scientific_policies import (
    PROFILE_CORRELATION_SHIFTED_UNIT_POLICY,
    PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY,
    CandidateSubstrateSelectionPolicy,
    KinaseProfileScoringPolicy,
    PreprocessingStageOrderPolicy,
    ScientificPolicyId,
    ScientificPolicyRecord,
    ScorePreconditioningPolicy,
    SignalomeMissingValueClusteringPolicy,
    build_signalome_module_candidate_score_policy,
    build_simplified_weighted_substrate_activity_policy,
    shift_correlation_to_unit_support,
)
from phospy.signalomes.clustering import derive_protein_modules


def test_scientific_policy_ids_are_stable() -> None:
    assert ScientificPolicyId.PROFILE_CORRELATION_SHIFTED_UNIT.value == (
        "profile_correlation_shifted_unit_v1"
    )
    assert ScientificPolicyId.KINASE_PROFILE_SCORING.value == (
        "kinase_profile_scoring_v1"
    )
    assert ScientificPolicyId.MOTIF_PROFILE_RANK_FUSION.value == (
        "motif_profile_rank_fusion_v1"
    )
    assert ScientificPolicyId.CANDIDATE_SUBSTRATE_SELECTION.value == (
        "candidate_substrate_selection_v1"
    )
    assert ScientificPolicyId.SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY.value == (
        "simplified_weighted_substrate_activity_v1"
    )
    assert ScientificPolicyId.SIGNALOME_MISSING_VALUE_CLUSTERING.value == (
        "signalome_missing_value_clustering_v1"
    )
    assert ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING.value == (
        "signalome_score_preconditioning_v1"
    )
    assert ScientificPolicyId.PREPROCESSING_STAGE_ORDER.value == (
        "preprocessing_stage_order_v1"
    )
    assert ScientificPolicyId.SIGNALOME_MODULE_CANDIDATE_SCORE.value == (
        "signalome_module_candidate_score_v1"
    )
    assert ScientificPolicyId.PROTEIN_MODULE_FROM_SITE_MEMBERSHIP.value == (
        "protein_module_from_site_membership_v1"
    )


def test_scientific_policy_record_payload_round_trip() -> None:
    policy = build_simplified_weighted_substrate_activity_policy(
        threshold=0.6,
        min_substrates=3,
        top_n_substrates=20,
    )
    restored = ScientificPolicyRecord.from_payload(policy.to_payload())

    assert restored.id == policy.id
    assert restored.name == policy.name
    assert restored.version == policy.version
    assert restored.parameters == policy.parameters
    assert restored.assumptions == policy.assumptions
    assert restored.output_scale == policy.output_scale


def test_shift_correlation_to_unit_support_respects_nan_and_bounds() -> None:
    correlation = np.asarray([-1.0, -0.5, 0.0, 0.5, 1.0, np.nan, 2.0], dtype=float)
    scores = shift_correlation_to_unit_support(correlation)

    assert scores[:5].tolist() == pytest.approx([0.0, 0.25, 0.5, 0.75, 1.0])
    assert np.isnan(scores[5])
    assert scores[6] == pytest.approx(1.0)
    assert PROFILE_CORRELATION_SHIFTED_UNIT_POLICY.id == (
        ScientificPolicyId.PROFILE_CORRELATION_SHIFTED_UNIT
    )


def test_motif_profile_rank_fusion_policy_exposes_metadata() -> None:
    policy = MotifProfileRankFusionPolicy(
        allow_profile_only_fallback=True,
        emit_weights=False,
    )

    record = policy.record
    assert record.id == ScientificPolicyId.MOTIF_PROFILE_RANK_FUSION
    assert record.parameters["allow_profile_only_fallback"] is True
    assert record.parameters["emit_weights"] is False


def test_signalome_candidate_score_policy_exposes_runtime_parameters() -> None:
    record = build_signalome_module_candidate_score_policy(
        requested_policy="sampled",
        candidate_scoring_policy="sampled",
        candidate_scoring_mode="sampled",
        max_exact_tree_sites=2000,
        max_full_candidate_scoring_sites=2000,
        candidate_scoring_evaluated=True,
        candidate_scoring_skip_reason=None,
    )

    assert record.id == ScientificPolicyId.SIGNALOME_MODULE_CANDIDATE_SCORE
    assert record.parameters["requested_policy"] == "sampled"
    assert record.parameters["candidate_scoring_policy"] == "sampled"
    assert record.parameters["candidate_scoring_mode"] == "sampled"
    assert (
        record.parameters["candidate_scoring_scope"]
        == "candidate_module_count_evaluation_only"
    )
    assert record.parameters["tree_generation_mode"] == "full_exact_tree_construction"
    assert record.parameters["tree_generation_is_approximate"] is False
    assert record.parameters["candidate_scoring_evaluated"] is True


def test_candidate_substrate_selection_policy_matches_selection_behavior() -> None:
    scores = pd.DataFrame(
        {"K1": [0.8, 0.5, 0.2]},
        index=pd.Index(["S1", "S2", "S3"], name="site_id"),
        dtype=float,
    )
    policy = CandidateSubstrateSelectionPolicy(
        top_k=2,
        score_threshold=0.5,
        inclusion=1,
    )
    selected = build_candidate_substrate_list(
        scores=scores,
        top=policy.top_k,
        score_threshold=policy.score_threshold,
        inclusion=policy.inclusion,
    )

    assert selected == {"K1": ["S1"]}
    assert policy.record.id == ScientificPolicyId.CANDIDATE_SUBSTRATE_SELECTION
    assert policy.record.parameters["top_k"] == 2
    assert policy.record.parameters["score_threshold"] == pytest.approx(0.5)
    assert policy.record.parameters["inclusion"] == 1


def test_kinase_profile_scoring_policy_exposes_self_inclusion_semantics() -> None:
    policy = KinaseProfileScoringPolicy(
        profile_missing_value_strategy="strict",
        min_substrates_floor=2,
        requested_min_substrates=3,
    )

    assert policy.record.id == ScientificPolicyId.KINASE_PROFILE_SCORING
    assert policy.record.parameters["self_inclusion_behavior"] == "self_inclusion"
    assert policy.record.parameters["leave_one_out_enabled"] is False


def test_signalome_missing_value_and_score_preconditioning_policies_are_serializable() -> (
    None
):
    missing_value_policy = SignalomeMissingValueClusteringPolicy(
        missing_value_policy=(
            "column_median_imputation_with_zero_for_all_missing_columns"
        ),
        applies_to="clustering_distance_and_tree_construction_only",
        imputed_values_exposed_in_output_tables=False,
    )
    preconditioning_policy = ScorePreconditioningPolicy(
        policy="allow_and_report",
        input_row_count=100,
        dropped_all_missing_row_count=5,
        retained_row_count=95,
    )
    restored_missing = ScientificPolicyRecord.from_payload(
        missing_value_policy.record.to_payload()
    )
    restored_preconditioning = ScientificPolicyRecord.from_payload(
        preconditioning_policy.record.to_payload()
    )

    assert restored_missing.id == ScientificPolicyId.SIGNALOME_MISSING_VALUE_CLUSTERING
    assert restored_preconditioning.id == (
        ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING
    )
    assert restored_preconditioning.parameters["retained_row_count"] == 95


def test_preprocessing_stage_order_policy_records_stage_order_details() -> None:
    policy = PreprocessingStageOrderPolicy(
        configured_stage_order=("missing_data", "intensity_transform", "site_matrix"),
        default_stage_order=("missing_data",),
        supported_stage_order=(
            "missing_data",
            "intensity_transform",
            "total_protein_correction",
            "site_matrix",
            "normalisation",
            "comparisons",
        ),
    )

    assert policy.record.id == ScientificPolicyId.PREPROCESSING_STAGE_ORDER
    assert policy.record.parameters["configured_stage_order"] == (
        "missing_data -> intensity_transform -> site_matrix"
    )
    assert policy.record.parameters["configured_stage_count"] == 3


def test_protein_module_from_site_membership_policy_matches_derivation_behavior() -> (
    None
):
    site_clusters = pd.Series(
        [1, 2, 1, 2, 3],
        index=pd.Index(["S1", "S2", "S3", "S4", "S5"], name="site_id"),
        dtype="int64",
    )
    site_to_protein = pd.Series(
        ["P1", "P1", "P2", "P2", "P3"],
        index=site_clusters.index.copy(),
        dtype=str,
    )

    protein_modules = derive_protein_modules(
        site_clusters=site_clusters,
        site_to_protein=site_to_protein,
    )

    assert protein_modules.at["P1"] == protein_modules.at["P2"]
    assert protein_modules.at["P3"] != protein_modules.at["P1"]
    assert PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY.id == (
        ScientificPolicyId.PROTEIN_MODULE_FROM_SITE_MEMBERSHIP
    )
