from __future__ import annotations

import importlib

import pandas as pd
import pytest

from phospy.activities.scientific_policies import (
    build_ksea_zscore_activity_policy,
    build_simplified_weighted_substrate_activity_policy,
)
from phospy.datasets.preprocessing.scientific_policies import (
    DUPLICATE_SITE_RESOLUTION_AGGREGATE_MEAN_POLICY,
    PreprocessingStageOrderPolicy,
    build_duplicate_site_resolution_policy,
)
from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.policies import resolve_prediction_sampling_policy
from phospy.prediction.scientific_policies import (
    PROFILE_CORRELATION_SHIFTED_UNIT_POLICY,
    CandidateSubstrateSelectionPolicy,
    KinaseProfileScoringPolicy,
)
from phospy.prediction.scoring import (
    MotifProfileRankFusionPolicy,
    resolve_downstream_score_matrix,
)
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.signalomes.clustering import derive_protein_modules
from phospy.signalomes.clustering.policies import (
    resolve_candidate_scoring_policy_definition,
)
from phospy.signalomes.clustering.scientific_policies import (
    PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY,
    SignalomeMissingValueClusteringPolicy,
    build_signalome_module_candidate_score_policy,
)
from phospy.workflows.signalome.scientific_policies import (
    ScorePreconditioningPolicy,
    resolve_score_preconditioning_policy,
)


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
    assert ScientificPolicyId.KSEA_ZSCORE_ACTIVITY.value == "ksea_zscore_activity_v1"
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
    assert ScientificPolicyId.DUPLICATE_SITE_RESOLUTION.value == (
        "duplicate_site_resolution_v1"
    )
    assert ScientificPolicyId.ADAPTIVE_PREDICTION_SAMPLING.value == (
        "adaptive_prediction_sampling_v1"
    )
    assert ScientificPolicyId.SIGNALOME_DOWNSTREAM_SCORE_SELECTION.value == (
        "signalome_downstream_score_selection_v1"
    )
    assert ScientificPolicyId.SIGNALOME_CANDIDATE_SCORING.value == (
        "signalome_candidate_scoring_v1"
    )


def test_scientific_policy_model_ownership_is_explicit() -> None:
    assert ScientificPolicyId.__module__ == "phospy.provenance.scientific_policy_models"
    assert (
        ScientificPolicyRecord.__module__
        == "phospy.provenance.scientific_policy_models"
    )


def test_prediction_scientific_policy_ownership_is_explicit() -> None:
    assert (
        KinaseProfileScoringPolicy.__module__ == "phospy.prediction.scientific_policies"
    )
    assert (
        CandidateSubstrateSelectionPolicy.__module__
        == "phospy.prediction.scientific_policies"
    )


def test_activity_scientific_policy_ownership_is_explicit() -> None:
    assert (
        build_ksea_zscore_activity_policy.__module__
        == "phospy.activities.scientific_policies"
    )
    assert (
        build_simplified_weighted_substrate_activity_policy.__module__
        == "phospy.activities.scientific_policies"
    )


def test_preprocessing_scientific_policy_ownership_is_explicit() -> None:
    assert (
        PreprocessingStageOrderPolicy.__module__
        == "phospy.datasets.preprocessing.scientific_policies"
    )
    assert (
        build_duplicate_site_resolution_policy.__module__
        == "phospy.datasets.preprocessing.scientific_policies"
    )


def test_signalome_workflow_scientific_policy_ownership_is_explicit() -> None:
    assert (
        ScorePreconditioningPolicy.__module__
        == "phospy.workflows.signalome.scientific_policies"
    )
    assert (
        resolve_score_preconditioning_policy.__module__
        == "phospy.workflows.signalome.scientific_policies"
    )


def test_signalome_clustering_scientific_policy_ownership_is_explicit() -> None:
    assert (
        SignalomeMissingValueClusteringPolicy.__module__
        == "phospy.signalomes.clustering.scientific_policies"
    )
    assert (
        build_signalome_module_candidate_score_policy.__module__
        == "phospy.signalomes.clustering.scientific_policies"
    )


def test_root_scientific_policies_module_is_removed() -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("phospy.scientific_policies")


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


def test_ksea_activity_policy_record_is_serializable() -> None:
    policy = build_ksea_zscore_activity_policy(
        evidence_threshold=0.6,
        min_substrates=5,
        p_value_method="normal_approximation",
        adjust_p_values=True,
        q_value_method="benjamini_hochberg",
    )
    restored = ScientificPolicyRecord.from_payload(policy.to_payload())
    assert restored.id == ScientificPolicyId.KSEA_ZSCORE_ACTIVITY
    assert restored.parameters["evidence_threshold"] == pytest.approx(0.6)
    assert restored.parameters["min_substrates"] == 5


def test_profile_correlation_shifted_unit_policy_record_is_stable() -> None:
    assert PROFILE_CORRELATION_SHIFTED_UNIT_POLICY.id == (
        ScientificPolicyId.PROFILE_CORRELATION_SHIFTED_UNIT
    )
    assert PROFILE_CORRELATION_SHIFTED_UNIT_POLICY.parameters["transform"] == (
        "(r + 1) / 2"
    )
    assert (
        PROFILE_CORRELATION_SHIFTED_UNIT_POLICY.parameters["clip_to_unit_interval"]
        is True
    )
    assert (
        PROFILE_CORRELATION_SHIFTED_UNIT_POLICY.parameters["preserve_undefined_as_nan"]
        is True
    )
    assert PROFILE_CORRELATION_SHIFTED_UNIT_POLICY.quantitative_meaning == (
        "relative_support_score"
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


def test_public_config_modes_resolve_to_expected_versioned_policy_objects() -> None:
    stable_sampling = resolve_prediction_sampling_policy("stable")
    parity_sampling = resolve_prediction_sampling_policy("r_parity")
    sampled_candidate = resolve_candidate_scoring_policy_definition(
        candidate_scoring_policy="sampled"
    )
    strict_preconditioning = resolve_score_preconditioning_policy(
        policy="error_on_drop"
    )
    duplicate_aggregate_mean = build_duplicate_site_resolution_policy(
        duplicate_site_policy="aggregate_mean"
    )

    assert stable_sampling.name == "adaptive_prediction_sampling_stable_v1"
    assert parity_sampling.name == "adaptive_prediction_sampling_r_parity_v1"
    assert sampled_candidate.name == "signalome_candidate_scoring_sampled_v1"
    assert strict_preconditioning.name == "score_preconditioning_error_on_drop_v1"
    assert (
        duplicate_aggregate_mean.name
        == DUPLICATE_SITE_RESOLUTION_AGGREGATE_MEAN_POLICY.name
    )


def test_downstream_score_resolution_returns_rank_weighted_preferred_policy() -> None:
    profile_scores = pd.DataFrame({"K1": [0.1, 0.2]}, index=["S1", "S2"])
    rank_weighted = pd.DataFrame({"K1": [0.7, 0.8]}, index=["S1", "S2"])
    _selected, source, policy = resolve_downstream_score_matrix(
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted,
    )

    assert source == "rank_weighted_fusion_scores"
    assert policy.name == "signalome_downstream_score_rank_weighted_preferred_v1"


def test_policy_parameters_are_immutable_mappings() -> None:
    policy = resolve_score_preconditioning_policy(policy="allow_and_report")
    with pytest.raises(TypeError):
        policy.parameters["new_key"] = "new_value"  # type: ignore[index]
