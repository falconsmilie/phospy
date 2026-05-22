from __future__ import annotations

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyPhosphoDataset, SignalomeWorkflow
from phospy.api import (
    Organism,
    SignalomeWorkflowRequest,
)
from phospy.api.configs import (
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.provenance.hashing import hash_table
from phospy.science.activities.models import KinaseActivityInputs, PredMatOverlapSummary
from phospy.science.activities.scoring import compute_activity_from_inputs
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.stages.normalisation import (
    NormalisationStage,
)
from phospy.science.datasets.preprocessing.stages.site_matrix import (
    _apply_duplicate_site_policy,
)
from phospy.science.prediction.scoring import (
    fuse_profile_and_motif_scores_by_rank_weight,
)
from phospy.science.references.models import ReferenceBundle
from phospy.science.signalomes.assignments import build_module_assignments
from phospy.science.signalomes.clustering.candidate_scoring import (
    compute_candidate_cluster_scores,
    summarize_profile_degeneracy,
)
from phospy.science.signalomes.clustering.policies import (
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_TREE_ENGINE_EXACT,
)
from phospy.science.signalomes.clustering.tree_building import (
    prepare_scoring_values_for_clustering,
)
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    rank_kinases_for_prediction,
    score_profile_correlations,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config


def _normalisation_state(*, phospho: pd.DataFrame, policy: str) -> PreprocessingState:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C", "D"][: len(phospho.index)],
            "site": ["S1", "S2", "S3", "S4"][: len(phospho.index)],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C", "SEQ_D"][: len(phospho.index)],
        },
        index=phospho.index.copy(),
    )
    return PreprocessingState(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            normalisation_policy=policy,
            stage_order=("normalisation",),
        ),
    )


def _analysis_ready_dataset(site_ids: list[str]) -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0 + index for index, _ in enumerate(site_ids)],
            "sample_b": [2.0 + index for index, _ in enumerate(site_ids)],
        },
        index=pd.Index(site_ids, name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [site_id.split(";", 1)[0] for site_id in site_ids],
            "site": [f"S{index + 1}" for index, _ in enumerate(site_ids)],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in [f"S{index + 1}" for index, _ in enumerate(site_ids)]
            ],
            "protein_id": [f"P{index + 1}" for index, _ in enumerate(site_ids)],
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _reference_bundle(site_ids: list[str]) -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K2"],
                "substrate_site": [site_ids[0], site_ids[1]],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31 for _ in site_ids]},
            index=pd.Index(site_ids, name="site_id"),
        ),
    )


def _kinase_result(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    prediction_matrix: pd.DataFrame,
    score_matrix: pd.DataFrame,
) -> KinaseWorkflowResult:
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_reference_bundle(dataset.phospho.index.astype(str).tolist()),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def test_normalisation_noop_preserves_values_and_hash() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [3.0, 4.0]},
        index=pd.Index(["A;S1;", "B;S2;"], name="site_id"),
    )
    state = _normalisation_state(phospho=phospho, policy="none")
    result = NormalisationStage().run(state)

    assert result.state is state
    pdt.assert_frame_equal(result.state.phospho, phospho)
    assert hash_table(phospho, name="phospho") == hash_table(
        result.state.phospho, name="phospho"
    )


def test_normalisation_noop_records_method_none() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [3.0, 4.0]},
        index=pd.Index(["A;S1;", "B;S2;"], name="site_id"),
    )
    result = NormalisationStage().run(
        _normalisation_state(phospho=phospho, policy="none")
    )

    diagnostics = result.diagnostics["diagnostics"]
    assert diagnostics["method"] == "none"
    assert diagnostics["parameters"] == {"applied": False}
    assert diagnostics["input_matrix_shape"] == {"rows": 2, "columns": 2}
    assert diagnostics["output_matrix_shape"] == {"rows": 2, "columns": 2}


def test_median_centering_records_before_after_sample_summaries() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [4.0, 2.0, 0.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;"], name="site_id"),
    )
    result = NormalisationStage().run(
        _normalisation_state(phospho=phospho, policy="median_center")
    )

    diagnostics = result.diagnostics["diagnostics"]
    before = diagnostics["per_sample_summary_before"]["sample_a"]
    after = diagnostics["per_sample_summary_after"]["sample_a"]
    assert before["median"] == pytest.approx(2.0)
    assert after["median"] == pytest.approx(0.0)
    assert diagnostics["rows_dropped"] is False
    assert diagnostics["columns_dropped"] is False


def test_normalisation_row_reordering_is_semantically_invariant() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [5.0, 2.0, 3.0, 4.0],
            "sample_b": [4.0, 1.0, 2.5, 2.0],
            "sample_c": [8.0, 6.0, 7.0, 9.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;", "D;S4;"], name="site_id"),
    )
    reordered = phospho.iloc[[2, 0, 3, 1], :]

    base = (
        NormalisationStage()
        .run(_normalisation_state(phospho=phospho, policy="quantile"))
        .state.phospho
    )
    shuffled = (
        NormalisationStage()
        .run(_normalisation_state(phospho=reordered, policy="quantile"))
        .state.phospho
    )

    pdt.assert_frame_equal(shuffled.loc[phospho.index, phospho.columns], base)


def test_normalisation_preserves_column_order_for_supported_policies() -> None:
    phospho = pd.DataFrame(
        {
            "sample_c": [10.0, 20.0, 30.0],
            "sample_a": [3.0, 5.0, 7.0],
            "sample_b": [9.0, 1.0, 4.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;"], name="site_id"),
    )

    median_centered = (
        NormalisationStage()
        .run(_normalisation_state(phospho=phospho, policy="median_center"))
        .state.phospho
    )
    quantile = (
        NormalisationStage()
        .run(_normalisation_state(phospho=phospho, policy="quantile"))
        .state.phospho
    )

    assert list(median_centered.columns) == list(phospho.columns)
    assert list(quantile.columns) == list(phospho.columns)


def test_median_centering_matches_hand_calculated_medians() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [4.0, 2.0, 0.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;"], name="site_id"),
    )
    centered = (
        NormalisationStage()
        .run(_normalisation_state(phospho=phospho, policy="median_center"))
        .state.phospho
    )

    expected = pd.DataFrame(
        {
            "sample_a": [-1.0, 0.0, 1.0],
            "sample_b": [2.0, 0.0, -2.0],
        },
        index=phospho.index.copy(),
    )
    pdt.assert_frame_equal(centered, expected)
    assert centered.median(axis=0).tolist() == pytest.approx([0.0, 0.0])


def test_quantile_normalisation_equalises_sample_distributions() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [5.0, 2.0, 3.0, 4.0],
            "sample_b": [4.0, 1.0, 2.5, 2.0],
            "sample_c": [8.0, 6.0, 7.0, 9.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;", "D;S4;"], name="site_id"),
    )
    quantile = (
        NormalisationStage()
        .run(_normalisation_state(phospho=phospho, policy="quantile"))
        .state.phospho
    )

    expected_sorted_distribution = [3.0, 4.0, 4.833333333333333, 6.0]
    for column in quantile.columns:
        assert sorted(quantile.loc[:, column].tolist()) == pytest.approx(
            expected_sorted_distribution
        )


def test_quantile_normalisation_preserves_matrix_shape() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [5.0, 2.0, 3.0, 4.0],
            "sample_b": [4.0, 1.0, 2.5, 2.0],
            "sample_c": [8.0, 6.0, 7.0, 9.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;", "D;S4;"], name="site_id"),
    )
    result = NormalisationStage().run(
        _normalisation_state(phospho=phospho, policy="quantile")
    )

    assert result.state.phospho.shape == phospho.shape
    diagnostics = result.diagnostics["diagnostics"]
    assert diagnostics["input_matrix_shape"] == {"rows": 4, "columns": 3}
    assert diagnostics["output_matrix_shape"] == {"rows": 4, "columns": 3}


def test_normalisation_stage_does_not_alter_site_metadata_index() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [4.0, 2.0, 0.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;"], name="site_id"),
    )
    state = _normalisation_state(phospho=phospho, policy="median_center")
    before_index = state.site_metadata.index.copy()
    result = NormalisationStage().run(state)

    assert result.state.site_metadata.index.equals(before_index)


def test_duplicate_site_resolution_is_deterministic_for_max_mean_signal_policy() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 9.0, 5.0],
            "sample_b": [1.0, 9.0, 6.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "protein_id": ["P1", "P1", "P3"],
        },
        index=phospho.index.copy(),
    )
    site_ids = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;", "AKT1;T308;"],
        index=phospho.index.copy(),
        name="site_id",
    )

    first = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=site_ids,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
    )
    second = _apply_duplicate_site_policy(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=site_ids,
        duplicate_site_policy=DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
    )

    pdt.assert_frame_equal(first.phospho, second.phospho)
    pdt.assert_frame_equal(first.site_metadata, second.site_metadata)
    pdt.assert_frame_equal(
        first.duplicate_site_resolution, second.duplicate_site_resolution
    )


def test_profile_scoring_explicitly_marks_zero_variance_cases_as_nan() -> None:
    phospho = pd.DataFrame(
        {"a": [1.0, 2.0], "b": [2.0, 2.0], "c": [3.0, 2.0]},
        index=pd.Index(["SITE_VAR", "SITE_CONST"], name="site_id"),
    )
    profiles = pd.DataFrame(
        {"a": [1.0, 4.0], "b": [2.0, 4.0], "c": [3.0, 4.0]},
        index=pd.Index(["KINASE_VAR", "KINASE_CONST"], name="kinase"),
    )

    scores = score_profile_correlations(phospho=phospho, profile_matrix=profiles)

    assert scores.at["SITE_VAR", "KINASE_VAR"] == pytest.approx(1.0)
    assert pd.isna(scores.at["SITE_VAR", "KINASE_CONST"])
    assert pd.isna(scores.at["SITE_CONST", "KINASE_VAR"])
    assert pd.isna(scores.at["SITE_CONST", "KINASE_CONST"])


def test_profile_self_inclusion_behavior_is_intentional() -> None:
    phospho = pd.DataFrame(
        {"a": [1.0, 10.0], "b": [2.0, 10.0], "c": [3.0, 10.0]},
        index=pd.Index(["SITE1", "SITE2"], name="site_id"),
    )
    substrate_map = pd.DataFrame(
        {
            "kinase": ["K1"],
            "substrate_site": ["SITE1"],
        }
    )

    profiles = build_kinase_profiles(
        phospho=phospho,
        kinase_substrate_map=substrate_map,
        min_substrates=1,
        allow_single_substrate_profiles=True,
    )
    scores = score_profile_correlations(
        phospho=phospho.loc[["SITE1"]],
        profile_matrix=profiles.profile_matrix,
    )

    assert profiles.quantified_substrates["K1"] == ["SITE1"]
    assert scores.at["SITE1", "K1"] == pytest.approx(1.0)


def test_substrate_mean_activity_matches_hand_calculated_toy_example() -> None:
    pred_mat = pd.DataFrame(
        {"K1": [0.9, 0.7, 0.1]},
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;"], name="site_id"),
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_a": [10.0, 4.0, 1.0],
            "sample_b": [0.0, 2.0, 3.0],
        },
        index=pred_mat.index.copy(),
    )
    inputs = KinaseActivityInputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.5,
        min_substrates=2,
        top_n_substrates=3,
        overlap_summary=PredMatOverlapSummary(
            overlap_count=3,
            pred_mat_rows=3,
            phospho_rows=3,
        ),
    )

    result = compute_activity_from_inputs(inputs)

    assert result.thresholded_substrate_counts.to_dict() == {"K1": 2}
    assert result.thresholded_substrate_mean_activity.at[
        "K1", "sample_a"
    ] == pytest.approx(7.0)
    assert result.thresholded_substrate_mean_activity.at[
        "K1", "sample_b"
    ] == pytest.approx(1.0)


def test_activity_ranking_tie_behavior_is_deterministic() -> None:
    score_matrix = pd.DataFrame(
        {"K2": [0.5, 0.5], "K1": [0.5, 0.5]},
        index=pd.Index(["S1", "S2"], name="site_id"),
    )
    candidates = {"K2": ["S1", "S2"], "K1": ["S1", "S2"]}

    first = rank_kinases_for_prediction(
        prediction_score_matrix=score_matrix,
        candidate_substrates=candidates,
    )
    second = rank_kinases_for_prediction(
        prediction_score_matrix=score_matrix,
        candidate_substrates=candidates,
    )

    pdt.assert_series_equal(first, second)
    assert list(first.index) == ["K2", "K1"]


def test_motif_profile_fusion_has_predictable_missing_motif_behavior() -> None:
    profile_scores = pd.DataFrame(
        {"K1": [0.8, 0.6], "K2": [0.3, 0.7]},
        index=pd.Index(["S1", "S2"], name="site_id"),
    )
    motif_scores = pd.DataFrame(
        {"K1": [float("nan"), 0.2]},
        index=profile_scores.index.copy(),
    )
    motif_sizes = pd.Series({"K1": 10.0})
    profile_sizes = pd.Series({"K1": 10.0})

    fused, weights = fuse_profile_and_motif_scores_by_rank_weight(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
    )

    assert fused.at["S1", "K1"] == pytest.approx(0.8)
    assert fused.at["S2", "K1"] == pytest.approx(0.4)
    assert fused.loc[:, "K2"].tolist() == pytest.approx([0.3, 0.7])
    assert weights is not None
    assert weights.at["K2", "motif_weight"] == pytest.approx(0.0)
    assert weights.at["K2", "profile_weight"] == pytest.approx(1.0)


def test_signalome_module_assignments_are_row_order_invariant_by_site_identity() -> (
    None
):
    prediction = pd.DataFrame(
        {
            "K1": [0.9, 0.2, 0.8],
            "K2": [0.1, 0.9, 0.8],
        },
        index=pd.Index(["P1;S1;", "P2;S2;", "P3;S3;"], name="site_id"),
    )
    proteins = pd.Series(
        ["P1", "P2", "P3"],
        index=prediction.index.copy(),
        name="protein_id",
        dtype=str,
    )

    reordered_prediction = prediction.loc[["P3;S3;", "P1;S1;", "P2;S2;"], :]
    reordered_proteins = proteins.loc[reordered_prediction.index]

    base = build_module_assignments(
        prediction_matrix=prediction,
        site_to_protein=proteins,
    )
    reordered = build_module_assignments(
        prediction_matrix=reordered_prediction,
        site_to_protein=reordered_proteins,
    )

    pdt.assert_frame_equal(base.sort_index(), reordered.sort_index())


def test_signalome_missing_value_preconditioning_is_recorded_in_diagnostics_and_provenance() -> (
    None
):
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _analysis_ready_dataset(site_ids)
    prediction_matrix = pd.DataFrame(
        {"K1": [0.9, 0.2, 0.7], "K2": [0.1, 0.8, 0.3]},
        index=pd.Index(site_ids, name="site_id"),
    )
    score_matrix = pd.DataFrame(
        {"K1": [0.95, 0.05, float("nan")], "K2": [0.05, 0.95, float("nan")]},
        index=pd.Index(site_ids, name="site_id"),
    )
    kinase_result = _kinase_result(
        dataset=dataset,
        prediction_matrix=prediction_matrix,
        score_matrix=score_matrix,
    )

    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(
                module_count=2,
                score_preconditioning_policy=(
                    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
                ),
            ),
        )
    )

    diagnostics = result.provenance.workflow_parameters[
        "score_preconditioning_diagnostics"
    ]
    assert diagnostics["input_row_count"] == 3
    assert diagnostics["dropped_all_missing_row_count"] == 1
    assert diagnostics["retained_row_count"] == 2
    assert diagnostics["policy"] == "allow_and_report"

    policy_ids = {policy.id for policy in result.provenance.scientific_policies}
    assert "signalome_score_preconditioning_v1" in {
        policy_id.value for policy_id in policy_ids
    }


def test_signalome_backend_selection_metadata_is_recorded_clearly() -> None:
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;", "P4;S4;"]
    dataset = _analysis_ready_dataset(site_ids)
    prediction_matrix = pd.DataFrame(
        {
            "K1": [0.9, 0.2, 0.8, 0.7],
            "K2": [0.1, 0.8, 0.2, 0.3],
            "K3": [0.3, 0.4, 0.7, 0.9],
        },
        index=pd.Index(site_ids, name="site_id"),
    )
    score_matrix = pd.DataFrame(
        {
            "K1": [0.95, 0.1, 0.9, 0.7],
            "K2": [0.05, 0.9, 0.2, 0.3],
            "K3": [0.2, 0.4, 0.75, 0.85],
        },
        index=pd.Index(site_ids, name="site_id"),
    )
    kinase_result = _kinase_result(
        dataset=dataset,
        prediction_matrix=prediction_matrix,
        score_matrix=score_matrix,
    )

    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(
                module_count=None,
                clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
                candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
                max_full_candidate_scoring_sites=2,
            ),
        )
    )

    scale_guard = result.provenance.workflow_parameters["scale_guard"]
    assert scale_guard["clustering_engine"] == SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    assert scale_guard["candidate_scoring_mode"] == "sampled"
    assert scale_guard["candidate_scoring_is_approximate"] is True
    assert scale_guard["tree_generation_is_approximate"] is False


def test_signalome_module_candidate_scoring_is_deterministic_for_fixed_inputs() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 0.0, 0.4],
            [0.9, 0.1, 0.5],
            [0.1, 1.0, 0.4],
            [0.2, 0.9, 0.6],
        ],
        dtype=float,
    )
    clustering_values = prepare_scoring_values_for_clustering(scoring_values)
    profile_degeneracy = summarize_profile_degeneracy(scoring_values)

    first = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=scoring_values,
        candidate_range=range(2, 4),
        profile_degeneracy=profile_degeneracy,
        n_sites=scoring_values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=100,
        max_full_candidate_scoring_sites=2,
    )
    second = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=scoring_values,
        candidate_range=range(2, 4),
        profile_degeneracy=profile_degeneracy,
        n_sites=scoring_values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=100,
        max_full_candidate_scoring_sites=2,
    )

    assert first.candidate_scores == second.candidate_scores
    assert first.candidate_scoring_sampling == second.candidate_scoring_sampling
