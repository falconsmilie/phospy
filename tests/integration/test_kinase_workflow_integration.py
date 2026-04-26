from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.workflows.kinase.executor as kinase_executor
from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
)
from phospy.api import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    ReferenceBundle,
    ReferencePreset,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset

pytestmark = pytest.mark.integration


def test_kinase_workflow_runs_without_dataset_site_sequence_column() -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    dataset_without_sequence = AnalysisReadyPhosphoDataset(
        phospho=dataset.phospho,
        site_metadata=dataset.site_metadata.drop(columns=["site_sequence"]),
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=dataset.processing_state,
        sample_metadata=dataset.sample_metadata,
        total=dataset.total,
        organism=dataset.organism,
    )
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset_without_sequence,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=8),
            activity_config=None,
        )
    )
    assert not result.scoring_result.profile_scores.empty
    assert result.scoring_result.motif_scores is None
    assert result.scoring_result.score_fusion_weights is None
    assert not result.prediction_result.pred_mat.empty


def test_kinase_workflow_runs_dataset_to_kinase_path() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=8),
            activity_config=KinaseActivityConfig(
                enabled=True,
                threshold=0.6,
                min_substrates=3,
                top_n_substrates=20,
            ),
        )
    )
    assert result.scoring_result.profile_scores.shape[0] == dataset.phospho.shape[0]
    assert result.scoring_result.profile_scores.shape[1] > 0
    assert result.scoring_result.rank_weighted_fusion_scores is not None
    assert result.scoring_result.motif_scores is None
    assert result.scoring_result.score_fusion_weights is None
    sequence_sites = set(result.references.site_sequences.index.astype(str))
    assert set(result.scoring_result.profile_scores.index.astype(str)).issubset(
        sequence_sites
    )
    assert result.prediction_result.pred_mat.shape[1] <= 8
    pred_values = result.prediction_result.pred_mat.to_numpy(dtype=float)
    finite_values = pred_values[np.isfinite(pred_values)]
    assert (finite_values >= 0.0).all()
    assert result.activity_result is not None
    assert not result.activity_result.weighted_activity.empty
    assert not result.activity_result.thresholded_substrate_mean_activity.empty
    assert not result.activity_result.thresholded_substrate_counts.empty
    assert not result.activity_result.target_counts.empty
    assert {"site_id", "kinase", "score"} <= set(
        result.activity_result.target_table.columns
    )
    assert not hasattr(result, "profile_scores")
    assert not hasattr(result, "rank_weighted_fusion_scores")
    assert not hasattr(result, "weights")
    assert not hasattr(result, "substrate_list")


def test_kinase_workflow_activity_stage_is_optional() -> None:
    dataset = build_rat_l6_dataset(n_sites=180)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=6),
            activity_config=None,
        )
    )
    assert result.activity_result is None


def test_kinase_workflow_default_scoring_floor_supports_realistic_input() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=6),
            activity_config=None,
        )
    )
    assert result.scoring_result.profile_scores.shape[1] > 0


def test_prediction_changes_when_downstream_matrix_switches_profile_vs_combined(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
        activity_config=None,
    )

    combined_lane = KinaseWorkflow().run(request)

    def _force_profile_lane(*, profile_scores, rank_weighted_fusion_scores):
        _ = rank_weighted_fusion_scores
        return profile_scores, "profile_scores"

    monkeypatch.setattr(
        kinase_executor,
        "select_downstream_score_matrix",
        _force_profile_lane,
    )
    profile_lane = KinaseWorkflow().run(request)

    combined_pred = combined_lane.prediction_result.pred_mat
    profile_pred = profile_lane.prediction_result.pred_mat
    assert not combined_pred.equals(profile_pred)

    rank_weighted_fusion_scores = (
        combined_lane.scoring_result.rank_weighted_fusion_scores
    )
    assert rank_weighted_fusion_scores is not None
    profile_scores = combined_lane.scoring_result.profile_scores
    shared_kinases = pd.Index(combined_pred.columns).intersection(profile_pred.columns)
    assert not shared_kinases.empty

    matched_a_differing_kinase = False
    for kinase in shared_kinases.astype(str):
        if combined_pred.loc[:, kinase].dropna().empty:
            continue
        if profile_pred.loc[:, kinase].dropna().empty:
            continue
        combined_top = rank_weighted_fusion_scores.loc[:, kinase].astype(float).idxmax()
        profile_top = profile_scores.loc[:, kinase].astype(float).idxmax()
        if combined_top == profile_top:
            continue
        assert (
            combined_pred.loc[:, kinase].astype(float).dropna().idxmax() == combined_top
        )
        assert (
            profile_pred.loc[:, kinase].astype(float).dropna().idxmax() == profile_top
        )
        matched_a_differing_kinase = True
        break

    assert matched_a_differing_kinase


def test_diagnostic_scoring_tables_are_opt_in_without_changing_supported_lane() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    request_default = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
        activity_config=None,
    )
    request_with_diagnostics = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=True,
        ),
        prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
        activity_config=None,
    )

    default_result = KinaseWorkflow().run(request_default)
    diagnostics_result = KinaseWorkflow().run(request_with_diagnostics)

    assert default_result.scoring_result.motif_scores is None
    assert default_result.scoring_result.score_fusion_weights is None
    assert diagnostics_result.scoring_result.motif_scores is not None
    assert diagnostics_result.scoring_result.score_fusion_weights is not None

    pd.testing.assert_frame_equal(
        default_result.scoring_result.profile_scores,
        diagnostics_result.scoring_result.profile_scores,
        check_dtype=False,
    )
    assert default_result.scoring_result.rank_weighted_fusion_scores is not None
    assert diagnostics_result.scoring_result.rank_weighted_fusion_scores is not None
    pd.testing.assert_frame_equal(
        default_result.scoring_result.rank_weighted_fusion_scores,
        diagnostics_result.scoring_result.rank_weighted_fusion_scores,
        check_dtype=False,
    )
    pd.testing.assert_frame_equal(
        default_result.prediction_result.pred_mat,
        diagnostics_result.prediction_result.pred_mat,
        check_dtype=False,
    )


def test_profile_missing_value_strategy_flows_from_request_to_science(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    captured_strategies: list[str] = []
    original_build_kinase_profiles = kinase_executor.build_kinase_profiles

    def _capture_strategy(**kwargs):
        captured_strategies.append(kwargs["profile_missing_value_strategy"])
        return original_build_kinase_profiles(**kwargs)

    monkeypatch.setattr(
        kinase_executor,
        "build_kinase_profiles",
        _capture_strategy,
    )
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                profile_missing_value_strategy="median_skipna",
            ),
            prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=8),
            activity_config=None,
        )
    )

    assert captured_strategies == ["median_skipna"]
    assert not result.scoring_result.profile_scores.empty


def test_motif_library_build_is_limited_to_profile_eligible_kinases(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=40)
    overlap_sites = dataset.phospho.index.astype(str).tolist()[:2]
    offlane_sites = [f"OFFLANE{i};S{i};" for i in range(1, 9)]
    references = ReferenceBundle(
        organism=dataset.organism,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K_ELIGIBLE", "K_ELIGIBLE"]
                + ["K_OFFLANE" for _ in offlane_sites],
                "substrate_site": overlap_sites + offlane_sites,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31 for _ in overlap_sites + offlane_sites]},
            index=pd.Index(overlap_sites + offlane_sites, name="site_id"),
        ),
    )
    captured_kinases: list[str] = []
    original_build_motif_library = kinase_executor.build_motif_library

    def _capture_eligible_kinases(*, kinase_substrate_map, site_sequences, flank_size):
        nonlocal captured_kinases
        captured_kinases = sorted(
            kinase_substrate_map.loc[:, "kinase"].astype(str).unique().tolist()
        )
        return original_build_motif_library(
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            flank_size=flank_size,
        )

    monkeypatch.setattr(
        kinase_executor,
        "build_motif_library",
        _capture_eligible_kinases,
    )
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
            ),
            prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=1),
            activity_config=None,
        )
    )

    assert captured_kinases == ["K_ELIGIBLE"]
    assert list(result.scoring_result.profile_scores.columns) == ["K_ELIGIBLE"]
    assert result.scoring_result.motif_scores is not None
    assert list(result.scoring_result.motif_scores.columns) == ["K_ELIGIBLE"]
