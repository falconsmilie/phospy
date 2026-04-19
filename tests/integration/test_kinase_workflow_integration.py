from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.workflows.kinase.executor as kinase_executor
from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    ReferencePreset,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset

pytestmark = pytest.mark.integration


def test_kinase_workflow_runs_without_dataset_site_sequence_column() -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    dataset_without_sequence = AnalysisReadyPhosphoDataset(
        phospho=dataset.phospho,
        site_metadata=dataset.site_metadata.drop(columns=["site_sequence"]),
        transformation_state=dataset.transformation_state,
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
    assert result.scoring_result.motif_scores is not None
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
    assert result.scoring_result.motif_scores is not None
    assert result.scoring_result.combined_scores is not None
    assert result.scoring_result.weights is not None
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
    assert not result.activity_result.ksea_scores.empty
    assert not result.activity_result.ksea_counts.empty
    assert not result.activity_result.target_counts.empty
    assert {"site_id", "kinase", "score"} <= set(
        result.activity_result.target_table.columns
    )
    assert not hasattr(result, "profile_scores")
    assert not hasattr(result, "combined_scores")
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

    def _force_profile_lane(*, profile_scores, combined_scores):
        _ = combined_scores
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

    combined_scores = combined_lane.scoring_result.combined_scores
    assert combined_scores is not None
    profile_scores = combined_lane.scoring_result.profile_scores
    shared_kinases = pd.Index(combined_pred.columns).intersection(profile_pred.columns)
    assert not shared_kinases.empty

    matched_a_differing_kinase = False
    for kinase in shared_kinases.astype(str):
        if combined_pred.loc[:, kinase].dropna().empty:
            continue
        if profile_pred.loc[:, kinase].dropna().empty:
            continue
        combined_top = combined_scores.loc[:, kinase].astype(float).idxmax()
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
