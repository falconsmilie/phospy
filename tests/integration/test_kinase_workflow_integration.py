from __future__ import annotations

import pytest

from phospy import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    ReferencePreset,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset

pytestmark = pytest.mark.integration


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
    assert result.scoring_result.motif_scores is None
    assert result.scoring_result.combined_scores is None
    assert result.scoring_result.weights is None
    sequence_sites = set(result.references.site_sequences.index.astype(str))
    assert set(result.scoring_result.profile_scores.index.astype(str)).issubset(
        sequence_sites
    )
    assert result.prediction_result.pred_mat.shape[1] <= 8
    assert (result.prediction_result.pred_mat.to_numpy() >= 0.0).all()
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
