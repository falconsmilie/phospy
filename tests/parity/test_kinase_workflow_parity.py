from __future__ import annotations

import pytest

from phospy import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    ReferencePreset,
)
from tests.support.rewrite_fixture_data import (
    build_rat_l6_dataset,
    load_expected_profile_scores,
)

# Fixture provenance:
# tests/fixtures/rewrite_parity/r_reference_l6/native_profile_scores.csv
# (see PROVENANCE.md in the same directory).
pytestmark = pytest.mark.parity


def test_scoring_outputs_match_selected_reference_profile_values() -> None:
    dataset = build_rat_l6_dataset(n_sites=None)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=12),
            activity_config=None,
        )
    )
    expected = load_expected_profile_scores()
    points = [
        ("AAK1;S677;", "AKT1"),
        ("ABCC4;S604;", "MAPK1"),
        ("ABI2;S165;", "PRKAA1"),
    ]
    for site_id, kinase in points:
        assert result.scoring_result.profile_scores.at[
            site_id, kinase
        ] == pytest.approx(
            expected.at[site_id, kinase],
            rel=1e-6,
            abs=1e-8,
        )


def test_prediction_top_sites_align_with_reference_ranking_subset() -> None:
    dataset = build_rat_l6_dataset(n_sites=None)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=3, ensemble_size=200),
            activity_config=None,
        )
    )
    combined_scores = result.scoring_result.combined_scores
    assert combined_scores is not None
    profile_scores = result.scoring_result.profile_scores
    pred_mat = result.prediction_result.pred_mat

    differing_kinases = 0
    for kinase in pred_mat.columns.astype(str):
        observed_top = pred_mat.loc[:, kinase].astype(float).dropna().idxmax()
        combined_top = combined_scores.loc[:, kinase].astype(float).idxmax()
        assert observed_top == combined_top
        profile_top = profile_scores.loc[:, kinase].astype(float).idxmax()
        if profile_top != combined_top:
            differing_kinases += 1

    assert differing_kinases > 0


def test_scoring_outputs_include_motif_and_combined_tables() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=12),
            activity_config=None,
        )
    )

    assert result.scoring_result.motif_scores is not None
    assert result.scoring_result.combined_scores is not None
    assert result.scoring_result.weights is not None
    assert not result.scoring_result.motif_scores.empty
    assert not result.scoring_result.combined_scores.empty
    assert set(result.scoring_result.weights.columns) == {
        "motif_weight",
        "profile_weight",
        "motif_rank_weight",
        "profile_rank_weight",
    }
