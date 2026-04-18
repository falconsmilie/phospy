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

# Fixture provenance: PhosR-aligned reference exported in
# tests_legacy/fixtures/r_reference_l6/native_profile_scores.csv.
pytestmark = pytest.mark.parity


def test_scoring_outputs_match_selected_reference_profile_values() -> None:
    dataset = build_rat_l6_dataset(n_sites=None)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=1),
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
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=3, ensemble_size=200),
            activity_config=None,
        )
    )
    expected = load_expected_profile_scores()
    substrate_map = result.references.kinase_substrate_map
    for kinase in ("AKT1", "MAPK1"):
        candidates = [
            site_id
            for site_id in substrate_map.loc[
                substrate_map.loc[:, "kinase"] == kinase, "substrate_site"
            ].astype(str)
            if site_id in expected.index
            and site_id in result.scoring_result.profile_scores.index
        ]
        expected_top = expected.loc[candidates, kinase].astype(float).idxmax()
        observed_top = result.prediction_result.substrate_list.loc[
            (result.prediction_result.substrate_list.loc[:, "kinase"] == kinase)
            & (result.prediction_result.substrate_list.loc[:, "rank"] == 1),
            "substrate_site",
        ].iloc[0]
        assert observed_top == expected_top
