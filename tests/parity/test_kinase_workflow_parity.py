from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
)
from tests.support.parity_reporting import (
    format_bool,
    format_shape,
    record_parity_metrics,
)
from tests.support.rewrite_fixture_data import (
    build_rat_l6_dataset,
    load_expected_profile_scores,
)

# Fixture provenance:
# tests/fixtures/rewrite_parity/r_reference_l6/native_profile_scores.csv
# (see PROVENANCE.md in the same directory).
pytestmark = pytest.mark.parity


def test_scoring_outputs_match_selected_reference_profile_values(
    request: pytest.FixtureRequest,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=None)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
            ),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    expected = load_expected_profile_scores()
    points = [
        ("AAK1;S677;", "AKT1"),
        ("ABCC4;S604;", "MAPK1"),
        ("ABI2;S165;", "PRKAA1"),
    ]
    point_abs_deltas: list[float] = []
    for site_id, kinase in points:
        observed_value = float(result.scoring_result.profile_scores.at[site_id, kinase])
        expected_value = float(expected.at[site_id, kinase])
        assert observed_value == pytest.approx(
            expected_value,
            rel=1e-6,
            abs=1e-8,
        )
        point_abs_deltas.append(abs(observed_value - expected_value))

    record_parity_metrics(
        request.config,
        family="kinase_workflow",
        metrics=[
            ("dataset site count", dataset.phospho.shape[0]),
            (
                "profile score table shape",
                format_shape(*result.scoring_result.profile_scores.shape),
            ),
            (
                "kinases scored",
                int(result.scoring_result.profile_scores.shape[1]),
            ),
            ("selected reference-point count", len(points)),
            (
                "selected reference-point mean abs diff",
                float(pd.Series(point_abs_deltas).mean()),
            ),
            (
                "selected reference-point max abs diff",
                float(pd.Series(point_abs_deltas).max()),
            ),
        ],
    )


def test_prediction_top_sites_align_with_reference_ranking_subset(
    request: pytest.FixtureRequest,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=None)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=3,
                deterministic_max_selected_kinases=200,
                adaptive_ensemble_runs=200,
            ),
            activity_config=None,
        )
    )
    rank_weighted_fusion_scores = result.scoring_result.rank_weighted_fusion_scores
    assert rank_weighted_fusion_scores is not None
    profile_scores = result.scoring_result.profile_scores
    pred_mat = result.prediction_result.pred_mat

    differing_kinases = 0
    for kinase in pred_mat.columns.astype(str):
        observed_top = pred_mat.loc[:, kinase].astype(float).dropna().idxmax()
        rank_weighted_fusion_top = (
            rank_weighted_fusion_scores.loc[:, kinase].astype(float).idxmax()
        )
        assert observed_top == rank_weighted_fusion_top
        profile_top = profile_scores.loc[:, kinase].astype(float).idxmax()
        if profile_top != rank_weighted_fusion_top:
            differing_kinases += 1

    assert differing_kinases > 0
    record_parity_metrics(
        request.config,
        family="kinase_workflow",
        metrics=[
            ("prediction matrix shape", format_shape(*pred_mat.shape)),
            ("kinases predicted", int(pred_mat.shape[1])),
            (
                "top-site rank-weighted-fusion-score alignment",
                f"{pred_mat.shape[1]}/{pred_mat.shape[1]}",
            ),
            (
                "rank-weighted-fusion-vs-profile top-site divergences",
                differing_kinases,
            ),
        ],
    )


def test_scoring_outputs_include_motif_and_rank_weighted_fusion_tables(
    request: pytest.FixtureRequest,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
            ),
            prediction_config=KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )

    assert result.scoring_result.motif_scores is not None
    assert result.scoring_result.rank_weighted_fusion_scores is not None
    assert result.scoring_result.score_fusion_weights is not None
    assert not result.scoring_result.motif_scores.empty
    assert not result.scoring_result.rank_weighted_fusion_scores.empty
    assert set(result.scoring_result.score_fusion_weights.columns) == {
        "motif_weight",
        "profile_weight",
        "motif_rank_weight",
        "profile_rank_weight",
    }
    record_parity_metrics(
        request.config,
        family="kinase_workflow",
        metrics=[
            (
                "diagnostic motif table present",
                format_bool(result.scoring_result.motif_scores is not None),
            ),
            (
                "diagnostic rank-weighted-fusion table present",
                format_bool(
                    result.scoring_result.rank_weighted_fusion_scores is not None
                ),
            ),
            (
                "diagnostic weight table present",
                format_bool(result.scoring_result.score_fusion_weights is not None),
            ),
            (
                "diagnostic motif score shape",
                format_shape(*result.scoring_result.motif_scores.shape),
            ),
            (
                "diagnostic rank-weighted-fusion score shape",
                format_shape(*result.scoring_result.rank_weighted_fusion_scores.shape),
            ),
        ],
    )


def test_profile_missing_value_policy_changes_downstream_lane_for_mixed_missing_input() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0, 2.0, 4.0],
            "sample_b": [2.0, 4.0, 1.0, 3.0],
            "sample_c": [3.0, 5.0, 0.0, 2.0],
        },
        index=pd.Index(
            ["GENEA;S1;", "GENEA;S2;", "GENEB;S3;", "GENEB;S4;"],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["GENEA", "GENEA", "GENEB", "GENEB"],
            "site": ["S1", "S2", "S3", "S4"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["S1", "S2", "S3", "S4"]
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    dataset._phospho.loc["GENEA;S2;", "sample_c"] = float("nan")
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K_MISSING", "K_MISSING", "K_STABLE", "K_STABLE"],
                "substrate_site": [
                    "GENEA;S1;",
                    "GENEA;S2;",
                    "GENEB;S3;",
                    "GENEB;S4;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": site_metadata.loc[:, "site_sequence"].values},
            index=phospho.index.copy(),
        ),
    )

    strict = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                profile_missing_value_strategy="strict",
            ),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )
    median_skipna = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                profile_missing_value_strategy="median_skipna",
            ),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    assert strict.scoring_result.profile_scores.loc[:, "K_MISSING"].isna().all()
    assert median_skipna.scoring_result.profile_scores.loc[:, "K_MISSING"].notna().any()
    assert strict.scoring_result.rank_weighted_fusion_scores is not None
    assert median_skipna.scoring_result.rank_weighted_fusion_scores is not None
    assert (
        strict.scoring_result.rank_weighted_fusion_scores.loc[:, "K_MISSING"]
        .isna()
        .all()
    )
    assert (
        median_skipna.scoring_result.rank_weighted_fusion_scores.loc[:, "K_MISSING"]
        .notna()
        .any()
    )
    assert "K_MISSING" not in strict.prediction_result.pred_mat.columns
    assert "K_MISSING" in median_skipna.prediction_result.pred_mat.columns
