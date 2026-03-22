from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    KinasePredictionResult,
    KinasePredictor,
    KinaseScoringResult,
    build_candidate_substrate_list,
)


def make_combined_scores() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.93, 0.91, 0.89, 0.20, 0.18, 0.16, 0.14],
            "KINASE_B": [0.10, 0.12, 0.14, 0.16, 0.96, 0.94, 0.92, 0.90],
        },
        index=[f"SITE_{i}" for i in range(1, 9)],
    )


def test_build_candidate_substrate_list_applies_selection_rules() -> None:
    substrate_list = build_candidate_substrate_list(
        make_combined_scores(),
        top=4,
        score_threshold=0.9,
        inclusion=2,
    )

    assert substrate_list == {
        "KINASE_A": ["SITE_1", "SITE_2", "SITE_3"],
        "KINASE_B": ["SITE_5", "SITE_6", "SITE_7"],
    }


def test_kinase_predictor_returns_probability_matrix() -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=3,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
    )

    assert isinstance(result, KinasePredictionResult)
    assert list(result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]
    assert ((result.pred_matrix >= 0.0) & (result.pred_matrix <= 1.0)).all().all()
    assert (
        result.pred_matrix.loc[
            ["SITE_1", "SITE_2", "SITE_3", "SITE_4"], "KINASE_A"
        ].mean()
        > result.pred_matrix.loc[
            ["SITE_5", "SITE_6", "SITE_7", "SITE_8"], "KINASE_A"
        ].mean()
    )
    assert (
        result.pred_matrix.loc[
            ["SITE_5", "SITE_6", "SITE_7", "SITE_8"], "KINASE_B"
        ].mean()
        > result.pred_matrix.loc[
            ["SITE_1", "SITE_2", "SITE_3", "SITE_4"], "KINASE_B"
        ].mean()
    )


def test_predict_from_scoring_result_uses_combined_scores() -> None:
    predictor = KinasePredictor()
    scoring_result = KinaseScoringResult(
        profile_scores=make_combined_scores(),
        combined_scores=make_combined_scores(),
    )

    result = predictor.predict_from_scoring_result(
        scoring_result,
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=9,
    )

    assert set(result.substrate_list) == {"KINASE_A", "KINASE_B"}


def test_predict_from_scoring_result_requires_explicit_profile_fallback() -> None:
    predictor = KinasePredictor()
    scoring_result = KinaseScoringResult(profile_scores=make_combined_scores())

    with pytest.raises(ValueError, match="combined_scores"):
        predictor.predict_from_scoring_result(scoring_result)

    result = predictor.predict_from_scoring_result(
        scoring_result,
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=9,
        allow_profile_only_fallback=True,
    )

    assert list(result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]


def test_predict_returns_empty_matrix_when_no_kinases_pass_inclusion() -> None:
    predictor = KinasePredictor()
    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=2,
        score_threshold=0.99,
        inclusion=2,
        n_iterations=2,
        random_state=3,
    )

    assert result.substrate_list == {}
    assert result.pred_matrix.empty
    assert list(result.pred_matrix.index) == list(make_combined_scores().index)


def test_predict_can_capture_debug_trace_for_selected_kinase() -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=2,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
        debug_kinases=["KINASE_A"],
        debug_top_n=3,
    )

    assert result.debug_traces is not None
    assert set(result.debug_traces) == {"KINASE_A"}
    trace = result.debug_traces["KINASE_A"]
    assert trace.candidate_substrates == ["SITE_1", "SITE_2", "SITE_3", "SITE_4"]
    assert len(trace.negative_pool_sites) == 4
    assert len(trace.ensemble_traces) == 2
    ensemble_trace = trace.ensemble_traces[0]
    assert len(ensemble_trace.initial_negative_sites) == 4
    assert len(ensemble_trace.iterations) == 2
    assert list(ensemble_trace.final_prediction_probabilities.columns) == ["1", "2"]
    assert len(ensemble_trace.final_top_sites) == 3
    iteration_trace = ensemble_trace.iterations[0]
    assert set(iteration_trace.labels.unique()) == {1, 2}
    assert list(iteration_trace.probabilities.columns) == ["1", "2"]
    assert len(iteration_trace.sampled_positive_sites) == 4
    assert len(iteration_trace.sampled_negative_sites) == 4
    assert iteration_trace.positive_weights is not None
    assert iteration_trace.negative_weights is not None


def test_predict_can_capture_debug_trace_for_all_kinases() -> None:
    predictor = KinasePredictor()

    result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=5,
        capture_debug_trace=True,
    )

    assert result.debug_traces is not None
    assert set(result.debug_traces) == {"KINASE_A", "KINASE_B"}


def test_predict_from_scoring_result_passes_debug_options_through() -> None:
    predictor = KinasePredictor()
    scoring_result = KinaseScoringResult(
        profile_scores=make_combined_scores(),
        combined_scores=make_combined_scores(),
    )

    result = predictor.predict_from_scoring_result(
        scoring_result,
        ensemble_size=1,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=9,
        capture_debug_trace=True,
        debug_kinases=["KINASE_B"],
        debug_top_n=2,
    )

    assert result.debug_traces is not None
    assert set(result.debug_traces) == {"KINASE_B"}
    assert len(result.debug_traces["KINASE_B"].ensemble_traces[0].final_top_sites) == 2


def test_build_candidate_substrate_list_preserves_input_order_for_ties() -> None:
    scores = pd.DataFrame(
        {"KINASE_A": [0.95, 0.95, 0.90, 0.10]},
        index=["SITE_A", "SITE_B", "SITE_C", "SITE_D"],
    )

    substrate_list = build_candidate_substrate_list(
        scores,
        top=3,
        score_threshold=0.5,
        inclusion=2,
    )

    assert substrate_list == {"KINASE_A": ["SITE_A", "SITE_B", "SITE_C"]}
