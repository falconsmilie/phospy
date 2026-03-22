from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.prediction import (
    KinasePredictionResult,
    KinasePredictor,
    _build_coverage_negative_batches,
    build_candidate_substrate_list,
)
from phospy.scoring import KinaseScoringResult


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


def test_build_candidate_substrate_list_breaks_ties_deterministically() -> None:
    tied = pd.DataFrame(
        {
            "KINASE_A": [0.9, 0.9, 0.9, 0.8],
        },
        index=["SITE_3", "SITE_1", "SITE_2", "SITE_4"],
    )

    substrate_list = build_candidate_substrate_list(
        tied,
        top=3,
        score_threshold=0.1,
        inclusion=1,
    )

    assert substrate_list["KINASE_A"] == ["SITE_1", "SITE_2", "SITE_3"]


def test_coverage_negative_batches_expand_pool_exposure() -> None:
    rng = np.random.default_rng(7)
    batches = _build_coverage_negative_batches(
        negative_index=np.array([f"NEG_{idx}" for idx in range(1, 11)], dtype=object),
        batch_size=4,
        ensemble_size=3,
        rng=rng,
    )

    assert len(batches) == 3
    assert all(len(batch) == 4 for batch in batches)
    distinct_seen = {str(site) for batch in batches for site in batch.tolist()}
    assert len(distinct_seen) >= 10


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


def test_predict_supports_multiple_negative_sampling_strategies() -> None:
    predictor = KinasePredictor()

    random_result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=4,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=11,
        negative_sampling_strategy="random",
    )
    coverage_result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=4,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=11,
        negative_sampling_strategy="coverage",
    )
    hybrid_result = predictor.predict(
        combined_scores=make_combined_scores(),
        ensemble_size=4,
        top=4,
        score_threshold=0.85,
        inclusion=3,
        n_iterations=2,
        random_state=11,
        negative_sampling_strategy="hybrid",
    )

    assert list(random_result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]
    assert list(coverage_result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]
    assert list(hybrid_result.pred_matrix.columns) == ["KINASE_A", "KINASE_B"]


def test_predict_rejects_unknown_negative_sampling_strategy() -> None:
    predictor = KinasePredictor()

    with pytest.raises(ValueError, match="negative_sampling_strategy"):
        predictor.predict(
            combined_scores=make_combined_scores(),
            negative_sampling_strategy="nope",
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
