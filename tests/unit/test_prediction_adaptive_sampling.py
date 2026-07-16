from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.api.configs import KinasePredictionConfig
from phospy.errors import WorkflowStageError, WorkflowValidationError
from phospy.science.prediction import sampling_core
from phospy.science.prediction.execution import run_adaptive_ensemble_prediction
from phospy.science.prediction.policies import (
    resolve_prediction_sampling_policy,
)
from phospy.science.prediction.sampling_core import run_adaptive_sampling_ensemble
from phospy.science.prediction.sampling_runtime import (
    PredictionSamplingRandomSource,
    transform_resampling_probabilities,
)


def test_kinase_prediction_config_rejects_adaptive_mode_without_seed() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match=(
            "prediction_config.random_state must be provided when "
            "prediction_config.mode='adaptive_ensemble'"
        ),
    ):
        KinasePredictionConfig(mode="adaptive_ensemble", random_state=None)


def test_kinase_prediction_config_allows_adaptive_mode_with_seed() -> None:
    config = KinasePredictionConfig(mode="adaptive_ensemble", random_state=1)
    assert config.random_state == 1


def test_resolve_prediction_sampling_policy_maps_public_modes() -> None:
    stable = resolve_prediction_sampling_policy("stable")
    r_parity = resolve_prediction_sampling_policy("r_parity")

    assert stable.seed_strategy == "stable_by_kinase"
    assert stable.resampling_weight_mode == "default"
    assert stable.final_score_mode == "mean_probability"
    assert r_parity.seed_strategy == "global_parity"
    assert r_parity.resampling_weight_mode == "r_parity"
    assert r_parity.final_score_mode == "decision_sigmoid"


def test_prediction_sampling_random_source_stable_policy_is_order_invariant() -> None:
    policy = resolve_prediction_sampling_policy("stable")
    source_a = PredictionSamplingRandomSource(policy=policy, random_state=17)
    source_b = PredictionSamplingRandomSource(policy=policy, random_state=17)

    a_first = source_a.generators_for_kinase(kinase="KINASE_A")
    _ = source_a.generators_for_kinase(kinase="KINASE_B")
    _ = source_b.generators_for_kinase(kinase="KINASE_B")
    b_second = source_b.generators_for_kinase(kinase="KINASE_A")

    assert int(a_first[0].integers(0, 1000)) == int(b_second[0].integers(0, 1000))
    assert int(a_first[1].integers(0, 1000)) == int(b_second[1].integers(0, 1000))


def test_prediction_sampling_random_source_global_parity_tracks_call_order() -> None:
    policy = resolve_prediction_sampling_policy("r_parity")
    source_a = PredictionSamplingRandomSource(policy=policy, random_state=17)
    source_b = PredictionSamplingRandomSource(policy=policy, random_state=17)

    a_first = source_a.generators_for_kinase(kinase="KINASE_A")
    _ = source_a.generators_for_kinase(kinase="KINASE_B")
    _ = source_b.generators_for_kinase(kinase="KINASE_B")
    b_second = source_b.generators_for_kinase(kinase="KINASE_A")

    assert int(a_first[0].integers(0, 1000)) != int(b_second[0].integers(0, 1000))


def test_transform_resampling_probabilities_flattens_stable_policy() -> None:
    values = np.asarray([0.9, 0.1], dtype=float)
    transformed = transform_resampling_probabilities(values, adaptive_policy="stable")
    normalized_input = values / values.sum()
    normalized_transformed = transformed / transformed.sum()

    assert normalized_transformed[0] < normalized_input[0]
    assert normalized_transformed[1] > normalized_input[1]


def test_run_adaptive_sampling_ensemble_separates_positive_and_negative_sites() -> None:
    feature_values = np.asarray(
        [
            [0.95, 0.1],
            [0.95, 0.2],
            [0.88, 0.2],
            [0.2, 0.88],
            [0.2, 0.95],
            [0.1, 0.95],
        ],
        dtype=float,
    )
    train_values = np.concatenate(
        [feature_values[:3, :], feature_values[3:, :]],
        axis=0,
    )
    labels = np.asarray([1, 1, 1, 2, 2, 2], dtype=int)
    scores = run_adaptive_sampling_ensemble(
        train_values=train_values,
        train_labels=labels,
        test_values=feature_values,
        kernel="rbf",
        n_iterations=2,
        resampling_rng=np.random.default_rng(13),
        sampling_policy=resolve_prediction_sampling_policy("stable"),
    )

    assert scores.shape == (6,)
    assert ((scores >= 0.0) & (scores <= 1.0)).all()
    assert float(scores[:3].mean()) > float(scores[3:].mean())


def test_run_adaptive_sampling_ensemble_skips_decision_vector_for_mean_probability(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _ProbabilityModel:
        classes_ = np.asarray([1, 2], dtype=int)

        @staticmethod
        def predict_proba(values: np.ndarray) -> np.ndarray:
            return np.column_stack(
                [
                    np.linspace(0.2, 0.8, values.shape[0]),
                    np.linspace(0.8, 0.2, values.shape[0]),
                ]
            )

    monkeypatch.setattr(
        sampling_core,
        "_run_sampling_iterations",
        lambda **_kwargs: _ProbabilityModel(),
    )

    def fail_if_called(**_kwargs) -> np.ndarray:
        raise AssertionError("decision vector should not be computed")

    monkeypatch.setattr(
        sampling_core,
        "aligned_binary_decision_vector",
        fail_if_called,
    )

    scores = run_adaptive_sampling_ensemble(
        train_values=np.ones((2, 2), dtype=float),
        train_labels=np.asarray([1, 2], dtype=int),
        test_values=np.ones((3, 2), dtype=float),
        kernel="rbf",
        n_iterations=1,
        resampling_rng=np.random.default_rng(13),
        sampling_policy=resolve_prediction_sampling_policy("stable"),
    )

    assert scores.tolist() == pytest.approx([0.2, 0.5, 0.8])


def test_run_adaptive_sampling_ensemble_uses_decision_vector_for_r_parity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _ProbabilityModel:
        classes_ = np.asarray([1, 2], dtype=int)

        @staticmethod
        def predict_proba(values: np.ndarray) -> np.ndarray:
            return np.column_stack(
                [
                    np.linspace(0.2, 0.8, values.shape[0]),
                    np.linspace(0.8, 0.2, values.shape[0]),
                ]
            )

    monkeypatch.setattr(
        sampling_core,
        "_run_sampling_iterations",
        lambda **_kwargs: _ProbabilityModel(),
    )
    monkeypatch.setattr(
        sampling_core,
        "aligned_binary_decision_vector",
        lambda **_kwargs: np.asarray([-1.0, 0.0, 1.0], dtype=float),
    )

    scores = run_adaptive_sampling_ensemble(
        train_values=np.ones((2, 2), dtype=float),
        train_labels=np.asarray([1, 2], dtype=int),
        test_values=np.ones((3, 2), dtype=float),
        kernel="rbf",
        n_iterations=1,
        resampling_rng=np.random.default_rng(13),
        sampling_policy=resolve_prediction_sampling_policy("r_parity"),
    )

    assert scores.tolist() == pytest.approx(
        [1.0 / (1.0 + np.exp(1.0)), 0.5, 1.0 / (1.0 + np.exp(-1.0))]
    )


def test_run_adaptive_ensemble_prediction_averages_per_ensemble_scores(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    score_matrix = pd.DataFrame(
        {"K1": [0.9, 0.5, 0.2]},
        index=["S1", "S2", "S3"],
    )
    call_count = {"count": 0}

    def fake_sampling(**kwargs) -> np.ndarray:
        del kwargs
        call_count["count"] += 1
        if call_count["count"] == 1:
            return np.asarray([0.2, 0.4, 0.6], dtype=float)
        return np.asarray([0.4, 0.6, 0.8], dtype=float)

    monkeypatch.setattr(
        "phospy.science.prediction.execution.run_adaptive_sampling_ensemble",
        fake_sampling,
    )
    observed = run_adaptive_ensemble_prediction(
        prediction_score_matrix=score_matrix,
        candidate_substrates={"K1": ["S1"]},
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            mode="adaptive_ensemble",
            n_iterations=2,
            random_state=5,
        ),
        random_state=5,
    )

    assert call_count["count"] == 2
    assert observed.loc[:, "K1"].tolist() == pytest.approx([0.3, 0.5, 0.7])


def test_run_adaptive_ensemble_prediction_requires_negative_pool() -> None:
    score_matrix = pd.DataFrame(
        {"K1": [0.9, 0.8]},
        index=["S1", "S2"],
    )
    with pytest.raises(WorkflowStageError, match="prediction.adaptive_negative_pool"):
        run_adaptive_ensemble_prediction(
            prediction_score_matrix=score_matrix,
            candidate_substrates={"K1": ["S1", "S2"]},
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=1,
                adaptive_ensemble_runs=1,
                mode="adaptive_ensemble",
                n_iterations=1,
                random_state=0,
            ),
            random_state=0,
        )


def test_adaptive_prediction_is_deterministic_for_same_seed() -> None:
    score_matrix = pd.DataFrame(
        {
            "K1": [0.95, 0.91, 0.88, 0.22, 0.18, 0.11],
            "K2": [0.12, 0.2, 0.28, 0.8, 0.86, 0.93],
        },
        index=["S1", "S2", "S3", "S4", "S5", "S6"],
    )
    candidate_substrates = {
        "K1": ["S1", "S2", "S3"],
        "K2": ["S4", "S5", "S6"],
    }
    prediction_config = KinasePredictionConfig(
        top_k=3,
        deterministic_max_selected_kinases=2,
        adaptive_ensemble_runs=4,
        mode="adaptive_ensemble",
        n_iterations=3,
        random_state=17,
    )

    first = run_adaptive_ensemble_prediction(
        prediction_score_matrix=score_matrix,
        candidate_substrates=candidate_substrates,
        prediction_config=prediction_config,
        random_state=17,
    )
    second = run_adaptive_ensemble_prediction(
        prediction_score_matrix=score_matrix,
        candidate_substrates=candidate_substrates,
        prediction_config=prediction_config,
        random_state=17,
    )

    pd.testing.assert_frame_equal(first, second)


def test_adaptive_prediction_outputs_can_differ_for_different_seeds() -> None:
    score_matrix = pd.DataFrame(
        {
            "K1": [0.95, 0.91, 0.88, 0.22, 0.18, 0.11],
            "K2": [0.12, 0.2, 0.28, 0.8, 0.86, 0.93],
        },
        index=["S1", "S2", "S3", "S4", "S5", "S6"],
    )
    candidate_substrates = {
        "K1": ["S1", "S2", "S3"],
        "K2": ["S4", "S5", "S6"],
    }

    first = run_adaptive_ensemble_prediction(
        prediction_score_matrix=score_matrix,
        candidate_substrates=candidate_substrates,
        prediction_config=KinasePredictionConfig(
            top_k=3,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=4,
            mode="adaptive_ensemble",
            n_iterations=3,
            random_state=17,
        ),
        random_state=17,
    )
    second = run_adaptive_ensemble_prediction(
        prediction_score_matrix=score_matrix,
        candidate_substrates=candidate_substrates,
        prediction_config=KinasePredictionConfig(
            top_k=3,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=4,
            mode="adaptive_ensemble",
            n_iterations=3,
            random_state=23,
        ),
        random_state=23,
    )

    assert not first.equals(second)
