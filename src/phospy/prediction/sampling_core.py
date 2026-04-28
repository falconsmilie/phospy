"""Adaptive ensemble sampling core used by kinase prediction workflow."""

from __future__ import annotations

from typing import Protocol, cast

import numpy as np

from phospy.prediction.policies import PredictionSamplingPolicy
from phospy.prediction.sampling_runtime import (
    normalize_probabilities,
    transform_resampling_probabilities,
)
from phospy.prediction.svm import (
    aligned_binary_decision_vector,
    make_svm,
    require_sklearn,
)


class _SamplingModel(Protocol):
    classes_: np.ndarray

    def fit(self, x: np.ndarray, y: np.ndarray) -> None: ...

    def predict_proba(self, x: np.ndarray) -> np.ndarray: ...


def _positive_probability_vector(
    *,
    prob_mat: np.ndarray,
    classes: np.ndarray,
) -> np.ndarray | None:
    positive_idx = np.flatnonzero(classes == 1)
    if len(positive_idx) != 1:
        return None
    return prob_mat[:, int(positive_idx[0])]


def _compute_class_weights(
    *,
    model: _SamplingModel,
    prob_mat: np.ndarray,
    base_y: np.ndarray,
    sampling_policy: PredictionSamplingPolicy,
) -> dict[int, np.ndarray | None]:
    weights_by_class: dict[int, np.ndarray | None] = {}
    model_classes = np.asarray(model.classes_, dtype=int)  # type: ignore[attr-defined]
    for class_idx, class_label in enumerate(model_classes):
        class_mask = base_y == class_label
        class_prob = transform_resampling_probabilities(
            prob_mat[class_mask, class_idx],
            sampling_policy=sampling_policy,
        )
        weights_by_class[int(class_label)] = normalize_probabilities(class_prob)
    return weights_by_class


def _resample_training_rows(
    *,
    model: _SamplingModel,
    base_x: np.ndarray,
    base_y: np.ndarray,
    weights_by_class: dict[int, np.ndarray | None],
    resampling_rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    resampled_x: list[np.ndarray] = []
    resampled_y: list[np.ndarray] = []
    model_classes = np.asarray(model.classes_, dtype=int)  # type: ignore[attr-defined]

    for class_label in model_classes:
        class_mask = base_y == class_label
        class_x = base_x[class_mask]
        if class_x.shape[0] == 0:
            continue
        sample_prob = weights_by_class.get(int(class_label))
        sampled_idx = cast(
            np.ndarray,
            resampling_rng.choice(
                class_x.shape[0],
                size=class_x.shape[0],
                replace=True,
                p=sample_prob,
            ),
        )
        resampled_x.append(class_x[sampled_idx])
        resampled_y.append(np.repeat(class_label, class_x.shape[0]))

    return np.vstack(resampled_x), np.concatenate(resampled_y)


def _run_sampling_iterations(
    *,
    base_x: np.ndarray,
    base_y: np.ndarray,
    kernel: str,
    n_iterations: int,
    resampling_rng: np.random.Generator,
    sampling_policy: PredictionSamplingPolicy,
) -> _SamplingModel:
    if n_iterations < 1:
        raise ValueError("n_iterations must be at least 1")

    StandardScaler, SVC = require_sklearn()
    use_r_parity_scaler = sampling_policy.adaptive_policy == "r_parity"
    current_x = base_x
    current_y = base_y
    model: _SamplingModel | None = None

    for _ in range(n_iterations):
        model = cast(
            _SamplingModel,
            make_svm(
                StandardScaler=StandardScaler,
                SVC=SVC,
                kernel=kernel,
                use_r_parity_scaler=use_r_parity_scaler,
            ),
        )
        model.fit(current_x, current_y)
        prob_mat = np.asarray(model.predict_proba(base_x), dtype=float)
        weights_by_class = _compute_class_weights(
            model=model,
            prob_mat=prob_mat,
            base_y=base_y,
            sampling_policy=sampling_policy,
        )
        current_x, current_y = _resample_training_rows(
            model=model,
            base_x=base_x,
            base_y=base_y,
            weights_by_class=weights_by_class,
            resampling_rng=resampling_rng,
        )

    if model is None:
        raise RuntimeError("adaptive sampling failed to produce a fitted model")
    return model


def run_adaptive_sampling_ensemble(
    *,
    train_values: np.ndarray,
    train_labels: np.ndarray,
    test_values: np.ndarray,
    kernel: str,
    n_iterations: int,
    resampling_rng: np.random.Generator,
    sampling_policy: PredictionSamplingPolicy,
) -> np.ndarray:
    """Run one adaptive ensemble and return per-site scores."""

    base_x = np.asarray(train_values, dtype=float)
    base_y = np.asarray(train_labels, dtype=int)
    test_x = np.asarray(test_values, dtype=float)
    model = _run_sampling_iterations(
        base_x=base_x,
        base_y=base_y,
        kernel=kernel,
        n_iterations=n_iterations,
        resampling_rng=resampling_rng,
        sampling_policy=sampling_policy,
    )
    pred = np.asarray(model.predict_proba(test_x), dtype=float)
    model_classes = np.asarray(model.classes_, dtype=int)
    positive_probabilities = _positive_probability_vector(
        prob_mat=pred,
        classes=model_classes,
    )
    final_decision_vector = aligned_binary_decision_vector(
        model=model,
        values=test_x,
        positive_probabilities=positive_probabilities,
    )
    if sampling_policy.final_score_mode == "mean_probability":
        if positive_probabilities is None:
            msg = (
                "Expected exactly one positive class labelled 1; found "
                f"{model_classes.tolist()}"
            )
            raise ValueError(msg)
        return np.asarray(positive_probabilities, dtype=float)
    return np.asarray(1.0 / (1.0 + np.exp(-final_decision_vector)), dtype=float)


__all__ = ["run_adaptive_sampling_ensemble"]
