from __future__ import annotations

import builtins
import warnings
from typing import Any

import numpy as np
import pytest

from phospy.science.prediction.policies import resolve_prediction_sampling_policy
from phospy.science.prediction.sampling_core import run_adaptive_sampling_ensemble
from phospy.science.prediction.svm import (
    aligned_binary_decision_vector,
    require_sklearn,
)


def test_require_sklearn_error_reports_broken_standard_install(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_import = builtins.__import__

    def failing_import(name: str, *args: Any, **kwargs: Any) -> object:
        if name == "sklearn" or name.startswith("sklearn."):
            raise ImportError("simulated sklearn import failure")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", failing_import)

    with pytest.raises(
        ImportError,
        match="part of PhosPy's standard install.*unexpected environment problem",
    ):
        require_sklearn()


def test_aligned_binary_decision_vector_skips_corrcoef_on_constant_vectors() -> None:
    class _ConstantDecisionModel:
        @staticmethod
        def decision_function(values: np.ndarray) -> np.ndarray:
            return np.ones(values.shape[0], dtype=float)

    values = np.asarray([[1.0], [2.0], [3.0]], dtype=float)
    probabilities = np.asarray([0.2, 0.5, 0.8], dtype=float)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        decision = aligned_binary_decision_vector(
            model=_ConstantDecisionModel(),
            values=values,
            positive_probabilities=probabilities,
        )

    assert decision.tolist() == [1.0, 1.0, 1.0]
    assert not any(issubclass(entry.category, RuntimeWarning) for entry in caught)


def test_adaptive_sampling_suppresses_known_sklearn_probability_warning() -> None:
    feature_values = np.asarray(
        [
            [0.95, 0.1],
            [0.9, 0.2],
            [0.2, 0.9],
            [0.1, 0.95],
        ],
        dtype=float,
    )
    labels = np.asarray([1, 1, 2, 2], dtype=int)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        scores = run_adaptive_sampling_ensemble(
            train_values=feature_values,
            train_labels=labels,
            test_values=feature_values,
            kernel="rbf",
            n_iterations=2,
            resampling_rng=np.random.default_rng(13),
            sampling_policy=resolve_prediction_sampling_policy("stable"),
        )

    assert scores.shape == (4,)
    assert not any(
        issubclass(entry.category, FutureWarning)
        and "`probability` parameter was deprecated" in str(entry.message)
        for entry in caught
    )
