from __future__ import annotations

import builtins
import warnings

import numpy as np
import pytest

from phospy.science.prediction.svm import (
    aligned_binary_decision_vector,
    require_sklearn,
)


def test_require_sklearn_error_reports_broken_standard_install(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_import = builtins.__import__

    def failing_import(name: str, *args: object, **kwargs: object):
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
