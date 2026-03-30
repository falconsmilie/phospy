from __future__ import annotations

import numpy as np
import pandas as pd

from ..types import PredictionSvmMode
from ..validation.errors import InputCompatibilityError


class _RLikeStandardScaler:
    """Match R scale() semantics used by e1071::svm scaling."""

    def __init__(self) -> None:
        self.mean_: np.ndarray | None = None
        self.scale_: np.ndarray | None = None

    @staticmethod
    def __sklearn_tags__():
        from sklearn.utils import InputTags, Tags, TargetTags, TransformerTags

        return Tags(
            estimator_type="transformer",
            target_tags=TargetTags(required=False),
            transformer_tags=TransformerTags(),
            input_tags=InputTags(two_d_array=True),
        )

    @staticmethod
    def get_params(deep: bool = True) -> dict[str, object]:
        del deep
        return {}

    def set_params(self, **params: object) -> _RLikeStandardScaler:
        if params:
            msg = f"_RLikeStandardScaler does not accept parameters: {sorted(params)}"
            raise InputCompatibilityError(msg)
        return self

    def fit(
        self,
        values: np.ndarray,
        y: np.ndarray | None = None,
    ) -> _RLikeStandardScaler:
        del y
        x = np.asarray(values, dtype=float)
        self.mean_ = x.mean(axis=0)
        if x.shape[0] > 1:
            scale = x.std(axis=0, ddof=1)
        else:
            scale = np.ones(x.shape[1], dtype=float)
        self.scale_ = np.where(np.isfinite(scale) & (scale > 0.0), scale, 1.0)
        return self

    def transform(self, values: np.ndarray) -> np.ndarray:
        if self.mean_ is None or self.scale_ is None:
            msg = "_RLikeStandardScaler must be fitted before transform()"
            raise ValueError(msg)
        x = np.asarray(values, dtype=float)
        return (x - self.mean_) / self.scale_

    def fit_transform(
        self,
        values: np.ndarray,
        y: np.ndarray | None = None,
    ) -> np.ndarray:
        return self.fit(values, y=y).transform(values)


def aligned_binary_decision_values(
    *,
    model,
    values: np.ndarray,
    index: pd.Index,
    positive_probabilities: pd.Series | None,
) -> pd.Series:
    """Return binary decision values aligned so larger means more class 1-like."""

    decision_values = np.asarray(model.decision_function(values), dtype=float).reshape(
        -1
    )
    series = pd.Series(decision_values, index=index, dtype=float)
    if positive_probabilities is None:
        return series
    aligned_probabilities = positive_probabilities.reindex(index)
    corr = series.corr(aligned_probabilities, method="pearson")
    if pd.notna(corr) and corr < 0:
        return -series
    return series


def make_svm(
    *,
    StandardScaler: type,
    SVC: type,
    kernel: str,
    svm_mode: PredictionSvmMode,
):
    from sklearn.pipeline import make_pipeline

    scaler = StandardScaler() if svm_mode == "default" else _RLikeStandardScaler()
    gamma: str | float = "scale" if svm_mode == "default" else "auto"
    random_state = resolve_svm_probability_random_state()
    return make_pipeline(
        scaler,
        SVC(
            kernel=kernel,
            gamma=gamma,
            probability=True,
            random_state=random_state,
        ),
    )


def resolve_svm_probability_random_state() -> int:
    """Return the deterministic SVM probability-calibration seed."""

    return 1


def extract_svm_probability_parameters(model) -> pd.DataFrame | None:
    """Return libsvm Platt-scaling parameters from the fitted SVC step."""

    svc = model.steps[-1][1]
    prob_a = np.asarray(getattr(svc, "probA_", np.asarray([])), dtype=float)
    prob_b = np.asarray(getattr(svc, "probB_", np.asarray([])), dtype=float)
    if prob_a.size == 0 or prob_b.size == 0:
        return None

    classes = [str(class_label) for class_label in getattr(svc, "classes_", [])]
    if len(classes) >= 2:
        class_pairs = [
            f"{left}|{right}"
            for idx, left in enumerate(classes[:-1])
            for right in classes[idx + 1 :]
        ]
    else:
        class_pairs = []
    if len(class_pairs) != prob_a.size:
        class_pairs = [str(index + 1) for index in range(prob_a.size)]

    return pd.DataFrame(
        {
            "class_pair": class_pairs,
            "probA": prob_a.astype(float),
            "probB": prob_b.astype(float),
        }
    )


def require_sklearn() -> tuple[type, type]:
    try:
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
    except ImportError as exc:  # pragma: no cover - environment-dependent
        msg = (
            "KinasePredictor requires scikit-learn, but it could not be imported. "
            "Install the project dependencies or add scikit-learn manually."
        )
        raise ImportError(msg) from exc
    return StandardScaler, SVC


__all__ = [
    "aligned_binary_decision_values",
    "extract_svm_probability_parameters",
    "make_svm",
    "require_sklearn",
]
