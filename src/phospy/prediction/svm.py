"""SVM helpers for adaptive prediction sampling."""

from __future__ import annotations

import numpy as np


class _RLikeStandardScaler:
    """Match R `scale()` semantics used by e1071::svm scaling."""

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
            raise ValueError(msg)
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


def aligned_binary_decision_vector(
    *,
    model: object,
    values: np.ndarray,
    positive_probabilities: np.ndarray | None,
) -> np.ndarray:
    """Return binary decision values aligned so larger means more class-1-like."""

    decision_values = np.asarray(
        model.decision_function(values),  # type: ignore[attr-defined]
        dtype=float,
    ).reshape(-1)
    if positive_probabilities is None:
        return decision_values

    aligned_probabilities = np.asarray(positive_probabilities, dtype=float).reshape(-1)
    finite_mask = np.isfinite(decision_values) & np.isfinite(aligned_probabilities)
    if finite_mask.sum() < 2:
        return decision_values

    corr = np.corrcoef(
        decision_values[finite_mask],
        aligned_probabilities[finite_mask],
    )[0, 1]
    if np.isfinite(corr) and corr < 0:
        return -decision_values
    return decision_values


def make_svm(
    *,
    StandardScaler: type,
    SVC: type,
    kernel: str,
    use_r_parity_scaler: bool,
) -> object:
    """Build one configured SVM pipeline used in adaptive sampling."""

    from sklearn.pipeline import make_pipeline

    scaler = StandardScaler() if not use_r_parity_scaler else _RLikeStandardScaler()
    gamma: str | float = "scale" if not use_r_parity_scaler else "auto"
    return make_pipeline(
        scaler,
        SVC(
            kernel=kernel,
            gamma=gamma,
            probability=True,
            random_state=1,
        ),
    )


def require_sklearn() -> tuple[type, type]:
    """Import and return scikit-learn classes needed by adaptive sampling."""

    try:
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
    except ImportError as exc:  # pragma: no cover - environment-dependent
        msg = (
            "adaptive ensemble prediction requires scikit-learn, but it could not "
            "be imported. Install scikit-learn to use prediction mode "
            "'adaptive_ensemble'."
        )
        raise ImportError(msg) from exc
    return StandardScaler, SVC


__all__ = [
    "aligned_binary_decision_vector",
    "make_svm",
    "require_sklearn",
]
