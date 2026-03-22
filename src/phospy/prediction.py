from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from .scoring import KinaseScoringResult

NegativeSamplingStrategy = Literal["random", "coverage", "hybrid"]


@dataclass(slots=True)
class KinasePredictionResult:
    pred_matrix: pd.DataFrame
    substrate_list: dict[str, list[str]]


class KinasePredictor:
    """Predict kinase-substrate relationships from phosphosite score matrices.

    This is a narrow native Python port of PhosR's score-to-prediction seam.
    It mirrors the broad structure of ``kinaseSubstratePred()``: candidate
    substrates are selected from the combined score matrix, then an ensemble of
    adaptive SVM models is used to produce a kinase prediction matrix.
    """

    def __init__(self, kernel: str = "rbf") -> None:
        self.kernel = kernel

    def predict(
        self,
        combined_scores: pd.DataFrame,
        ensemble_size: int = 10,
        top: int = 50,
        score_threshold: float = 0.8,
        inclusion: int = 20,
        n_iterations: int = 5,
        random_state: int | None = None,
        negative_sampling_strategy: NegativeSamplingStrategy = "hybrid",
    ) -> KinasePredictionResult:
        _validate_positive_int(ensemble_size, name="ensemble_size")
        _validate_positive_int(top, name="top")
        _validate_positive_int(inclusion, name="inclusion")
        _validate_positive_int(n_iterations, name="n_iterations")
        _validate_negative_sampling_strategy(negative_sampling_strategy)

        substrate_list = build_candidate_substrate_list(
            combined_scores,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
        )
        if not substrate_list:
            empty = pd.DataFrame(index=combined_scores.index.copy(), dtype=float)
            return KinasePredictionResult(pred_matrix=empty, substrate_list={})

        rng = np.random.default_rng(random_state)
        feature_mat = combined_scores.astype(float)
        pred_matrix = pd.DataFrame(
            0.0,
            index=feature_mat.index.copy(),
            columns=list(substrate_list),
        )

        for kinase, substrates in substrate_list.items():
            positive_train = feature_mat.loc[substrates, :]
            negative_pool = feature_mat.loc[
                ~feature_mat.index.isin(substrates),
                :,
            ]
            if negative_pool.empty:
                msg = (
                    f"No negative pool sites available to train predictor for {kinase}"
                )
                raise ValueError(msg)

            negative_batches = _build_initial_negative_batches(
                negative_index=negative_pool.index.to_numpy(),
                batch_size=len(positive_train),
                ensemble_size=ensemble_size,
                rng=rng,
                strategy=negative_sampling_strategy,
            )
            for negative_indices in negative_batches:
                negative_train = negative_pool.loc[list(negative_indices), :]
                train_mat = pd.concat([positive_train, negative_train], axis=0)
                labels = np.concatenate(
                    [
                        np.repeat(1, len(positive_train)),
                        np.repeat(2, len(negative_train)),
                    ]
                )
                pred_matrix.loc[:, kinase] += _multi_ada_sampling(
                    train_mat=train_mat,
                    test_mat=feature_mat,
                    labels=labels,
                    kernel=self.kernel,
                    n_iterations=n_iterations,
                    rng=rng,
                )

        pred_matrix /= float(ensemble_size)
        return KinasePredictionResult(
            pred_matrix=pred_matrix,
            substrate_list=substrate_list,
        )

    def predict_from_scoring_result(
        self,
        scoring_result: KinaseScoringResult,
        ensemble_size: int = 10,
        top: int = 50,
        score_threshold: float = 0.8,
        inclusion: int = 20,
        n_iterations: int = 5,
        random_state: int | None = None,
        allow_profile_only_fallback: bool = False,
        negative_sampling_strategy: NegativeSamplingStrategy = "hybrid",
    ) -> KinasePredictionResult:
        if scoring_result.combined_scores is not None:
            feature_mat = scoring_result.combined_scores
        elif allow_profile_only_fallback:
            feature_mat = scoring_result.profile_scores
        else:
            msg = (
                "scoring_result does not contain combined_scores; pass "
                "allow_profile_only_fallback=True to use profile_scores instead"
            )
            raise ValueError(msg)

        return self.predict(
            combined_scores=feature_mat,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            negative_sampling_strategy=negative_sampling_strategy,
        )


def build_candidate_substrate_list(
    combined_scores: pd.DataFrame,
    top: int = 50,
    score_threshold: float = 0.8,
    inclusion: int = 20,
) -> dict[str, list[str]]:
    """Select candidate kinase substrates from the combined score matrix."""

    _validate_positive_int(top, name="top")
    _validate_positive_int(inclusion, name="inclusion")

    substrate_list: dict[str, list[str]] = {}
    for kinase in combined_scores.columns:
        series = combined_scores.loc[:, kinase]
        selected = (
            pd.DataFrame(
                {"site_id": series.index.astype(str), "score": series.to_numpy()}
            )
            .sort_values(
                ["score", "site_id"], ascending=[False, True], kind="mergesort"
            )
            .head(top)
        )
        sites = selected.loc[selected["score"] > score_threshold, "site_id"].tolist()
        if len(sites) >= inclusion:
            substrate_list[kinase] = sites
    return substrate_list


def _build_initial_negative_batches(
    *,
    negative_index: np.ndarray,
    batch_size: int,
    ensemble_size: int,
    rng: np.random.Generator,
    strategy: NegativeSamplingStrategy,
) -> list[np.ndarray]:
    if strategy == "random":
        return [
            rng.choice(
                negative_index,
                size=batch_size,
                replace=len(negative_index) < batch_size,
            )
            for _ in range(ensemble_size)
        ]
    if strategy == "coverage":
        return _build_coverage_negative_batches(
            negative_index=negative_index,
            batch_size=batch_size,
            ensemble_size=ensemble_size,
            rng=rng,
        )
    if strategy == "hybrid":
        n_coverage = ensemble_size // 2
        n_random = ensemble_size - n_coverage
        coverage_batches = _build_coverage_negative_batches(
            negative_index=negative_index,
            batch_size=batch_size,
            ensemble_size=n_coverage,
            rng=rng,
        )
        random_batches = [
            rng.choice(
                negative_index,
                size=batch_size,
                replace=len(negative_index) < batch_size,
            )
            for _ in range(n_random)
        ]
        return [*coverage_batches, *random_batches]
    raise ValueError(f"Unsupported negative_sampling_strategy: {strategy}")


def _build_coverage_negative_batches(
    *,
    negative_index: np.ndarray,
    batch_size: int,
    ensemble_size: int,
    rng: np.random.Generator,
) -> list[np.ndarray]:
    shuffled = rng.permutation(np.asarray(negative_index, dtype=object))
    if shuffled.size == 0:
        return []

    cursor = 0
    batches: list[np.ndarray] = []
    for _ in range(ensemble_size):
        pieces: list[np.ndarray] = []
        needed = batch_size
        while needed > 0:
            if cursor >= shuffled.size:
                shuffled = rng.permutation(np.asarray(negative_index, dtype=object))
                cursor = 0
            take_n = min(needed, shuffled.size - cursor)
            pieces.append(shuffled[cursor : cursor + take_n])
            cursor += take_n
            needed -= take_n
        batches.append(np.concatenate(pieces))
    return batches


def _multi_ada_sampling(
    train_mat: pd.DataFrame,
    test_mat: pd.DataFrame,
    labels: np.ndarray,
    kernel: str,
    n_iterations: int,
    rng: np.random.Generator,
) -> pd.Series:
    StandardScaler, SVC = _require_sklearn()

    base_x = train_mat.to_numpy(dtype=float)
    base_y = np.asarray(labels, dtype=int)
    current_x = base_x
    current_y = base_y
    model = None

    for _ in range(n_iterations):
        model = _make_svm(
            StandardScaler=StandardScaler, SVC=SVC, kernel=kernel, rng=rng
        )
        model.fit(current_x, current_y)
        prob_mat = model.predict_proba(base_x)

        resampled_x: list[np.ndarray] = []
        resampled_y: list[np.ndarray] = []
        for class_idx, class_label in enumerate(model.classes_):
            class_mask = base_y == class_label
            class_x = base_x[class_mask]
            class_prob = prob_mat[class_mask, class_idx]
            sample_prob = _normalize_probabilities(class_prob)
            sampled_idx = rng.choice(
                class_x.shape[0],
                size=class_x.shape[0],
                replace=True,
                p=sample_prob,
            )
            resampled_x.append(class_x[sampled_idx])
            resampled_y.append(np.repeat(class_label, class_x.shape[0]))

        current_x = np.vstack(resampled_x)
        current_y = np.concatenate(resampled_y)

    if model is None:
        msg = "n_iterations must be at least 1"
        raise ValueError(msg)

    pred = model.predict_proba(test_mat.to_numpy(dtype=float))
    positive_idx = np.flatnonzero(model.classes_ == 1)
    if len(positive_idx) != 1:
        msg = f"Expected exactly one positive class labelled 1; found {model.classes_.tolist()}"
        raise ValueError(msg)
    return pd.Series(pred[:, positive_idx[0]], index=test_mat.index.copy(), dtype=float)


def _make_svm(
    *, StandardScaler: type, SVC: type, kernel: str, rng: np.random.Generator
):
    from sklearn.pipeline import make_pipeline

    return make_pipeline(
        StandardScaler(),
        SVC(
            kernel=kernel,
            probability=True,
            random_state=int(rng.integers(0, 2**31 - 1)),
        ),
    )


def _normalize_probabilities(values: np.ndarray) -> np.ndarray | None:
    total = float(np.nansum(values))
    if total <= 0.0 or not np.isfinite(total):
        return None
    return values / total


def _require_sklearn() -> tuple[type, type]:
    try:
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
    except ImportError as exc:  # pragma: no cover - environment-dependent
        msg = (
            "KinasePredictor requires scikit-learn. Install the package with "
            "the 'ml' extra, for example: pip install .[ml]"
        )
        raise ImportError(msg) from exc
    return StandardScaler, SVC


def _validate_positive_int(value: int, name: str) -> None:
    if value < 1:
        raise ValueError(f"{name} must be at least 1")


def _validate_negative_sampling_strategy(value: str) -> None:
    if value not in {"random", "coverage", "hybrid"}:
        raise ValueError(
            "negative_sampling_strategy must be one of: random, coverage, hybrid"
        )
