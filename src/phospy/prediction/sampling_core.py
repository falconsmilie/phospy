from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..types import PredictionSvmMode, PredictionTraceLevel
from .models import AdaptiveSamplingEnsembleTrace, SamplingTraceOverrideEnsemble
from .sampling_runtime import (
    normalize_probabilities,
    resolve_sampled_site_positions,
    transform_resampling_probabilities,
)
from .sampling_trace_writer import write_final_trace_rows, write_iteration_trace_rows
from .svm import (
    aligned_binary_decision_values,
    extract_svm_probability_parameters,
    make_svm,
    require_sklearn,
)
from .trace_runtime import TraceSink


def _fit_sampling_model(
    *,
    current_x: np.ndarray,
    current_y: np.ndarray,
    kernel: str,
    svm_mode: PredictionSvmMode,
    StandardScaler: type[Any],
    SVC: type[Any],
) -> Any:
    model = make_svm(
        StandardScaler=StandardScaler,
        SVC=SVC,
        kernel=kernel,
        svm_mode=svm_mode,
    )
    model.fit(current_x, current_y)
    return model


def _extract_iteration_probabilities(
    *,
    model: Any,
    base_x: np.ndarray,
    base_y: np.ndarray,
    base_index: pd.Index,
) -> tuple[np.ndarray, pd.DataFrame, pd.Series, pd.DataFrame | None, pd.Series]:
    prob_mat = model.predict_proba(base_x)
    prob_df = pd.DataFrame(
        prob_mat,
        index=base_index,
        columns=[str(class_label) for class_label in model.classes_],
    )
    decision_series = aligned_binary_decision_values(
        model=model,
        values=base_x,
        index=base_index,
        positive_probabilities=prob_df.get("1"),
    )
    probability_parameters = extract_svm_probability_parameters(model)
    label_series = pd.Series(base_y, index=base_index, dtype=int)
    return prob_mat, prob_df, decision_series, probability_parameters, label_series


def _compute_class_weights(
    *,
    model: Any,
    prob_mat: np.ndarray,
    base_y: np.ndarray,
    base_index: pd.Index,
    svm_mode: PredictionSvmMode,
) -> dict[int, pd.Series | None]:
    weights_by_class: dict[int, pd.Series | None] = {}
    for class_idx, class_label in enumerate(model.classes_):
        class_mask = base_y == class_label
        class_index = base_index[class_mask]
        class_prob = transform_resampling_probabilities(
            prob_mat[class_mask, class_idx],
            svm_mode=svm_mode,
        )
        sample_prob = normalize_probabilities(class_prob)
        weights_by_class[int(class_label)] = (
            pd.Series(sample_prob, index=class_index, dtype=float)
            if sample_prob is not None
            else None
        )
    return weights_by_class


def _resample_training_rows(
    *,
    model: Any,
    base_x: np.ndarray,
    base_y: np.ndarray,
    base_index: pd.Index,
    weights_by_class: dict[int, pd.Series | None],
    resampling_rng: np.random.Generator,
    sampling_override: SamplingTraceOverrideEnsemble | None,
    iteration_index: int,
    ensemble_index: int,
) -> tuple[np.ndarray, np.ndarray, dict[int, list[str]]]:
    resampled_x: list[np.ndarray] = []
    resampled_y: list[np.ndarray] = []
    sampled_sites_by_class: dict[int, list[str]] = {}

    for class_label in model.classes_:
        class_mask = base_y == class_label
        class_x = base_x[class_mask]
        class_index = base_index[class_mask]
        sample_weights = weights_by_class.get(int(class_label))
        sample_prob = (
            None if sample_weights is None else sample_weights.to_numpy(dtype=float)
        )

        override_sites = None
        if sampling_override is not None:
            override_sites = sampling_override.get_iteration_sample_sites(
                iteration_index=iteration_index,
                class_label=int(class_label),
            )
        if override_sites is not None:
            sampled_idx = resolve_sampled_site_positions(
                available_sites=class_index,
                sampled_sites=override_sites,
                expected_size=class_x.shape[0],
                context=(
                    f"iteration={iteration_index}, ensemble={ensemble_index}, "
                    f"class_label={int(class_label)}"
                ),
            )
        else:
            sampled_idx = resampling_rng.choice(
                class_x.shape[0],
                size=class_x.shape[0],
                replace=True,
                p=sample_prob,
            )

        sampled_sites_by_class[int(class_label)] = class_index[sampled_idx].tolist()
        resampled_x.append(class_x[sampled_idx])
        resampled_y.append(np.repeat(class_label, class_x.shape[0]))

    return np.vstack(resampled_x), np.concatenate(resampled_y), sampled_sites_by_class


def multi_ada_sampling(
    train_mat: pd.DataFrame,
    test_mat: pd.DataFrame,
    labels: np.ndarray,
    kernel: str,
    n_iterations: int,
    resampling_rng: np.random.Generator,
    capture_trace: bool,
    trace_level: PredictionTraceLevel,
    trace_sink: TraceSink | None,
    kinase: str,
    ensemble_index: int,
    initial_negative_sites: list[str],
    debug_top_n: int,
    svm_mode: PredictionSvmMode,
    sampling_override: SamplingTraceOverrideEnsemble | None,
) -> tuple[pd.Series, AdaptiveSamplingEnsembleTrace | None]:
    StandardScaler, SVC = require_sklearn()

    base_x = train_mat.to_numpy(dtype=float)
    base_y = np.asarray(labels, dtype=int)
    base_index = train_mat.index.copy()
    current_x = base_x
    current_y = base_y
    if capture_trace and trace_level == "full" and trace_sink is None:
        msg = "full trace capture requires a trace sink"
        raise ValueError(msg)

    model = None

    for iteration_index in range(1, n_iterations + 1):
        model = _fit_sampling_model(
            current_x=current_x,
            current_y=current_y,
            kernel=kernel,
            svm_mode=svm_mode,
            StandardScaler=StandardScaler,
            SVC=SVC,
        )
        (
            prob_mat,
            prob_df,
            decision_series,
            probability_parameters,
            label_series,
        ) = _extract_iteration_probabilities(
            model=model,
            base_x=base_x,
            base_y=base_y,
            base_index=base_index,
        )
        weights_by_class = _compute_class_weights(
            model=model,
            prob_mat=prob_mat,
            base_y=base_y,
            base_index=base_index,
            svm_mode=svm_mode,
        )
        current_x, current_y, sampled_sites_by_class = _resample_training_rows(
            model=model,
            base_x=base_x,
            base_y=base_y,
            base_index=base_index,
            weights_by_class=weights_by_class,
            resampling_rng=resampling_rng,
            sampling_override=sampling_override,
            iteration_index=iteration_index,
            ensemble_index=ensemble_index,
        )

        if capture_trace and trace_level == "full" and trace_sink is not None:
            write_iteration_trace_rows(
                trace_sink=trace_sink,
                kinase=kinase,
                ensemble_index=ensemble_index,
                iteration_index=iteration_index,
                labels=label_series,
                probabilities=prob_df,
                probability_parameters=probability_parameters,
                decision_values=decision_series,
                weights_by_class=weights_by_class,
                sampled_sites_by_class=sampled_sites_by_class,
            )

    if model is None:
        msg = "n_iterations must be at least 1"
        raise ValueError(msg)

    test_x = test_mat.to_numpy(dtype=float)
    pred = model.predict_proba(test_x)
    pred_df = pd.DataFrame(
        pred,
        index=test_mat.index.copy(),
        columns=[str(class_label) for class_label in model.classes_],
    )
    positive_idx = np.flatnonzero(model.classes_ == 1)
    if len(positive_idx) != 1:
        msg = (
            "Expected exactly one positive class labelled 1; found "
            f"{model.classes_.tolist()}"
        )
        raise ValueError(msg)
    positive_series = pd.Series(
        pred[:, positive_idx[0]], index=test_mat.index.copy(), dtype=float
    )

    ensemble_trace = None
    if capture_trace:
        final_top_sites = (
            positive_series.sort_values(ascending=False)
            .head(debug_top_n)
            .index.tolist()
        )
        final_decision_values = aligned_binary_decision_values(
            model=model,
            values=test_x,
            index=test_mat.index.copy(),
            positive_probabilities=pred_df.get("1"),
        )
        if trace_level == "full":
            write_final_trace_rows(
                trace_sink=trace_sink,
                kinase=kinase,
                ensemble_index=ensemble_index,
                final_probabilities=pred_df,
                final_decision_values=final_decision_values,
                final_top_sites=final_top_sites,
            )

        ensemble_trace = AdaptiveSamplingEnsembleTrace(
            ensemble_index=ensemble_index,
            initial_negative_sites=list(initial_negative_sites),
            final_top_sites=final_top_sites,
        )

    return positive_series, ensemble_trace


__all__ = [
    "multi_ada_sampling",
]
