from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..internal.types import (
    PREDICTION_TRACE_LEVEL_FULL,
    PredictionSvmMode,
    PredictionTraceLevel,
)
from .policies import PredictionSamplingPolicy, resolve_prediction_sampling_policy
from .results import AdaptiveSamplingEnsembleTrace, SamplingTraceOverrideEnsemble
from .sampling_runtime import (
    normalize_probabilities,
    resolve_sampled_site_positions,
    transform_resampling_probabilities,
)
from .sampling_trace_writer import write_final_trace_rows, write_iteration_trace_rows
from .svm import (
    aligned_binary_decision_vector,
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


def _build_probability_frame(
    *,
    prob_mat: np.ndarray,
    index: pd.Index,
    classes: np.ndarray,
) -> pd.DataFrame:
    return pd.DataFrame(
        prob_mat,
        index=index,
        columns=[str(class_label) for class_label in classes],
    )


def _resolve_matrix_values(
    *,
    frame: pd.DataFrame | None,
    values: np.ndarray | None,
    index: pd.Index | None,
    context: str,
) -> tuple[np.ndarray, pd.Index]:
    """Resolve numeric matrix values and aligned row index for sampling.

    The public DataFrame path remains supported, but internal prediction callers
    can pass precomputed NumPy arrays and indexes to avoid repeated DataFrame
    materialisation inside the ensemble loop.
    """

    if values is None:
        if frame is None:
            msg = f"{context} requires either a DataFrame or precomputed values"
            raise ValueError(msg)
        return frame.to_numpy(dtype=float), frame.index.copy()

    if index is None:
        msg = f"{context} precomputed values require a matching index"
        raise ValueError(msg)

    resolved_values = np.asarray(values, dtype=float)
    if resolved_values.ndim != 2:
        msg = f"{context} precomputed values must be two-dimensional"
        raise ValueError(msg)
    if resolved_values.shape[0] != len(index):
        msg = f"{context} precomputed values row count must match the provided index"
        raise ValueError(msg)
    return resolved_values, index.copy()


def _positive_probability_vector(
    *,
    prob_mat: np.ndarray,
    classes: np.ndarray,
) -> np.ndarray | None:
    positive_idx = np.flatnonzero(classes == 1)
    if len(positive_idx) != 1:
        return None
    return prob_mat[:, int(positive_idx[0])]


def _extract_iteration_sampling_state(
    *,
    model: Any,
    base_x: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    prob_mat = model.predict_proba(base_x)
    decision_values = aligned_binary_decision_vector(
        model=model,
        values=base_x,
        positive_probabilities=_positive_probability_vector(
            prob_mat=prob_mat,
            classes=np.asarray(model.classes_),
        ),
    )
    return prob_mat, decision_values


def _extract_iteration_trace_payload(
    *,
    model: Any,
    prob_mat: np.ndarray,
    decision_values: np.ndarray,
    base_y: np.ndarray,
    base_index: pd.Index,
) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame | None, pd.Series]:
    prob_df = _build_probability_frame(
        prob_mat=prob_mat,
        index=base_index,
        classes=np.asarray(model.classes_),
    )
    decision_series = pd.Series(decision_values, index=base_index, dtype=float)
    probability_parameters = extract_svm_probability_parameters(model)
    label_series = pd.Series(base_y, index=base_index, dtype=int)
    return prob_df, decision_series, probability_parameters, label_series


def _compute_class_weights(
    *,
    model: Any,
    prob_mat: np.ndarray,
    base_y: np.ndarray,
    base_index: pd.Index,
    svm_mode: PredictionSvmMode,
    sampling_policy: PredictionSamplingPolicy | None,
) -> dict[int, pd.Series | None]:
    resolved_sampling_policy = sampling_policy or resolve_prediction_sampling_policy(
        svm_mode
    )
    weights_by_class: dict[int, pd.Series | None] = {}
    for class_idx, class_label in enumerate(model.classes_):
        class_mask = base_y == class_label
        class_index = base_index[class_mask]
        class_prob = transform_resampling_probabilities(
            prob_mat[class_mask, class_idx],
            sampling_policy=resolved_sampling_policy,
        )
        sample_prob = normalize_probabilities(class_prob)
        weights_by_class[int(class_label)] = (
            pd.Series(sample_prob, index=class_index, dtype=float)
            if sample_prob is not None
            else None
        )
    return weights_by_class


def _resolve_final_score_series(
    *,
    pred_df: pd.DataFrame,
    final_decision_values: pd.Series,
    sampling_policy: PredictionSamplingPolicy,
) -> pd.Series:
    """Return the final per-site score series for the configured sampling mode.

    This helper remains the public review seam for tests and diagnostics. The
    hot prediction path now bypasses the surrounding DataFrame materialisation
    when full trace capture is disabled.
    """

    positive_probabilities = pd.Series(
        pred_df.loc[:, "1"].to_numpy(dtype=float, copy=False),
        index=pred_df.index.copy(),
        dtype=float,
    )
    if sampling_policy.final_score_mode == "mean_probability":
        return positive_probabilities

    decision_vector = final_decision_values.to_numpy(dtype=float, copy=False)
    return pd.Series(
        1.0 / (1.0 + np.exp(-decision_vector)),
        index=pred_df.index.copy(),
        dtype=float,
    )


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


def _validate_trace_configuration(
    *,
    capture_trace: bool,
    trace_level: PredictionTraceLevel,
    trace_sink: TraceSink | None,
) -> None:
    if (
        capture_trace
        and trace_level == PREDICTION_TRACE_LEVEL_FULL
        and trace_sink is None
    ):
        msg = "full trace capture requires a trace sink"
        raise ValueError(msg)


def _write_iteration_trace_if_requested(
    *,
    capture_trace: bool,
    trace_level: PredictionTraceLevel,
    trace_sink: TraceSink | None,
    model: Any,
    prob_mat: np.ndarray,
    decision_values: np.ndarray,
    base_y: np.ndarray,
    base_index: pd.Index,
    weights_by_class: dict[int, pd.Series | None],
    sampled_sites_by_class: dict[int, list[str]],
    kinase: str,
    ensemble_index: int,
    iteration_index: int,
) -> None:
    if (
        not capture_trace
        or trace_level != PREDICTION_TRACE_LEVEL_FULL
        or trace_sink is None
    ):
        return

    prob_df, decision_series, probability_parameters, label_series = (
        _extract_iteration_trace_payload(
            model=model,
            prob_mat=prob_mat,
            decision_values=decision_values,
            base_y=base_y,
            base_index=base_index,
        )
    )
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


def _run_sampling_iterations(
    *,
    base_x: np.ndarray,
    base_y: np.ndarray,
    base_index: pd.Index,
    kernel: str,
    n_iterations: int,
    resampling_rng: np.random.Generator,
    capture_trace: bool,
    trace_level: PredictionTraceLevel,
    trace_sink: TraceSink | None,
    kinase: str,
    ensemble_index: int,
    svm_mode: PredictionSvmMode,
    sampling_policy: PredictionSamplingPolicy,
    sampling_override: SamplingTraceOverrideEnsemble | None,
    StandardScaler: type[Any],
    SVC: type[Any],
) -> Any:
    current_x = base_x
    current_y = base_y
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
        prob_mat, decision_values = _extract_iteration_sampling_state(
            model=model,
            base_x=base_x,
        )
        weights_by_class = _compute_class_weights(
            model=model,
            prob_mat=prob_mat,
            base_y=base_y,
            base_index=base_index,
            svm_mode=svm_mode,
            sampling_policy=sampling_policy,
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
        _write_iteration_trace_if_requested(
            capture_trace=capture_trace,
            trace_level=trace_level,
            trace_sink=trace_sink,
            model=model,
            prob_mat=prob_mat,
            decision_values=decision_values,
            base_y=base_y,
            base_index=base_index,
            weights_by_class=weights_by_class,
            sampled_sites_by_class=sampled_sites_by_class,
            kinase=kinase,
            ensemble_index=ensemble_index,
            iteration_index=iteration_index,
        )

    if model is None:
        msg = "n_iterations must be at least 1"
        raise ValueError(msg)
    return model


def _resolve_final_score_values(
    *,
    model: Any,
    test_mat: pd.DataFrame | None,
    test_values: np.ndarray | None,
    test_index: pd.Index | None,
    sampling_policy: PredictionSamplingPolicy,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.Index]:
    test_x, resolved_test_index = _resolve_matrix_values(
        frame=test_mat,
        values=test_values,
        index=test_index,
        context="test_mat",
    )
    pred = model.predict_proba(test_x)
    positive_probabilities = _positive_probability_vector(
        prob_mat=pred,
        classes=np.asarray(model.classes_),
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
                f"{np.asarray(model.classes_).tolist()}"
            )
            raise ValueError(msg)
        final_score_values = positive_probabilities
    else:
        final_score_values = 1.0 / (1.0 + np.exp(-final_decision_vector))
    return (
        np.asarray(final_score_values, dtype=float),
        final_decision_vector,
        pred,
        resolved_test_index,
    )


def _build_final_trace_artifacts(
    *,
    capture_trace: bool,
    trace_level: PredictionTraceLevel,
    trace_sink: TraceSink | None,
    kinase: str,
    ensemble_index: int,
    initial_negative_sites: list[str],
    debug_top_n: int,
    model: Any,
    pred: np.ndarray,
    final_decision_vector: np.ndarray,
    final_score_values: np.ndarray,
    resolved_test_index: pd.Index,
) -> tuple[pd.Series | None, AdaptiveSamplingEnsembleTrace | None]:
    if not capture_trace:
        return None, None

    final_score_series = pd.Series(
        final_score_values,
        index=resolved_test_index,
        dtype=float,
    )
    final_top_sites = (
        final_score_series.sort_values(ascending=False).head(debug_top_n).index.tolist()
    )

    if trace_level == PREDICTION_TRACE_LEVEL_FULL and trace_sink is not None:
        pred_df = _build_probability_frame(
            prob_mat=pred,
            index=resolved_test_index,
            classes=np.asarray(model.classes_),
        )
        final_decision_values = pd.Series(
            final_decision_vector,
            index=resolved_test_index,
            dtype=float,
        )
        write_final_trace_rows(
            trace_sink=trace_sink,
            kinase=kinase,
            ensemble_index=ensemble_index,
            final_probabilities=pred_df,
            final_decision_values=final_decision_values,
            final_top_sites=final_top_sites,
        )

    return final_score_series, AdaptiveSamplingEnsembleTrace(
        ensemble_index=ensemble_index,
        initial_negative_sites=list(initial_negative_sites),
        final_top_sites=final_top_sites,
    )


def multi_ada_sampling(
    train_mat: pd.DataFrame | None,
    test_mat: pd.DataFrame | None,
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
    sampling_policy: PredictionSamplingPolicy | None = None,
    sampling_override: SamplingTraceOverrideEnsemble | None = None,
    train_values: np.ndarray | None = None,
    train_index: pd.Index | None = None,
    test_values: np.ndarray | None = None,
    test_index: pd.Index | None = None,
    return_values: bool = False,
) -> tuple[np.ndarray | pd.Series, AdaptiveSamplingEnsembleTrace | None]:
    StandardScaler, SVC = require_sklearn()

    base_x, base_index = _resolve_matrix_values(
        frame=train_mat,
        values=train_values,
        index=train_index,
        context="train_mat",
    )
    base_y = np.asarray(labels, dtype=int)
    resolved_sampling_policy = sampling_policy or resolve_prediction_sampling_policy(
        svm_mode
    )
    _validate_trace_configuration(
        capture_trace=capture_trace,
        trace_level=trace_level,
        trace_sink=trace_sink,
    )
    model = _run_sampling_iterations(
        base_x=base_x,
        base_y=base_y,
        base_index=base_index,
        kernel=kernel,
        n_iterations=n_iterations,
        resampling_rng=resampling_rng,
        capture_trace=capture_trace,
        trace_level=trace_level,
        trace_sink=trace_sink,
        kinase=kinase,
        ensemble_index=ensemble_index,
        svm_mode=svm_mode,
        sampling_policy=resolved_sampling_policy,
        sampling_override=sampling_override,
        StandardScaler=StandardScaler,
        SVC=SVC,
    )
    final_score_values, final_decision_vector, pred, resolved_test_index = (
        _resolve_final_score_values(
            model=model,
            test_mat=test_mat,
            test_values=test_values,
            test_index=test_index,
            sampling_policy=resolved_sampling_policy,
        )
    )
    final_score_series, ensemble_trace = _build_final_trace_artifacts(
        capture_trace=capture_trace,
        trace_level=trace_level,
        trace_sink=trace_sink,
        kinase=kinase,
        ensemble_index=ensemble_index,
        initial_negative_sites=initial_negative_sites,
        debug_top_n=debug_top_n,
        model=model,
        pred=pred,
        final_decision_vector=final_decision_vector,
        final_score_values=final_score_values,
        resolved_test_index=resolved_test_index,
    )

    if return_values:
        return final_score_values, ensemble_trace

    if final_score_series is None:
        final_score_series = pd.Series(
            final_score_values,
            index=resolved_test_index,
            dtype=float,
        )
    return final_score_series, ensemble_trace


__all__ = [
    "multi_ada_sampling",
]
