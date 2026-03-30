from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..types import PredictionSvmMode, PredictionTraceLevel
from .models import (
    AdaptiveSamplingEnsembleTrace,
    AdaptiveSamplingIterationTrace,
    SamplingTraceOverrideEnsemble,
)
from .svm import (
    aligned_binary_decision_values,
    extract_svm_probability_parameters,
    make_svm,
    require_sklearn,
)
from .traces import PredictionSamplingTrace, TraceSink


def make_prediction_random_generators(
    rng: np.random.Generator,
) -> tuple[np.random.Generator, np.random.Generator]:
    """Create independent RNG streams for prediction sampling steps."""

    return (
        np.random.default_rng(int(rng.integers(0, 2**31 - 1))),
        np.random.default_rng(int(rng.integers(0, 2**31 - 1))),
    )


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


def _build_iteration_trace(
    *,
    iteration_index: int,
    labels: pd.Series,
    probabilities: pd.DataFrame,
    probability_parameters: pd.DataFrame | None,
    decision_values: pd.Series,
    weights_by_class: dict[int, pd.Series | None],
    sampled_sites_by_class: dict[int, list[str]],
) -> AdaptiveSamplingIterationTrace:
    return AdaptiveSamplingIterationTrace(
        iteration_index=iteration_index,
        labels=labels,
        probabilities=probabilities,
        probability_parameters=probability_parameters,
        decision_values=decision_values,
        positive_weights=weights_by_class.get(1),
        negative_weights=weights_by_class.get(2),
        sampled_positive_sites=sampled_sites_by_class.get(1, []),
        sampled_negative_sites=sampled_sites_by_class.get(2, []),
    )


def _write_iteration_trace_rows(
    *,
    trace_sink: TraceSink,
    kinase: str,
    ensemble_index: int,
    iteration_index: int,
    labels: pd.Series,
    probabilities: pd.DataFrame,
    probability_parameters: pd.DataFrame | None,
    decision_values: pd.Series,
    weights_by_class: dict[int, pd.Series | None],
    sampled_sites_by_class: dict[int, list[str]],
) -> None:
    label_rows: list[dict[str, object]] = []
    probability_rows: list[dict[str, object]] = []
    decision_rows: list[dict[str, object]] = []
    for position, site in enumerate(probabilities.index.tolist()):
        label_value = int(labels.iloc[position])
        prob_row = probabilities.iloc[position]
        decision_value = float(decision_values.iloc[position])
        label_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "site": site,
                "label": label_value,
            }
        )
        probability_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "site": site,
                "label": label_value,
                "prob_class_1": float(prob_row.get("1", float("nan"))),
                "prob_class_2": float(prob_row.get("2", float("nan"))),
            }
        )
        decision_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "site": site,
                "label": label_value,
                "decision_value_class_1": decision_value,
            }
        )
    trace_sink.write_rows("trace_iteration_labels", label_rows)
    trace_sink.write_rows("trace_iteration_probabilities", probability_rows)
    trace_sink.write_rows("trace_iteration_decision_values", decision_rows)

    if probability_parameters is not None:
        trace_sink.write_rows(
            "trace_iteration_probability_parameters",
            [
                {
                    "kinase": kinase,
                    "ensemble": ensemble_index,
                    "iteration": iteration_index,
                    "class_pair": str(row["class_pair"]),
                    "probA": float(row["probA"]),
                    "probB": float(row["probB"]),
                }
                for _, row in probability_parameters.iterrows()
            ],
        )

    weight_rows: list[dict[str, object]] = []
    for class_label, weights in (
        (1, weights_by_class.get(1)),
        (2, weights_by_class.get(2)),
    ):
        if weights is None:
            continue
        for site, weight in weights.items():
            weight_rows.append(
                {
                    "kinase": kinase,
                    "ensemble": ensemble_index,
                    "iteration": iteration_index,
                    "class_label": class_label,
                    "site": site,
                    "normalized_weight": float(weight),
                }
            )
    trace_sink.write_rows("trace_iteration_resampling_weights", weight_rows)

    sample_rows: list[dict[str, object]] = []
    for class_label in (1, 2):
        for draw, site in enumerate(
            sampled_sites_by_class.get(class_label, []), start=1
        ):
            sample_rows.append(
                {
                    "kinase": kinase,
                    "ensemble": ensemble_index,
                    "iteration": iteration_index,
                    "class_label": class_label,
                    "draw": draw,
                    "site": site,
                }
            )
    trace_sink.write_rows("trace_iteration_samples", sample_rows)


def _write_final_trace_rows(
    *,
    trace_sink: TraceSink,
    kinase: str,
    ensemble_index: int,
    final_probabilities: pd.DataFrame,
    final_decision_values: pd.Series,
    final_top_sites: list[str],
) -> None:
    prediction_rows: list[dict[str, object]] = []
    decision_rows: list[dict[str, object]] = []
    top_rows: list[dict[str, object]] = []
    for site in final_probabilities.index:
        prediction_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "site": site,
                "prob_class_1": float(final_probabilities.loc[site, "1"])
                if "1" in final_probabilities.columns
                else float("nan"),
                "prob_class_2": float(final_probabilities.loc[site, "2"])
                if "2" in final_probabilities.columns
                else float("nan"),
            }
        )
        decision_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "site": site,
                "decision_value_class_1": float(final_decision_values.loc[site]),
            }
        )
    for rank, site in enumerate(final_top_sites, start=1):
        top_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "rank": rank,
                "site": site,
                "prob_class_1": float(final_probabilities.loc[site, "1"])
                if "1" in final_probabilities.columns
                else float("nan"),
            }
        )
    trace_sink.write_rows("trace_final_ensemble_predictions", prediction_rows)
    trace_sink.write_rows("trace_final_ensemble_decision_values", decision_rows)
    trace_sink.write_rows("trace_final_ensemble_top", top_rows)


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
    model = None
    iteration_traces: list[AdaptiveSamplingIterationTrace] = []

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
            _write_iteration_trace_rows(
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
        if trace_level == "full" and trace_sink is not None:
            _write_final_trace_rows(
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
        elif trace_level == "summary":
            ensemble_trace = AdaptiveSamplingEnsembleTrace(
                ensemble_index=ensemble_index,
                initial_negative_sites=list(initial_negative_sites),
                final_top_sites=final_top_sites,
            )
        else:
            ensemble_trace = AdaptiveSamplingEnsembleTrace(
                ensemble_index=ensemble_index,
                initial_negative_sites=list(initial_negative_sites),
                iterations=iteration_traces,
                final_prediction_probabilities=pred_df,
                final_decision_values=final_decision_values,
                final_top_sites=final_top_sites,
            )

    return positive_series, ensemble_trace


def normalize_probabilities(values: np.ndarray) -> np.ndarray | None:
    total = float(np.nansum(values))
    if total <= 0.0 or not np.isfinite(total):
        return None
    return values / total


def transform_resampling_probabilities(
    values: np.ndarray,
    *,
    svm_mode: PredictionSvmMode,
) -> np.ndarray:
    """Adjust class resampling weights before normalization.

    In the native non-replayed path, scikit-learn's probability calibration can
    be slightly peakier than the e1071 path used by PhosR. Flattening the
    resampling weights a little reduces learner-side jitter in the adaptive
    sampling loop without changing candidate selection or the replay override
    seam.
    """

    weights = np.asarray(values, dtype=float)
    if svm_mode == "default":
        return np.power(weights, 0.8)
    return weights


def coerce_sampling_trace(
    sampling_trace: PredictionSamplingTrace | str | Path | None,
) -> PredictionSamplingTrace | None:
    if sampling_trace is None:
        return None
    if isinstance(sampling_trace, PredictionSamplingTrace):
        return sampling_trace
    return PredictionSamplingTrace.from_trace_directory(sampling_trace)


def resolve_sampled_site_positions(
    *,
    available_sites: pd.Index,
    sampled_sites: list[str],
    expected_size: int,
    context: str,
) -> np.ndarray:
    sampled_site_list = validate_override_sites(
        available_sites=available_sites,
        sampled_sites=sampled_sites,
        expected_size=expected_size,
        context=context,
    )
    position_lookup: dict[str, int] = {}
    for position, site in enumerate(available_sites.astype(str).tolist()):
        position_lookup.setdefault(site, position)
    return np.asarray([position_lookup[site] for site in sampled_site_list], dtype=int)


def validate_override_sites(
    *,
    available_sites: pd.Index,
    sampled_sites: list[str],
    expected_size: int,
    context: str,
) -> list[str]:
    sampled_site_list = [str(site) for site in sampled_sites]
    if len(sampled_site_list) != expected_size:
        msg = (
            f"Sampling override for {context} has {len(sampled_site_list)} rows; "
            f"expected {expected_size}"
        )
        raise ValueError(msg)
    available_site_set = set(available_sites.astype(str).tolist())
    invalid_sites = [
        site for site in sampled_site_list if site not in available_site_set
    ]
    if invalid_sites:
        unique_invalid_sites = sorted(set(invalid_sites))
        msg = (
            f"Sampling override for {context} contains sites outside the "
            f"available training rows: {', '.join(unique_invalid_sites)}"
        )
        raise ValueError(msg)
    return sampled_site_list


__all__ = [
    "coerce_sampling_trace",
    "make_prediction_random_generators",
    "multi_ada_sampling",
]
