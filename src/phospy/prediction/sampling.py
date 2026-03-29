from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..types import PredictionSvmMode
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
from .traces import PredictionSamplingTrace


def make_prediction_random_generators(
    rng: np.random.Generator,
) -> tuple[np.random.Generator, np.random.Generator]:
    """Create independent RNG streams for prediction sampling steps."""

    return (
        np.random.default_rng(int(rng.integers(0, 2**31 - 1))),
        np.random.default_rng(int(rng.integers(0, 2**31 - 1))),
    )


def multi_ada_sampling(
    train_mat: pd.DataFrame,
    test_mat: pd.DataFrame,
    labels: np.ndarray,
    kernel: str,
    n_iterations: int,
    resampling_rng: np.random.Generator,
    capture_trace: bool,
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
        model = make_svm(
            StandardScaler=StandardScaler,
            SVC=SVC,
            kernel=kernel,
            svm_mode=svm_mode,
        )
        model.fit(current_x, current_y)
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

        resampled_x: list[np.ndarray] = []
        resampled_y: list[np.ndarray] = []
        sampled_sites_by_class: dict[int, list[str]] = {}
        weights_by_class: dict[int, pd.Series | None] = {}
        for class_idx, class_label in enumerate(model.classes_):
            class_mask = base_y == class_label
            class_x = base_x[class_mask]
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
            sampled_sites = class_index[sampled_idx].tolist()
            sampled_sites_by_class[int(class_label)] = sampled_sites
            resampled_x.append(class_x[sampled_idx])
            resampled_y.append(np.repeat(class_label, class_x.shape[0]))

        if capture_trace:
            iteration_traces.append(
                AdaptiveSamplingIterationTrace(
                    iteration_index=iteration_index,
                    labels=label_series,
                    probabilities=prob_df,
                    probability_parameters=probability_parameters,
                    decision_values=decision_series,
                    positive_weights=weights_by_class.get(1),
                    negative_weights=weights_by_class.get(2),
                    sampled_positive_sites=sampled_sites_by_class.get(1, []),
                    sampled_negative_sites=sampled_sites_by_class.get(2, []),
                )
            )

        current_x = np.vstack(resampled_x)
        current_y = np.concatenate(resampled_y)

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
