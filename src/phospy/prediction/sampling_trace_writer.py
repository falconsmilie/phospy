from __future__ import annotations

import numpy as np
import pandas as pd

from .trace_runtime import TraceSink


def probability_column(frame: pd.DataFrame, class_label: str) -> np.ndarray:
    if class_label in frame.columns:
        return frame.loc[:, class_label].to_numpy(dtype=float, copy=False)
    return np.full(len(frame.index), np.nan, dtype=float)


def write_iteration_trace_rows(
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
    sites = probabilities.index.tolist()
    label_values = labels.to_numpy(dtype=int, copy=False)
    class_1_probs = probability_column(probabilities, "1")
    class_2_probs = probability_column(probabilities, "2")
    decision_vector = decision_values.to_numpy(dtype=float, copy=False)

    label_rows: list[dict[str, object]] = []
    probability_rows: list[dict[str, object]] = []
    decision_rows: list[dict[str, object]] = []
    for site, label_value, class_1_prob, class_2_prob, decision_value in zip(
        sites,
        label_values,
        class_1_probs,
        class_2_probs,
        decision_vector,
        strict=True,
    ):
        normalized_label = int(label_value)
        label_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "site": site,
                "label": normalized_label,
            }
        )
        probability_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "site": site,
                "label": normalized_label,
                "prob_class_1": float(class_1_prob),
                "prob_class_2": float(class_2_prob),
            }
        )
        decision_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "site": site,
                "label": normalized_label,
                "decision_value_class_1": float(decision_value),
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
                    "class_pair": str(row.class_pair),
                    "probA": float(row.probA),
                    "probB": float(row.probB),
                }
                for row in probability_parameters.itertuples(index=False)
            ],
        )

    weight_rows: list[dict[str, object]] = []
    for class_label, weights in (
        (1, weights_by_class.get(1)),
        (2, weights_by_class.get(2)),
    ):
        if weights is None:
            continue
        for site, weight in zip(
            weights.index.tolist(),
            weights.to_numpy(dtype=float, copy=False),
            strict=True,
        ):
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


def write_final_trace_rows(
    *,
    trace_sink: TraceSink,
    kinase: str,
    ensemble_index: int,
    final_probabilities: pd.DataFrame,
    final_decision_values: pd.Series,
    final_top_sites: list[str],
) -> None:
    sites = final_probabilities.index.tolist()
    class_1_probs = probability_column(final_probabilities, "1")
    class_2_probs = probability_column(final_probabilities, "2")
    decision_vector = final_decision_values.to_numpy(dtype=float, copy=False)
    prediction_rows: list[dict[str, object]] = []
    decision_rows: list[dict[str, object]] = []
    prob_class_1_by_site: dict[object, float] = {}
    for site, class_1_prob, class_2_prob, decision_value in zip(
        sites,
        class_1_probs,
        class_2_probs,
        decision_vector,
        strict=True,
    ):
        normalized_prob_class_1 = float(class_1_prob)
        prob_class_1_by_site[site] = normalized_prob_class_1
        prediction_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "site": site,
                "prob_class_1": normalized_prob_class_1,
                "prob_class_2": float(class_2_prob),
            }
        )
        decision_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "site": site,
                "decision_value_class_1": float(decision_value),
            }
        )
    top_rows = [
        {
            "kinase": kinase,
            "ensemble": ensemble_index,
            "rank": rank,
            "site": site,
            "prob_class_1": prob_class_1_by_site.get(site, float("nan")),
        }
        for rank, site in enumerate(final_top_sites, start=1)
    ]
    trace_sink.write_rows("trace_final_ensemble_predictions", prediction_rows)
    trace_sink.write_rows("trace_final_ensemble_decision_values", decision_rows)
    trace_sink.write_rows("trace_final_ensemble_top", top_rows)


__all__ = [
    "probability_column",
    "write_final_trace_rows",
    "write_iteration_trace_rows",
]
