from __future__ import annotations

import numpy as np
import pandas as pd

from ..internal.types import PREDICTION_TRACE_LEVEL_FULL
from .results import KinasePredictionResult
from .trace_runtime import TRACE_TABLE_NAMES


def _empty_trace_tables() -> dict[str, pd.DataFrame]:
    return {name: pd.DataFrame() for name in TRACE_TABLE_NAMES}


def _flatten_candidate_rows(result: KinasePredictionResult) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for kinase, trace in (result.debug_traces or {}).items():
        for rank, site in enumerate(trace.candidate_substrates, start=1):
            rows.append({"kinase": kinase, "candidate_rank": rank, "site": site})
    return rows


def _flatten_negative_pool_rows(
    result: KinasePredictionResult,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for kinase, trace in (result.debug_traces or {}).items():
        for pool_index, site in enumerate(trace.negative_pool_sites, start=1):
            rows.append({"kinase": kinase, "pool_index": pool_index, "site": site})
    return rows


def _flatten_initial_rows(result: KinasePredictionResult) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for kinase, trace in (result.debug_traces or {}).items():
        for ensemble_trace in trace.ensemble_traces:
            for draw, site in enumerate(ensemble_trace.initial_negative_sites, start=1):
                rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": ensemble_trace.ensemble_index,
                        "draw": draw,
                        "site": site,
                    }
                )
    return rows


def _probability_column(frame: pd.DataFrame, class_label: str) -> np.ndarray:
    if class_label in frame.columns:
        return frame.loc[:, class_label].to_numpy(dtype=float, copy=False)
    return np.full(len(frame.index), np.nan, dtype=float)


def _flatten_iteration_rows(
    result: KinasePredictionResult,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    label_rows: list[dict[str, object]] = []
    probability_rows: list[dict[str, object]] = []
    probability_parameter_rows: list[dict[str, object]] = []
    decision_rows: list[dict[str, object]] = []
    for kinase, trace in (result.debug_traces or {}).items():
        for ensemble_trace in trace.ensemble_traces:
            for iteration_trace in ensemble_trace.iterations:
                labels = iteration_trace.labels.to_numpy(dtype=int, copy=False)
                probs = iteration_trace.probabilities
                sites = probs.index.tolist()
                class_1_probs = _probability_column(probs, "1")
                class_2_probs = _probability_column(probs, "2")
                decision_values = iteration_trace.decision_values.to_numpy(
                    dtype=float,
                    copy=False,
                )
                for (
                    site,
                    label_value,
                    class_1_prob,
                    class_2_prob,
                    decision_value,
                ) in zip(
                    sites,
                    labels,
                    class_1_probs,
                    class_2_probs,
                    decision_values,
                    strict=True,
                ):
                    normalized_label = int(label_value)
                    label_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "site": site,
                            "label": normalized_label,
                        }
                    )
                    probability_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "site": site,
                            "label": normalized_label,
                            "prob_class_1": float(class_1_prob),
                            "prob_class_2": float(class_2_prob),
                        }
                    )
                    decision_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "site": site,
                            "label": normalized_label,
                            "decision_value_class_1": float(decision_value),
                        }
                    )
                if iteration_trace.probability_parameters is not None:
                    for row in iteration_trace.probability_parameters.itertuples(
                        index=False
                    ):
                        probability_parameter_rows.append(
                            {
                                "kinase": kinase,
                                "ensemble": ensemble_trace.ensemble_index,
                                "iteration": iteration_trace.iteration_index,
                                "class_pair": str(row.class_pair),
                                "probA": float(row.probA),
                                "probB": float(row.probB),
                            }
                        )
    return label_rows, probability_rows, probability_parameter_rows, decision_rows


def _flatten_weight_rows(result: KinasePredictionResult) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for kinase, trace in (result.debug_traces or {}).items():
        for ensemble_trace in trace.ensemble_traces:
            for iteration_trace in ensemble_trace.iterations:
                if iteration_trace.positive_weights is not None:
                    positive_weights = iteration_trace.positive_weights
                    for site, weight in zip(
                        positive_weights.index.tolist(),
                        positive_weights.to_numpy(dtype=float, copy=False),
                        strict=True,
                    ):
                        rows.append(
                            {
                                "kinase": kinase,
                                "ensemble": ensemble_trace.ensemble_index,
                                "iteration": iteration_trace.iteration_index,
                                "class_label": 1,
                                "site": site,
                                "normalized_weight": float(weight),
                            }
                        )
                if iteration_trace.negative_weights is not None:
                    negative_weights = iteration_trace.negative_weights
                    for site, weight in zip(
                        negative_weights.index.tolist(),
                        negative_weights.to_numpy(dtype=float, copy=False),
                        strict=True,
                    ):
                        rows.append(
                            {
                                "kinase": kinase,
                                "ensemble": ensemble_trace.ensemble_index,
                                "iteration": iteration_trace.iteration_index,
                                "class_label": 2,
                                "site": site,
                                "normalized_weight": float(weight),
                            }
                        )
    return rows


def _flatten_sample_rows(result: KinasePredictionResult) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for kinase, trace in (result.debug_traces or {}).items():
        for ensemble_trace in trace.ensemble_traces:
            for iteration_trace in ensemble_trace.iterations:
                for draw, site in enumerate(
                    iteration_trace.sampled_positive_sites, start=1
                ):
                    rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "class_label": 1,
                            "draw": draw,
                            "site": site,
                        }
                    )
                for draw, site in enumerate(
                    iteration_trace.sampled_negative_sites, start=1
                ):
                    rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "class_label": 2,
                            "draw": draw,
                            "site": site,
                        }
                    )
    return rows


def _flatten_final_prediction_rows(
    result: KinasePredictionResult,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    final_prediction_rows: list[dict[str, object]] = []
    final_decision_rows: list[dict[str, object]] = []
    final_top_rows: list[dict[str, object]] = []
    for kinase, trace in (result.debug_traces or {}).items():
        for ensemble_trace in trace.ensemble_traces:
            final_probs = ensemble_trace.final_prediction_probabilities
            final_decisions = ensemble_trace.final_decision_values
            prob_class_1_by_site: dict[object, float] = {}
            if final_probs is not None:
                sites = final_probs.index.tolist()
                class_1_probs = _probability_column(final_probs, "1")
                class_2_probs = _probability_column(final_probs, "2")
                decision_values = (
                    final_decisions.to_numpy(dtype=float, copy=False)
                    if final_decisions is not None
                    else np.empty(0, dtype=float)
                )
                prob_class_1_by_site = {
                    site: float(probability)
                    for site, probability in zip(sites, class_1_probs, strict=True)
                }
                for position, site in enumerate(sites):
                    final_prediction_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "site": site,
                            "prob_class_1": float(class_1_probs[position]),
                            "prob_class_2": float(class_2_probs[position]),
                        }
                    )
                    if final_decisions is not None:
                        final_decision_rows.append(
                            {
                                "kinase": kinase,
                                "ensemble": ensemble_trace.ensemble_index,
                                "site": site,
                                "decision_value_class_1": float(
                                    decision_values[position]
                                ),
                            }
                        )
            for rank, site in enumerate(ensemble_trace.final_top_sites, start=1):
                final_top_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": ensemble_trace.ensemble_index,
                        "rank": rank,
                        "site": site,
                        "prob_class_1": prob_class_1_by_site.get(site, float("nan")),
                    }
                )
    return final_prediction_rows, final_decision_rows, final_top_rows


def prediction_debug_trace_tables(
    result: KinasePredictionResult,
) -> dict[str, pd.DataFrame]:
    if (
        result.trace_level == PREDICTION_TRACE_LEVEL_FULL
        and result.trace_sink is not None
    ):
        return result.trace_sink.read_all_tables()

    tables = _empty_trace_tables()
    tables["trace_selected_candidates"] = pd.DataFrame(_flatten_candidate_rows(result))
    tables["trace_negative_pool"] = pd.DataFrame(_flatten_negative_pool_rows(result))
    tables["trace_initial_negatives"] = pd.DataFrame(_flatten_initial_rows(result))
    (
        label_rows,
        probability_rows,
        probability_parameter_rows,
        decision_rows,
    ) = _flatten_iteration_rows(result)
    tables["trace_iteration_labels"] = pd.DataFrame(label_rows)
    tables["trace_iteration_probabilities"] = pd.DataFrame(probability_rows)
    tables["trace_iteration_probability_parameters"] = pd.DataFrame(
        probability_parameter_rows
    )
    tables["trace_iteration_decision_values"] = pd.DataFrame(decision_rows)
    tables["trace_iteration_resampling_weights"] = pd.DataFrame(
        _flatten_weight_rows(result)
    )
    tables["trace_iteration_samples"] = pd.DataFrame(_flatten_sample_rows(result))
    (
        final_prediction_rows,
        final_decision_rows,
        final_top_rows,
    ) = _flatten_final_prediction_rows(result)
    tables["trace_final_ensemble_predictions"] = pd.DataFrame(final_prediction_rows)
    tables["trace_final_ensemble_decision_values"] = pd.DataFrame(final_decision_rows)
    tables["trace_final_ensemble_top"] = pd.DataFrame(final_top_rows)
    return tables


__all__ = [
    "prediction_debug_trace_tables",
]
