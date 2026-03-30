from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..validation.errors import TableSchemaError
from .models import (
    AdaptiveSamplingEnsembleTrace,
    AdaptiveSamplingIterationTrace,
    KinasePredictionDebugTrace,
    KinasePredictionResult,
    SamplingTraceOverrideEnsemble,
)


class PredictionSamplingTrace:
    def __init__(
        self, ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]]
    ) -> None:
        self.ensembles_by_kinase = ensembles_by_kinase

    @classmethod
    def from_trace_directory(cls, trace_dir: str | Path) -> PredictionSamplingTrace:
        path = Path(trace_dir)
        initial_path = path / "trace_initial_negatives.csv"
        samples_path = path / "trace_iteration_samples.csv"
        if not initial_path.exists() and not samples_path.exists():
            msg = (
                "sampling trace directory must contain trace_initial_negatives.csv "
                "and/or trace_iteration_samples.csv"
            )
            raise TableSchemaError(msg)

        ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]] = {}

        if initial_path.exists():
            initial_df = pd.read_csv(initial_path)
            required_initial_cols = {"kinase", "ensemble", "draw", "site"}
            if not required_initial_cols.issubset(initial_df.columns):
                missing = sorted(required_initial_cols.difference(initial_df.columns))
                msg = (
                    "trace_initial_negatives.csv is missing required columns: "
                    + ", ".join(missing)
                )
                raise TableSchemaError(msg)
            initial_df = initial_df.sort_values(
                ["kinase", "ensemble", "draw"], kind="mergesort"
            )
            for (kinase, ensemble), group in initial_df.groupby(
                ["kinase", "ensemble"], sort=False
            ):
                ensemble_map = ensembles_by_kinase.setdefault(str(kinase), {})
                ensemble_map[int(ensemble)] = SamplingTraceOverrideEnsemble(
                    initial_negative_sites=group.loc[:, "site"].astype(str).tolist(),
                    iteration_sample_sites={},
                )

        if samples_path.exists():
            samples_df = pd.read_csv(samples_path)
            required_sample_cols = {
                "kinase",
                "ensemble",
                "iteration",
                "class_label",
                "draw",
                "site",
            }
            if not required_sample_cols.issubset(samples_df.columns):
                missing = sorted(required_sample_cols.difference(samples_df.columns))
                msg = (
                    "trace_iteration_samples.csv is missing required columns: "
                    + ", ".join(missing)
                )
                raise TableSchemaError(msg)
            samples_df = samples_df.sort_values(
                ["kinase", "ensemble", "iteration", "class_label", "draw"],
                kind="mergesort",
            )
            for (kinase, ensemble, iteration, class_label), group in samples_df.groupby(
                ["kinase", "ensemble", "iteration", "class_label"], sort=False
            ):
                ensemble_map = ensembles_by_kinase.setdefault(str(kinase), {})
                ensemble_override = ensemble_map.setdefault(
                    int(ensemble),
                    SamplingTraceOverrideEnsemble(
                        initial_negative_sites=None,
                        iteration_sample_sites={},
                    ),
                )
                iteration_map = ensemble_override.iteration_sample_sites.setdefault(
                    int(iteration), {}
                )
                iteration_map[int(class_label)] = (
                    group.loc[:, "site"].astype(str).tolist()
                )

        return cls(ensembles_by_kinase=ensembles_by_kinase)

    def get_ensemble_override(
        self, kinase: str, ensemble_index: int
    ) -> SamplingTraceOverrideEnsemble | None:
        return self.ensembles_by_kinase.get(kinase, {}).get(ensemble_index)

    def subset_kinases(self, kinases: list[str] | set[str]) -> PredictionSamplingTrace:
        kinase_set = {str(kinase) for kinase in kinases}
        return PredictionSamplingTrace(
            ensembles_by_kinase={
                kinase: ensemble_map
                for kinase, ensemble_map in self.ensembles_by_kinase.items()
                if kinase in kinase_set
            }
        )


def _get_probability_value(
    probabilities: pd.DataFrame, site: str, column: str
) -> float:
    if column not in probabilities.columns:
        return float("nan")
    return float(probabilities.loc[site, column])


def _flatten_candidate_rows(
    kinase: str,
    trace: KinasePredictionDebugTrace,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    candidate_rows: list[dict[str, object]] = []
    negative_pool_rows: list[dict[str, object]] = []

    for rank, site in enumerate(trace.candidate_substrates, start=1):
        candidate_rows.append(
            {
                "kinase": kinase,
                "candidate_rank": rank,
                "site": site,
            }
        )

    for pool_index, site in enumerate(trace.negative_pool_sites, start=1):
        negative_pool_rows.append(
            {
                "kinase": kinase,
                "pool_index": pool_index,
                "site": site,
            }
        )

    return candidate_rows, negative_pool_rows


def _flatten_weight_rows(
    *,
    kinase: str,
    ensemble_index: int,
    iteration_index: int,
    class_label: int,
    weights: pd.Series | None,
) -> list[dict[str, object]]:
    if weights is None:
        return []

    rows: list[dict[str, object]] = []
    for site, weight in weights.items():
        rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "class_label": class_label,
                "site": site,
                "normalized_weight": float(weight),
            }
        )
    return rows


def _flatten_sample_rows(
    *,
    kinase: str,
    ensemble_index: int,
    iteration_index: int,
    class_label: int,
    sampled_sites: list[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for draw, site in enumerate(sampled_sites, start=1):
        rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "class_label": class_label,
                "draw": draw,
                "site": site,
            }
        )
    return rows


def _flatten_iteration_rows(
    kinase: str,
    ensemble_trace: AdaptiveSamplingEnsembleTrace,
    iteration_trace: AdaptiveSamplingIterationTrace,
) -> dict[str, list[dict[str, object]]]:
    label_rows: list[dict[str, object]] = []
    probability_rows: list[dict[str, object]] = []
    probability_parameter_rows: list[dict[str, object]] = []
    decision_rows: list[dict[str, object]] = []

    labels = iteration_trace.labels
    probabilities = iteration_trace.probabilities
    ensemble_index = ensemble_trace.ensemble_index
    iteration_index = iteration_trace.iteration_index

    for site in probabilities.index:
        label_value = int(labels.loc[site])
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
                "prob_class_1": _get_probability_value(probabilities, site, "1"),
                "prob_class_2": _get_probability_value(probabilities, site, "2"),
            }
        )
        decision_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "iteration": iteration_index,
                "site": site,
                "label": label_value,
                "decision_value_class_1": float(
                    iteration_trace.decision_values.loc[site]
                ),
            }
        )

    if iteration_trace.probability_parameters is not None:
        for _, row in iteration_trace.probability_parameters.iterrows():
            probability_parameter_rows.append(
                {
                    "kinase": kinase,
                    "ensemble": ensemble_index,
                    "iteration": iteration_index,
                    "class_pair": str(row["class_pair"]),
                    "probA": float(row["probA"]),
                    "probB": float(row["probB"]),
                }
            )

    return {
        "labels": label_rows,
        "probabilities": probability_rows,
        "probability_parameters": probability_parameter_rows,
        "decision_values": decision_rows,
        "weights": [
            *_flatten_weight_rows(
                kinase=kinase,
                ensemble_index=ensemble_index,
                iteration_index=iteration_index,
                class_label=1,
                weights=iteration_trace.positive_weights,
            ),
            *_flatten_weight_rows(
                kinase=kinase,
                ensemble_index=ensemble_index,
                iteration_index=iteration_index,
                class_label=2,
                weights=iteration_trace.negative_weights,
            ),
        ],
        "samples": [
            *_flatten_sample_rows(
                kinase=kinase,
                ensemble_index=ensemble_index,
                iteration_index=iteration_index,
                class_label=1,
                sampled_sites=iteration_trace.sampled_positive_sites,
            ),
            *_flatten_sample_rows(
                kinase=kinase,
                ensemble_index=ensemble_index,
                iteration_index=iteration_index,
                class_label=2,
                sampled_sites=iteration_trace.sampled_negative_sites,
            ),
        ],
    }


def _flatten_initial_negative_rows(
    kinase: str,
    ensemble_trace: AdaptiveSamplingEnsembleTrace,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
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


def _flatten_final_prediction_rows(
    kinase: str,
    ensemble_trace: AdaptiveSamplingEnsembleTrace,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    final_prediction_rows: list[dict[str, object]] = []
    final_decision_rows: list[dict[str, object]] = []
    final_top_rows: list[dict[str, object]] = []

    final_probabilities = ensemble_trace.final_prediction_probabilities
    ensemble_index = ensemble_trace.ensemble_index

    for site in final_probabilities.index:
        final_prediction_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "site": site,
                "prob_class_1": _get_probability_value(final_probabilities, site, "1"),
                "prob_class_2": _get_probability_value(final_probabilities, site, "2"),
            }
        )
        final_decision_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "site": site,
                "decision_value_class_1": float(
                    ensemble_trace.final_decision_values.loc[site]
                ),
            }
        )

    for rank, site in enumerate(ensemble_trace.final_top_sites, start=1):
        final_top_rows.append(
            {
                "kinase": kinase,
                "ensemble": ensemble_index,
                "rank": rank,
                "site": site,
                "prob_class_1": _get_probability_value(final_probabilities, site, "1"),
            }
        )

    return final_prediction_rows, final_decision_rows, final_top_rows


def prediction_debug_trace_tables(
    result: KinasePredictionResult,
) -> dict[str, pd.DataFrame]:
    initial_rows: list[dict[str, object]] = []
    label_rows: list[dict[str, object]] = []
    probability_rows: list[dict[str, object]] = []
    probability_parameter_rows: list[dict[str, object]] = []
    decision_rows: list[dict[str, object]] = []
    weight_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    final_prediction_rows: list[dict[str, object]] = []
    final_decision_rows: list[dict[str, object]] = []
    final_top_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    negative_pool_rows: list[dict[str, object]] = []

    debug_traces = result.debug_traces or {}
    for kinase, trace in debug_traces.items():
        candidate_trace_rows, negative_pool_trace_rows = _flatten_candidate_rows(
            kinase, trace
        )
        candidate_rows.extend(candidate_trace_rows)
        negative_pool_rows.extend(negative_pool_trace_rows)

        for ensemble_trace in trace.ensemble_traces:
            initial_rows.extend(_flatten_initial_negative_rows(kinase, ensemble_trace))

            for iteration_trace in ensemble_trace.iterations:
                iteration_rows = _flatten_iteration_rows(
                    kinase, ensemble_trace, iteration_trace
                )
                label_rows.extend(iteration_rows["labels"])
                probability_rows.extend(iteration_rows["probabilities"])
                probability_parameter_rows.extend(
                    iteration_rows["probability_parameters"]
                )
                decision_rows.extend(iteration_rows["decision_values"])
                weight_rows.extend(iteration_rows["weights"])
                sample_rows.extend(iteration_rows["samples"])

            (
                ensemble_prediction_rows,
                ensemble_decision_rows,
                ensemble_top_rows,
            ) = _flatten_final_prediction_rows(kinase, ensemble_trace)
            final_prediction_rows.extend(ensemble_prediction_rows)
            final_decision_rows.extend(ensemble_decision_rows)
            final_top_rows.extend(ensemble_top_rows)

    return {
        "trace_selected_candidates": pd.DataFrame(candidate_rows),
        "trace_negative_pool": pd.DataFrame(negative_pool_rows),
        "trace_initial_negatives": pd.DataFrame(initial_rows),
        "trace_iteration_labels": pd.DataFrame(label_rows),
        "trace_iteration_probabilities": pd.DataFrame(probability_rows),
        "trace_iteration_probability_parameters": pd.DataFrame(
            probability_parameter_rows
        ),
        "trace_iteration_decision_values": pd.DataFrame(decision_rows),
        "trace_iteration_resampling_weights": pd.DataFrame(weight_rows),
        "trace_iteration_samples": pd.DataFrame(sample_rows),
        "trace_final_ensemble_predictions": pd.DataFrame(final_prediction_rows),
        "trace_final_ensemble_decision_values": pd.DataFrame(final_decision_rows),
        "trace_final_ensemble_top": pd.DataFrame(final_top_rows),
    }


__all__ = ["PredictionSamplingTrace", "prediction_debug_trace_tables"]
