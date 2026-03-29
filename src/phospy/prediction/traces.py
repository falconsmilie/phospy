from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..validation.errors import TableSchemaError
from .models import KinasePredictionResult, SamplingTraceOverrideEnsemble


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

        for ensemble_trace in trace.ensemble_traces:
            for draw, site in enumerate(ensemble_trace.initial_negative_sites, start=1):
                initial_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": ensemble_trace.ensemble_index,
                        "draw": draw,
                        "site": site,
                    }
                )

            for iteration_trace in ensemble_trace.iterations:
                labels = iteration_trace.labels
                probs = iteration_trace.probabilities
                for site in probs.index:
                    label_value = int(labels.loc[site])
                    label_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "site": site,
                            "label": label_value,
                        }
                    )
                    probability_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "site": site,
                            "label": label_value,
                            "prob_class_1": float(probs.loc[site, "1"])
                            if "1" in probs.columns
                            else float("nan"),
                            "prob_class_2": float(probs.loc[site, "2"])
                            if "2" in probs.columns
                            else float("nan"),
                        }
                    )
                    decision_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
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
                                "ensemble": ensemble_trace.ensemble_index,
                                "iteration": iteration_trace.iteration_index,
                                "class_pair": str(row["class_pair"]),
                                "probA": float(row["probA"]),
                                "probB": float(row["probB"]),
                            }
                        )
                if iteration_trace.positive_weights is not None:
                    for site, weight in iteration_trace.positive_weights.items():
                        weight_rows.append(
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
                    for site, weight in iteration_trace.negative_weights.items():
                        weight_rows.append(
                            {
                                "kinase": kinase,
                                "ensemble": ensemble_trace.ensemble_index,
                                "iteration": iteration_trace.iteration_index,
                                "class_label": 2,
                                "site": site,
                                "normalized_weight": float(weight),
                            }
                        )

                for draw, site in enumerate(
                    iteration_trace.sampled_positive_sites, start=1
                ):
                    sample_rows.append(
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
                    sample_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "class_label": 2,
                            "draw": draw,
                            "site": site,
                        }
                    )

            final_probs = ensemble_trace.final_prediction_probabilities
            for site in final_probs.index:
                final_prediction_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": ensemble_trace.ensemble_index,
                        "site": site,
                        "prob_class_1": float(final_probs.loc[site, "1"])
                        if "1" in final_probs.columns
                        else float("nan"),
                        "prob_class_2": float(final_probs.loc[site, "2"])
                        if "2" in final_probs.columns
                        else float("nan"),
                    }
                )
                final_decision_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": ensemble_trace.ensemble_index,
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
                        "ensemble": ensemble_trace.ensemble_index,
                        "rank": rank,
                        "site": site,
                        "prob_class_1": float(final_probs.loc[site, "1"])
                        if "1" in final_probs.columns
                        else float("nan"),
                    }
                )

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
