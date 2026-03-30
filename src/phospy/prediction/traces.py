from __future__ import annotations

import tempfile
from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd

from ..types import PredictionTraceFormat
from ..validation.errors import TableSchemaError
from .models import KinasePredictionResult, SamplingTraceOverrideEnsemble

TRACE_TABLE_NAMES: tuple[str, ...] = (
    "trace_selected_candidates",
    "trace_negative_pool",
    "trace_initial_negatives",
    "trace_iteration_labels",
    "trace_iteration_probabilities",
    "trace_iteration_probability_parameters",
    "trace_iteration_decision_values",
    "trace_iteration_resampling_weights",
    "trace_iteration_samples",
    "trace_final_ensemble_predictions",
    "trace_final_ensemble_decision_values",
    "trace_final_ensemble_top",
)


class TraceSink(ABC):
    @abstractmethod
    def write_rows(self, table_name: str, rows: list[dict[str, object]]) -> None:
        raise NotImplementedError

    @abstractmethod
    def read_table(self, table_name: str) -> pd.DataFrame:
        raise NotImplementedError

    def read_all_tables(self) -> dict[str, pd.DataFrame]:
        return {name: self.read_table(name) for name in TRACE_TABLE_NAMES}


class DirectoryTraceSink(TraceSink):
    def __init__(
        self,
        output_dir: str | Path,
        *,
        fmt: PredictionTraceFormat = "csv",
    ) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fmt = fmt
        self._part_counters: dict[str, int] = {name: 0 for name in TRACE_TABLE_NAMES}

    def write_rows(self, table_name: str, rows: list[dict[str, object]]) -> None:
        if table_name not in TRACE_TABLE_NAMES:
            msg = f"Unsupported trace table: {table_name}"
            raise ValueError(msg)
        if not rows:
            return
        frame = pd.DataFrame(rows)
        if self.fmt == "csv":
            path = self.output_dir / f"{table_name}.csv"
            frame.to_csv(path, mode="a", header=not path.exists(), index=False)
            return

        part_index = self._part_counters[table_name]
        self._part_counters[table_name] = part_index + 1
        path = self.output_dir / f"{table_name}.part-{part_index:06d}.parquet"
        try:
            frame.to_parquet(path, index=False)
        except Exception as error:  # pragma: no cover - engine availability varies
            msg = (
                "Unable to write parquet trace output. Install a supported parquet "
                "engine such as 'pyarrow', or use trace_sink_format='csv'."
            )
            raise RuntimeError(msg) from error

    def read_table(self, table_name: str) -> pd.DataFrame:
        if table_name not in TRACE_TABLE_NAMES:
            msg = f"Unsupported trace table: {table_name}"
            raise ValueError(msg)
        if self.fmt == "csv":
            path = self.output_dir / f"{table_name}.csv"
            return pd.read_csv(path) if path.exists() else pd.DataFrame()

        parts = sorted(self.output_dir.glob(f"{table_name}.part-*.parquet"))
        if not parts:
            return pd.DataFrame()
        try:
            return pd.concat(
                [pd.read_parquet(part) for part in parts], ignore_index=True
            )
        except Exception as error:  # pragma: no cover - engine availability varies
            msg = (
                "Unable to read parquet trace output. Install a supported parquet "
                "engine such as 'pyarrow', or use trace_sink_format='csv'."
            )
            raise RuntimeError(msg) from error


def create_trace_sink(
    trace_sink: TraceSink | str | Path | None,
    *,
    fmt: PredictionTraceFormat,
) -> TraceSink:
    if isinstance(trace_sink, TraceSink):
        return trace_sink
    if trace_sink is None:
        temp_dir = Path(tempfile.mkdtemp(prefix="phospy_prediction_trace_"))
        return DirectoryTraceSink(temp_dir, fmt=fmt)
    return DirectoryTraceSink(trace_sink, fmt=fmt)


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
                labels = iteration_trace.labels
                probs = iteration_trace.probabilities
                for position, site in enumerate(probs.index.tolist()):
                    label_value = int(labels.iloc[position])
                    prob_row = probs.iloc[position]
                    decision_value = float(
                        iteration_trace.decision_values.iloc[position]
                    )
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
                            "prob_class_1": float(prob_row.get("1", float("nan"))),
                            "prob_class_2": float(prob_row.get("2", float("nan"))),
                        }
                    )
                    decision_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "site": site,
                            "label": label_value,
                            "decision_value_class_1": decision_value,
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
    return label_rows, probability_rows, probability_parameter_rows, decision_rows


def _flatten_weight_rows(result: KinasePredictionResult) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for kinase, trace in (result.debug_traces or {}).items():
        for ensemble_trace in trace.ensemble_traces:
            for iteration_trace in ensemble_trace.iterations:
                if iteration_trace.positive_weights is not None:
                    for site, weight in iteration_trace.positive_weights.items():
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
                    for site, weight in iteration_trace.negative_weights.items():
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
            if final_probs is not None:
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
                    if final_decisions is not None:
                        final_decision_rows.append(
                            {
                                "kinase": kinase,
                                "ensemble": ensemble_trace.ensemble_index,
                                "site": site,
                                "decision_value_class_1": float(
                                    final_decisions.loc[site]
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
                        "prob_class_1": (
                            float(final_probs.loc[site, "1"])
                            if final_probs is not None and "1" in final_probs.columns
                            else float("nan")
                        ),
                    }
                )
    return final_prediction_rows, final_decision_rows, final_top_rows


def prediction_debug_trace_tables(
    result: KinasePredictionResult,
) -> dict[str, pd.DataFrame]:
    if result.trace_level == "full" and result.trace_sink is not None:
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
    "DirectoryTraceSink",
    "PredictionSamplingTrace",
    "TRACE_TABLE_NAMES",
    "TraceSink",
    "create_trace_sink",
    "prediction_debug_trace_tables",
]
