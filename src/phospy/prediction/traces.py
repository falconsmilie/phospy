from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from ..types import PredictionTraceFormat, PredictionTraceLevel
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

    def flush(self) -> None:
        return None

    def close(self) -> None:
        return None

    def __enter__(self) -> TraceSink:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> None:
        self.close()

    def read_all_tables(self) -> dict[str, pd.DataFrame]:
        self.flush()
        return {name: self.read_table(name) for name in TRACE_TABLE_NAMES}


class DirectoryTraceSink(TraceSink):
    def __init__(
        self,
        output_dir: str | Path,
        *,
        fmt: PredictionTraceFormat = "csv",
        max_buffer_rows: int = 1000,
        owned_temp_dir: TemporaryDirectory[str] | None = None,
    ) -> None:
        if max_buffer_rows < 1:
            msg = "max_buffer_rows must be greater than or equal to 1"
            raise ValueError(msg)
        self._owned_temp_dir = owned_temp_dir
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fmt = fmt
        self.max_buffer_rows = max_buffer_rows
        self._part_counters: dict[str, int] = {name: 0 for name in TRACE_TABLE_NAMES}
        self._buffered_rows: dict[str, list[dict[str, object]]] = {
            name: [] for name in TRACE_TABLE_NAMES
        }
        self._closed = False

    def write_rows(self, table_name: str, rows: list[dict[str, object]]) -> None:
        if table_name not in TRACE_TABLE_NAMES:
            msg = f"Unsupported trace table: {table_name}"
            raise ValueError(msg)
        if not rows:
            return
        self._buffered_rows[table_name].extend(rows)
        if len(self._buffered_rows[table_name]) >= self.max_buffer_rows:
            self._flush_table(table_name)

    def _flush_table(self, table_name: str) -> None:
        rows = self._buffered_rows[table_name]
        if not rows:
            return

        frame = pd.DataFrame(rows)
        rows.clear()
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

    def flush(self) -> None:
        for table_name in TRACE_TABLE_NAMES:
            self._flush_table(table_name)

    def close(self) -> None:
        if self._closed:
            return
        self.flush()
        if self._owned_temp_dir is not None:
            self._owned_temp_dir.cleanup()
            self._owned_temp_dir = None
        self._closed = True

    def read_table(self, table_name: str) -> pd.DataFrame:
        if table_name not in TRACE_TABLE_NAMES:
            msg = f"Unsupported trace table: {table_name}"
            raise ValueError(msg)
        self._flush_table(table_name)
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
    max_buffer_rows: int = 1000,
) -> TraceSink:
    if isinstance(trace_sink, TraceSink):
        return trace_sink
    if trace_sink is None:
        temp_dir = TemporaryDirectory(prefix="phospy_prediction_trace_")
        return DirectoryTraceSink(
            temp_dir.name,
            fmt=fmt,
            max_buffer_rows=max_buffer_rows,
            owned_temp_dir=temp_dir,
        )
    return DirectoryTraceSink(trace_sink, fmt=fmt, max_buffer_rows=max_buffer_rows)


@dataclass(frozen=True, slots=True)
class PredictionExecutionContext:
    sampling_trace: PredictionSamplingTrace | None
    trace_sink: TraceSink | None


def build_prediction_execution_context(
    *,
    sampling_trace: PredictionSamplingTrace | str | Path | None,
    trace_level: PredictionTraceLevel,
    trace_sink: TraceSink | str | Path | None,
    trace_sink_format: PredictionTraceFormat,
    trace_sink_max_buffer_rows: int = 1000,
) -> PredictionExecutionContext:
    """Resolve runtime-only prediction execution resources.

    Sampling trace replay input is normalized first so invalid replay sources do
    not leave behind an owned temporary trace sink.
    """

    from .sampling import coerce_sampling_trace

    resolved_sampling_trace = coerce_sampling_trace(sampling_trace)
    resolved_trace_sink = (
        create_trace_sink(
            trace_sink,
            fmt=trace_sink_format,
            max_buffer_rows=trace_sink_max_buffer_rows,
        )
        if trace_level == "full"
        else None
    )
    return PredictionExecutionContext(
        sampling_trace=resolved_sampling_trace,
        trace_sink=resolved_trace_sink,
    )


def _read_trace_table_from_directory(
    trace_dir: Path,
    table_name: str,
) -> tuple[pd.DataFrame | None, str | None]:
    csv_path = trace_dir / f"{table_name}.csv"
    if csv_path.exists():
        return pd.read_csv(csv_path), csv_path.name

    parquet_path = trace_dir / f"{table_name}.parquet"
    parquet_parts = sorted(trace_dir.glob(f"{table_name}.part-*.parquet"))
    if parquet_path.exists() or parquet_parts:
        parquet_sources = [parquet_path] if parquet_path.exists() else parquet_parts
        try:
            frame = pd.concat(
                [pd.read_parquet(source) for source in parquet_sources],
                ignore_index=True,
            )
        except Exception as error:  # pragma: no cover - engine availability varies
            msg = (
                "Unable to read parquet trace replay input. Install a supported "
                "parquet engine such as 'pyarrow', or provide CSV trace files."
            )
            raise RuntimeError(msg) from error
        label = (
            parquet_path.name
            if parquet_path.exists()
            else f"{table_name}.part-*.parquet"
        )
        return frame, label

    return None, None


class PredictionSamplingTrace:
    def __init__(
        self, ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]]
    ) -> None:
        self.ensembles_by_kinase = ensembles_by_kinase

    @classmethod
    def from_trace_directory(cls, trace_dir: str | Path) -> PredictionSamplingTrace:
        path = Path(trace_dir)
        initial_df, initial_label = _read_trace_table_from_directory(
            path, "trace_initial_negatives"
        )
        samples_df, samples_label = _read_trace_table_from_directory(
            path, "trace_iteration_samples"
        )
        if initial_df is None and samples_df is None:
            msg = (
                "sampling trace directory must contain trace_initial_negatives "
                "and/or trace_iteration_samples in CSV or parquet format"
            )
            raise TableSchemaError(msg)

        ensembles_by_kinase: dict[str, dict[int, SamplingTraceOverrideEnsemble]] = {}

        if initial_df is not None:
            required_initial_cols = {"kinase", "ensemble", "draw", "site"}
            if not required_initial_cols.issubset(initial_df.columns):
                missing = sorted(required_initial_cols.difference(initial_df.columns))
                msg = f"{initial_label} is missing required columns: " + ", ".join(
                    missing
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

        if samples_df is not None:
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
                msg = f"{samples_label} is missing required columns: " + ", ".join(
                    missing
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
    "PredictionExecutionContext",
    "PredictionSamplingTrace",
    "TRACE_TABLE_NAMES",
    "TraceSink",
    "build_prediction_execution_context",
    "create_trace_sink",
    "prediction_debug_trace_tables",
]
