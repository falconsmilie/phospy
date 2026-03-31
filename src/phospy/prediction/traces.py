from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
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

    def write_frame(self, table_name: str, frame: pd.DataFrame) -> None:
        if frame.empty:
            return
        self.write_rows(
            table_name,
            frame.to_dict(orient="records"),
        )

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
        owned_temp_dir: TemporaryDirectory[str] | None = None,
        flush_row_threshold: int = 1000,
    ) -> None:
        if flush_row_threshold < 1:
            msg = "flush_row_threshold must be at least 1"
            raise ValueError(msg)
        self._owned_temp_dir = owned_temp_dir
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fmt = fmt
        self.flush_row_threshold = flush_row_threshold
        self._part_counters: dict[str, int] = {name: 0 for name in TRACE_TABLE_NAMES}
        self._buffered_rows: dict[str, list[dict[str, object]]] = {
            name: [] for name in TRACE_TABLE_NAMES
        }
        self._closed = False

    @staticmethod
    def _validate_table_name(table_name: str) -> None:
        if table_name not in TRACE_TABLE_NAMES:
            msg = f"Unsupported trace table: {table_name}"
            raise ValueError(msg)

    def write_rows(self, table_name: str, rows: list[dict[str, object]]) -> None:
        self._validate_table_name(table_name)
        if not rows:
            return
        buffered_rows = self._buffered_rows[table_name]
        buffered_rows.extend(rows)
        if len(buffered_rows) >= self.flush_row_threshold:
            self._flush_table(table_name)

    def write_frame(self, table_name: str, frame: pd.DataFrame) -> None:
        self._validate_table_name(table_name)
        if frame.empty:
            return
        self._write_frame(table_name, frame)

    def _write_frame(self, table_name: str, frame: pd.DataFrame) -> None:
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

    def _flush_table(self, table_name: str) -> None:
        rows = self._buffered_rows[table_name]
        if not rows:
            return

        frame = pd.DataFrame.from_records(rows)
        rows.clear()
        self._write_frame(table_name, frame)

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
) -> TraceSink:
    if isinstance(trace_sink, TraceSink):
        return trace_sink
    if trace_sink is None:
        temp_dir = TemporaryDirectory(prefix="phospy_prediction_trace_")
        return DirectoryTraceSink(temp_dir.name, fmt=fmt, owned_temp_dir=temp_dir)
    return DirectoryTraceSink(trace_sink, fmt=fmt)


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


def _raise_missing_trace_columns(
    frame: pd.DataFrame,
    *,
    label: str,
    required_columns: set[str],
) -> None:
    if required_columns.issubset(frame.columns):
        return
    missing = sorted(required_columns.difference(frame.columns))
    msg = f"{label} is missing required columns: " + ", ".join(missing)
    raise TableSchemaError(msg)


def _format_trace_group(
    group_values: tuple[object, ...],
    group_columns: tuple[str, ...],
) -> str:
    return ", ".join(
        f"{column}={value}"
        for column, value in zip(group_columns, group_values, strict=True)
    )


def _normalize_trace_draw_column(frame: pd.DataFrame, *, label: str) -> pd.DataFrame:
    normalized = frame.copy()
    draw_values = pd.to_numeric(normalized.loc[:, "draw"], errors="coerce")
    valid_draws = draw_values.notna() & np.isfinite(draw_values)
    valid_draws &= np.equal(draw_values % 1, 0)
    valid_draws &= draw_values >= 1
    if not bool(valid_draws.all()):
        invalid_values = normalized.loc[~valid_draws, "draw"].astype(str).tolist()[:5]
        invalid_display = ", ".join(invalid_values)
        msg = (
            f"{label} draw values must be positive integers; invalid values: "
            f"{invalid_display}"
        )
        raise TableSchemaError(msg)
    normalized.loc[:, "draw"] = draw_values.astype(int)
    return normalized


def _validate_duplicate_trace_keys(
    frame: pd.DataFrame,
    *,
    label: str,
    key_columns: tuple[str, ...],
) -> None:
    duplicated = frame.duplicated(subset=list(key_columns), keep=False)
    if not bool(duplicated.any()):
        return
    first_duplicate = frame.loc[duplicated, list(key_columns)].iloc[0].tolist()
    key_display = _format_trace_group(tuple(first_duplicate), key_columns)
    joined_columns = ", ".join(key_columns)
    msg = f"{label} contains duplicate ({joined_columns}) entries: {key_display}"
    raise TableSchemaError(msg)


def _validate_contiguous_trace_draws(
    frame: pd.DataFrame,
    *,
    label: str,
    group_columns: tuple[str, ...],
) -> None:
    sorted_frame = frame.sort_values([*group_columns, "draw"], kind="mergesort")
    for group_values, group in sorted_frame.groupby(list(group_columns), sort=False):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        draws = group.loc[:, "draw"].astype(int).tolist()
        expected = list(range(1, len(draws) + 1))
        if draws == expected:
            continue
        key_display = _format_trace_group(group_values, group_columns)
        draw_display = ", ".join(str(draw) for draw in draws)
        msg = (
            f"{label} draw values must be contiguous within each "
            f"({', '.join(group_columns)}) group; {key_display} has draws: {draw_display}"
        )
        raise TableSchemaError(msg)


def _validate_initial_negative_replay_table(
    frame: pd.DataFrame,
    *,
    label: str,
) -> pd.DataFrame:
    required_columns = {"kinase", "ensemble", "draw", "site"}
    _raise_missing_trace_columns(frame, label=label, required_columns=required_columns)
    normalized = _normalize_trace_draw_column(frame, label=label)
    _validate_duplicate_trace_keys(
        normalized,
        label=label,
        key_columns=("kinase", "ensemble", "draw"),
    )
    _validate_contiguous_trace_draws(
        normalized,
        label=label,
        group_columns=("kinase", "ensemble"),
    )
    return normalized


def _validate_iteration_sample_replay_table(
    frame: pd.DataFrame,
    *,
    label: str,
) -> pd.DataFrame:
    required_columns = {
        "kinase",
        "ensemble",
        "iteration",
        "class_label",
        "draw",
        "site",
    }
    _raise_missing_trace_columns(frame, label=label, required_columns=required_columns)
    normalized = _normalize_trace_draw_column(frame, label=label)
    _validate_duplicate_trace_keys(
        normalized,
        label=label,
        key_columns=("kinase", "ensemble", "iteration", "class_label", "draw"),
    )
    _validate_contiguous_trace_draws(
        normalized,
        label=label,
        group_columns=("kinase", "ensemble", "iteration", "class_label"),
    )
    return normalized


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
            if initial_label is None:  # pragma: no cover - defensive guard
                initial_label = "trace_initial_negatives"
            initial_df = _validate_initial_negative_replay_table(
                initial_df,
                label=initial_label,
            )
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
            if samples_label is None:  # pragma: no cover - defensive guard
                samples_label = "trace_iteration_samples"
            samples_df = _validate_iteration_sample_replay_table(
                samples_df,
                label=samples_label,
            )
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
    "PredictionSamplingTrace",
    "TRACE_TABLE_NAMES",
    "TraceSink",
    "create_trace_sink",
    "prediction_debug_trace_tables",
]
