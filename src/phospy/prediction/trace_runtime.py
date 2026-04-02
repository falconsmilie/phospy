from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

import pandas as pd

from ..types import PredictionTraceFormat, PredictionTraceLevel

if TYPE_CHECKING:
    from ..validation.prediction import PredictionRequest

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
        except (
            ImportError,
            ModuleNotFoundError,
            ValueError,
        ) as error:  # pragma: no cover - engine availability varies
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
        except (
            ImportError,
            ModuleNotFoundError,
            ValueError,
        ) as error:  # pragma: no cover - engine availability varies
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
    sampling_trace: object | None
    trace_sink: TraceSink | None
    owns_trace_sink: bool = False

    def close_owned_trace_sink(self) -> None:
        if self.owns_trace_sink and self.trace_sink is not None:
            self.trace_sink.close()


@dataclass(slots=True)
class PredictionRuntimeSession:
    runtime_request: PredictionRequest
    execution_context: PredictionExecutionContext

    def __enter__(self) -> PredictionRuntimeSession:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> None:
        self.execution_context.close_owned_trace_sink()


def build_prediction_execution_context(
    *,
    sampling_trace: object | str | Path | None,
    trace_level: PredictionTraceLevel,
    trace_sink: TraceSink | str | Path | None,
    trace_sink_format: PredictionTraceFormat,
    trace_sink_max_buffer_rows: int = 1000,
) -> PredictionExecutionContext:
    """Resolve runtime-only prediction execution resources.

    Sampling trace replay input is normalized first so invalid replay sources do
    not leave behind an owned temporary trace sink.
    """

    from .sampling_runtime import coerce_sampling_trace

    resolved_sampling_trace = coerce_sampling_trace(sampling_trace)
    owns_trace_sink = False
    resolved_trace_sink = None
    if trace_level == "full":
        owns_trace_sink = not isinstance(trace_sink, TraceSink)
        resolved_trace_sink = create_trace_sink(
            trace_sink,
            fmt=trace_sink_format,
            max_buffer_rows=trace_sink_max_buffer_rows,
        )
    return PredictionExecutionContext(
        sampling_trace=resolved_sampling_trace,
        trace_sink=resolved_trace_sink,
        owns_trace_sink=owns_trace_sink,
    )


def build_prediction_runtime_session(
    request: PredictionRequest,
    *,
    trace_sink_max_buffer_rows: int = 1000,
) -> PredictionRuntimeSession:
    execution_context = build_prediction_execution_context(
        sampling_trace=request.sampling_trace,
        trace_level=request.trace_level,
        trace_sink=request.trace_sink,
        trace_sink_format=request.trace_sink_format,
        trace_sink_max_buffer_rows=trace_sink_max_buffer_rows,
    )
    runtime_request = request.model_copy(
        update={
            "sampling_trace": execution_context.sampling_trace,
            "trace_sink": execution_context.trace_sink,
        }
    )
    return PredictionRuntimeSession(
        runtime_request=runtime_request,
        execution_context=execution_context,
    )


__all__ = [
    "DirectoryTraceSink",
    "PredictionExecutionContext",
    "PredictionRuntimeSession",
    "TRACE_TABLE_NAMES",
    "TraceSink",
    "build_prediction_execution_context",
    "build_prediction_runtime_session",
    "create_trace_sink",
]
