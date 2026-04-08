from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..types import PredictionTraceLevel

if TYPE_CHECKING:
    from .traces import TraceSink


@dataclass(slots=True)
class AdaptiveSamplingIterationTrace:
    iteration_index: int
    labels: pd.Series
    probabilities: pd.DataFrame
    probability_parameters: pd.DataFrame | None
    decision_values: pd.Series
    positive_weights: pd.Series | None
    negative_weights: pd.Series | None
    sampled_positive_sites: list[str]
    sampled_negative_sites: list[str]


@dataclass(slots=True)
class AdaptiveSamplingEnsembleTrace:
    ensemble_index: int
    initial_negative_sites: list[str]
    iterations: list[AdaptiveSamplingIterationTrace] = field(default_factory=list)
    final_prediction_probabilities: pd.DataFrame | None = None
    final_decision_values: pd.Series | None = None
    final_top_sites: list[str] = field(default_factory=list)


@dataclass(slots=True)
class KinasePredictionDebugTrace:
    kinase: str
    candidate_substrates: list[str]
    negative_pool_sites: list[str]
    ensemble_traces: list[AdaptiveSamplingEnsembleTrace]


@dataclass(frozen=True, slots=True)
class PredMatResult:
    """Stable predMat contract with explicit in-memory and CSV access.

    Rows are phosphosite identifiers, columns are kinase identifiers, and each
    value is the final predicted score for that phosphosite-kinase pair. The
    wrapped DataFrame is the owned in-memory predMat produced by the package.
    Use ``to_frame(copy=True)`` when you need a detached copy.
    """

    data_frame: pd.DataFrame

    def __post_init__(self) -> None:
        if self.data_frame.index.name is None:
            self.data_frame.index.name = "phosphosite"

    @property
    def phosphosite_ids(self) -> pd.Index:
        """Return the phosphosite identifiers stored on the predMat rows."""

        return self.data_frame.index

    @property
    def kinase_names(self) -> pd.Index:
        """Return the kinase identifiers stored on the predMat columns."""

        return self.data_frame.columns

    def to_frame(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the in-memory predMat as a pandas DataFrame.

        By default, this returns a detached deep copy so callers can mutate it
        without affecting the owning result object. Pass ``copy=False`` to work
        with the owned in-memory frame directly.
        """

        if copy:
            return self.data_frame.copy(deep=True)
        return self.data_frame

    def to_csv(
        self,
        path: str | Path,
        *,
        index_label: str = "phosphosite",
    ) -> Path:
        """Persist the predMat to a canonical CSV representation."""

        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        self.data_frame.to_csv(
            target,
            index=True,
            index_label=index_label,
            encoding="utf-8",
            float_format="%.17g",
            lineterminator="\n",
        )
        return target


@dataclass(slots=True)
class KinasePredictionResult:
    """Prediction outputs and optional traces for one prediction run.

    The canonical predMat contract is exposed through ``pred_mat_result``.
    The lower-level prediction matrix remains available as ``pred_matrix``.
    When ``owns_trace_sink`` is ``True``, call ``close()`` or use a context
    manager to release the trace sink deterministically.
    """

    pred_matrix: pd.DataFrame
    substrate_list: dict[str, list[str]]
    debug_traces: dict[str, KinasePredictionDebugTrace] | None = None
    trace_level: PredictionTraceLevel = "none"
    trace_sink: TraceSink | None = None
    owns_trace_sink: bool = False
    pred_mat_result: PredMatResult = field(init=False)

    def __post_init__(self) -> None:
        self.pred_mat_result = PredMatResult(self.pred_matrix)

    def close(self) -> None:
        if not self.owns_trace_sink or self.trace_sink is None:
            return
        self.trace_sink.close()
        self.owns_trace_sink = False

    def __enter__(self) -> KinasePredictionResult:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> None:
        self.close()


@dataclass(slots=True)
class SamplingTraceOverrideEnsemble:
    initial_negative_sites: list[str] | None
    iteration_sample_sites: dict[int, dict[int, list[str]]]

    def get_iteration_sample_sites(
        self, iteration_index: int, class_label: int
    ) -> list[str] | None:
        return self.iteration_sample_sites.get(iteration_index, {}).get(class_label)


__all__ = [
    "AdaptiveSamplingEnsembleTrace",
    "AdaptiveSamplingIterationTrace",
    "KinasePredictionDebugTrace",
    "KinasePredictionResult",
    "PredMatResult",
    "SamplingTraceOverrideEnsemble",
]
