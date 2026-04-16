from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..internal.pandas_copy import detached_frame_copy
from ..internal.types import PREDICTION_TRACE_LEVEL_NONE, PredictionTraceLevel

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


@dataclass(slots=True)
class PredMatResult:
    """Stable predMat contract with explicit ownership-aware access.

    Rows are phosphosite identifiers, columns are kinase identifiers, and each
    value is the final predicted score for that phosphosite-kinase pair. The
    wrapped DataFrame is the owned in-memory predMat produced by the package.
    This wrapper is intentionally a plain mutable container around that owned
    frame rather than a frozen value object.

    Ownership convention:
    - ``to_frame()`` returns a detached safe copy
    - ``to_owned_frame()`` returns cheap shared owned state
    - ``to_mutable_frame_unsafe()`` returns explicit mutable shared state
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

    def to_frame(self, *, copy: bool | None = None) -> pd.DataFrame:
        """Return a detached predMat frame.

        The deprecated ``copy`` parameter is still accepted for compatibility.
        ``copy=False`` routes to ``to_owned_frame()``.
        """

        if copy is False:
            return self.to_owned_frame()
        return detached_frame_copy(self.data_frame)

    def to_owned_frame(self) -> pd.DataFrame:
        """Return the cheap shared owned predMat frame (no copy)."""

        return self.data_frame

    def to_mutable_frame_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared predMat state.

        Warning: mutating this frame mutates the owning prediction result.
        """

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
    trace_level: PredictionTraceLevel = PREDICTION_TRACE_LEVEL_NONE
    trace_sink: TraceSink | None = None
    owns_trace_sink: bool = False
    pred_mat_result: PredMatResult = field(init=False)

    def __post_init__(self) -> None:
        self.pred_mat_result = PredMatResult(self.pred_matrix)

    def to_pred_matrix(self) -> pd.DataFrame:
        """Return a detached prediction matrix."""

        return self.pred_mat_result.to_frame()

    def to_owned_pred_matrix(self) -> pd.DataFrame:
        """Return the cheap shared owned prediction matrix (no copy)."""

        return self.pred_mat_result.to_owned_frame()

    def to_mutable_pred_matrix_unsafe(self) -> pd.DataFrame:
        """Return explicit mutable shared prediction matrix state.

        Warning: mutating this frame mutates the owning prediction result.
        """

        return self.pred_mat_result.to_mutable_frame_unsafe()

    def to_substrate_list(self) -> dict[str, list[str]]:
        """Return a detached substrate-list mapping."""

        return {
            kinase: list(substrates)
            for kinase, substrates in self.substrate_list.items()
        }

    def to_owned_substrate_list(self) -> dict[str, list[str]]:
        """Return the cheap shared owned substrate-list mapping (no copy)."""

        return self.substrate_list

    def to_mutable_substrate_list_unsafe(self) -> dict[str, list[str]]:
        """Return explicit mutable shared substrate-list state.

        Warning: mutating this mapping mutates the owning prediction result.
        """

        return self.substrate_list

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
