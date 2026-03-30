from __future__ import annotations

from dataclasses import dataclass, field
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


@dataclass(slots=True)
class KinasePredictionResult:
    pred_matrix: pd.DataFrame
    substrate_list: dict[str, list[str]]
    debug_traces: dict[str, KinasePredictionDebugTrace] | None = None
    trace_level: PredictionTraceLevel = "none"
    trace_sink: TraceSink | None = None

    def close(self) -> None:
        if self.trace_sink is not None:
            self.trace_sink.close()

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
    "SamplingTraceOverrideEnsemble",
]
