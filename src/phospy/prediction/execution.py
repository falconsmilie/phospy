from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..errors import PredictionConfigurationError
from ..internal.types import PredictionTraceLevel
from ..validation.requests import PredictionRequest
from .policies import (
    PredictionSamplingPolicy,
    PredictionSamplingRandomSource,
    resolve_prediction_sampling_policy,
)
from .results import KinasePredictionDebugTrace
from .sampling import (
    multi_ada_sampling,
    validate_override_sites,
)
from .traces import TraceSink


@dataclass(frozen=True, slots=True)
class KinasePredictionBatch:
    kinase: str
    scores: pd.Series


@dataclass(slots=True)
class PredictionTraceState:
    trace_level: PredictionTraceLevel
    trace_sink: TraceSink | None
    traced_kinases: set[str]
    debug_traces: dict[str, KinasePredictionDebugTrace] | None


@dataclass(frozen=True, slots=True)
class PredictionSamplingSession:
    policy: PredictionSamplingPolicy
    random_source: PredictionSamplingRandomSource

    @classmethod
    def from_request(cls, request: PredictionRequest) -> PredictionSamplingSession:
        policy = resolve_prediction_sampling_policy(request.svm_mode)
        return cls(
            policy=policy,
            random_source=PredictionSamplingRandomSource(
                policy=policy,
                random_state=request.random_state,
            ),
        )


class NegativePoolSampler:
    @staticmethod
    def sample_initial_negative_sites(
        *,
        negative_pool: pd.DataFrame,
        positive_size: int,
        negative_sampling_rng: np.random.Generator,
        ensemble_override: object | None,
        kinase: str,
        ensemble_index: int,
    ) -> list[str]:
        override_sites = getattr(ensemble_override, "initial_negative_sites", None)
        if override_sites is not None:
            return validate_override_sites(
                available_sites=negative_pool.index,
                sampled_sites=override_sites,
                expected_size=positive_size,
                context=(
                    f"initial negatives for kinase={kinase}, ensemble={ensemble_index}"
                ),
            )

        negative_indices = negative_sampling_rng.choice(
            negative_pool.index.to_numpy(),
            size=positive_size,
            replace=len(negative_pool) < positive_size,
        )
        return list(negative_indices.tolist())


class TraceRecorder:
    @staticmethod
    def create_state(
        *,
        substrate_list: dict[str, list[str]],
        trace_level: PredictionTraceLevel,
        debug_kinases: list[str] | None,
        trace_sink: TraceSink | None,
    ) -> PredictionTraceState:
        traced_kinases = (
            set(substrate_list)
            if trace_level != "none" and debug_kinases is None
            else set(debug_kinases or [])
        )
        debug_traces: dict[str, KinasePredictionDebugTrace] | None = (
            {} if trace_level != "none" else None
        )
        return PredictionTraceState(
            trace_level=trace_level,
            trace_sink=trace_sink,
            traced_kinases=traced_kinases,
            debug_traces=debug_traces,
        )

    def is_traced_kinase(
        self,
        *,
        trace_state: PredictionTraceState,
        kinase: str,
    ) -> bool:
        return (
            trace_state.trace_level != "none" and kinase in trace_state.traced_kinases
        )

    def start_kinase(
        self,
        *,
        trace_state: PredictionTraceState,
        kinase: str,
        substrates: list[str],
        negative_pool: pd.DataFrame,
    ) -> None:
        if not self.is_traced_kinase(trace_state=trace_state, kinase=kinase):
            return
        if trace_state.debug_traces is not None:
            trace_state.debug_traces[kinase] = KinasePredictionDebugTrace(
                kinase=kinase,
                candidate_substrates=list(substrates),
                negative_pool_sites=negative_pool.index.tolist(),
                ensemble_traces=[],
            )
        if trace_state.trace_sink is not None:
            trace_state.trace_sink.write_rows(
                "trace_selected_candidates",
                [
                    {
                        "kinase": kinase,
                        "candidate_rank": rank,
                        "site": site,
                    }
                    for rank, site in enumerate(substrates, start=1)
                ],
            )
            trace_state.trace_sink.write_rows(
                "trace_negative_pool",
                [
                    {
                        "kinase": kinase,
                        "pool_index": pool_index,
                        "site": site,
                    }
                    for pool_index, site in enumerate(
                        negative_pool.index.tolist(), start=1
                    )
                ],
            )

    def record_initial_negatives(
        self,
        *,
        trace_state: PredictionTraceState,
        kinase: str,
        ensemble_index: int,
        negative_sites: list[str],
    ) -> None:
        if not self.is_traced_kinase(trace_state=trace_state, kinase=kinase):
            return
        if trace_state.trace_sink is None:
            return
        trace_state.trace_sink.write_rows(
            "trace_initial_negatives",
            [
                {
                    "kinase": kinase,
                    "ensemble": ensemble_index,
                    "draw": draw,
                    "site": site,
                }
                for draw, site in enumerate(negative_sites, start=1)
            ],
        )

    @staticmethod
    def record_ensemble_trace(
        *,
        trace_state: PredictionTraceState,
        kinase: str,
        ensemble_trace: object | None,
    ) -> None:
        if ensemble_trace is None or trace_state.debug_traces is None:
            return
        if kinase not in trace_state.debug_traces:
            return
        trace_state.debug_traces[kinase].ensemble_traces.append(ensemble_trace)

    def flush_kinase(self, *, trace_state: PredictionTraceState, kinase: str) -> None:
        if (
            self.is_traced_kinase(trace_state=trace_state, kinase=kinase)
            and trace_state.trace_sink is not None
        ):
            trace_state.trace_sink.flush()

    @staticmethod
    def flush_final(*, trace_state: PredictionTraceState) -> None:
        if trace_state.trace_sink is not None:
            trace_state.trace_sink.flush()


class EnsemblePredictor:
    def __init__(
        self,
        *,
        kernel: str,
        negative_pool_sampler: NegativePoolSampler,
        trace_recorder: TraceRecorder,
    ) -> None:
        self.kernel = kernel
        self.negative_pool_sampler = negative_pool_sampler
        self.trace_recorder = trace_recorder

    def predict_kinase(
        self,
        *,
        kinase: str,
        substrates: list[str],
        feature_mat: pd.DataFrame,
        request: PredictionRequest,
        trace_state: PredictionTraceState,
        sampling_session: PredictionSamplingSession,
    ) -> KinasePredictionBatch:
        negative_sampling_rng, resampling_rng = (
            sampling_session.random_source.generators_for_kinase(kinase=kinase)
        )
        positive_train = feature_mat.loc[substrates, :]
        negative_pool = feature_mat.loc[
            ~feature_mat.index.isin(substrates),
            :,  # noqa: E203
        ]
        if negative_pool.empty:
            msg = f"No negative pool sites available to train predictor for {kinase}"
            raise PredictionConfigurationError(msg)

        self.trace_recorder.start_kinase(
            trace_state=trace_state,
            kinase=kinase,
            substrates=substrates,
            negative_pool=negative_pool,
        )

        kinase_scores = pd.Series(0.0, index=feature_mat.index.copy(), dtype=float)
        is_traced_kinase = self.trace_recorder.is_traced_kinase(
            trace_state=trace_state,
            kinase=kinase,
        )

        for ensemble_idx in range(request.ensemble_size):
            ensemble_index = ensemble_idx + 1
            ensemble_override = None
            if request.sampling_trace is not None:
                ensemble_override = request.sampling_trace.get_ensemble_override(
                    kinase=kinase,
                    ensemble_index=ensemble_index,
                )

            negative_sites = self.negative_pool_sampler.sample_initial_negative_sites(
                negative_pool=negative_pool,
                positive_size=len(positive_train),
                negative_sampling_rng=negative_sampling_rng,
                ensemble_override=ensemble_override,
                kinase=kinase,
                ensemble_index=ensemble_index,
            )
            negative_train = negative_pool.loc[negative_sites, :]
            train_mat = pd.concat([positive_train, negative_train], axis=0)
            labels = np.concatenate(
                [
                    np.repeat(1, len(positive_train)),
                    np.repeat(2, len(negative_train)),
                ]
            )

            self.trace_recorder.record_initial_negatives(
                trace_state=trace_state,
                kinase=kinase,
                ensemble_index=ensemble_index,
                negative_sites=negative_sites,
            )

            if is_traced_kinase and trace_state.debug_traces is not None:
                series, ensemble_trace = multi_ada_sampling(
                    train_mat=train_mat,
                    test_mat=feature_mat,
                    labels=labels,
                    kernel=self.kernel,
                    n_iterations=request.n_iterations,
                    resampling_rng=resampling_rng,
                    capture_trace=True,
                    trace_level=request.trace_level,
                    trace_sink=trace_state.trace_sink,
                    kinase=kinase,
                    ensemble_index=ensemble_index,
                    initial_negative_sites=negative_sites,
                    debug_top_n=request.debug_top_n,
                    svm_mode=request.svm_mode,
                    sampling_policy=sampling_session.policy,
                    sampling_override=ensemble_override,
                )
                self.trace_recorder.record_ensemble_trace(
                    trace_state=trace_state,
                    kinase=kinase,
                    ensemble_trace=ensemble_trace,
                )
            else:
                series, _ = multi_ada_sampling(
                    train_mat=train_mat,
                    test_mat=feature_mat,
                    labels=labels,
                    kernel=self.kernel,
                    n_iterations=request.n_iterations,
                    resampling_rng=resampling_rng,
                    capture_trace=False,
                    trace_level="none",
                    trace_sink=None,
                    kinase=kinase,
                    ensemble_index=ensemble_index,
                    initial_negative_sites=negative_sites,
                    debug_top_n=request.debug_top_n,
                    svm_mode=request.svm_mode,
                    sampling_policy=sampling_session.policy,
                    sampling_override=ensemble_override,
                )
            kinase_scores += series

        self.trace_recorder.flush_kinase(trace_state=trace_state, kinase=kinase)
        return KinasePredictionBatch(kinase=kinase, scores=kinase_scores)


__all__ = [
    "EnsemblePredictor",
    "KinasePredictionBatch",
    "NegativePoolSampler",
    "PredictionTraceState",
    "PredictionSamplingSession",
    "TraceRecorder",
]
