from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd

from ..errors import PredictionConfigurationError
from ..internal.types import (
    PREDICTION_TRACE_LEVEL_NONE,
    PredictionTraceLevel,
)
from ..validation.requests.prediction import PredictionRequest
from .contracts import EnsemblePredictorContract
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
    """Validated per-kinase score batch passed to prediction aggregation."""

    kinase: str
    score_values: np.ndarray
    score_index: pd.Index

    @property
    def scores(self) -> pd.Series:
        """Compatibility view for callers that still expect a pandas Series."""

        return pd.Series(
            self.score_values,
            index=self.score_index,
            dtype=float,
        )


@dataclass(slots=True)
class PredictionTraceState:
    """Mutable trace bookkeeping shared across kinase prediction execution."""

    trace_level: PredictionTraceLevel
    trace_sink: TraceSink | None
    traced_kinases: set[str]
    debug_traces: dict[str, KinasePredictionDebugTrace] | None


@dataclass(frozen=True, slots=True)
class PredictionSamplingSession:
    """Sampling policy and RNG source resolved for one prediction request."""

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


@dataclass(frozen=True, slots=True)
class IndexedFeatureMatrix:
    feature_mat: pd.DataFrame
    feature_values: np.ndarray
    feature_index: pd.Index
    feature_columns: pd.Index
    site_positions: Mapping[object, int]
    all_positions: np.ndarray

    @staticmethod
    def cache_signature_for(
        feature_mat: pd.DataFrame,
        *,
        feature_values: np.ndarray,
    ) -> tuple[int, int, int, tuple[int, int], int, tuple[int, ...]]:
        """Return a cheap cache signature for indexed feature reuse.

        The cached indexed matrix is safe to reuse only while the caller keeps
        the same frame object, row and column labels, shape, and underlying
        numeric buffer layout. Structural mutation or buffer replacement forces
        a rebuild, while value-only updates can reuse the cached numeric view.
        """

        array_interface = cast(
            tuple[int, bool],
            feature_values.__array_interface__["data"],
        )
        return (
            id(feature_mat),
            id(feature_mat.index),
            id(feature_mat.columns),
            feature_values.shape,
            int(array_interface[0]),
            tuple(feature_values.strides),
        )

    @classmethod
    def from_feature_matrix(
        cls,
        feature_mat: pd.DataFrame,
        *,
        feature_values: np.ndarray | None = None,
    ) -> IndexedFeatureMatrix:
        resolved_feature_values = (
            feature_mat.to_numpy(copy=False)
            if feature_values is None
            else np.asarray(feature_values)
        )
        return cls(
            feature_mat=feature_mat,
            feature_values=resolved_feature_values,
            feature_index=feature_mat.index,
            feature_columns=feature_mat.columns,
            site_positions={site: i for i, site in enumerate(feature_mat.index)},
            all_positions=np.arange(len(feature_mat.index), dtype=int),
        )

    def prepare_kinase(self, *, substrates: list[str]) -> PreparedKinaseTrainingData:
        positive_positions = np.fromiter(
            (self.site_positions[site] for site in substrates),
            dtype=int,
            count=len(substrates),
        )
        negative_mask = np.ones(len(self.all_positions), dtype=bool)
        negative_mask[positive_positions] = False
        negative_positions = self.all_positions[negative_mask]
        return PreparedKinaseTrainingData(
            positive_positions=positive_positions,
            negative_pool_positions=negative_positions,
            positive_sites=self.feature_index.take(positive_positions),
            negative_pool_sites=self.feature_index.take(negative_positions),
            positive_values=self.feature_values[positive_positions, :],
            labels=np.concatenate(
                [
                    np.repeat(1, len(positive_positions)),
                    np.repeat(2, len(positive_positions)),
                ]
            ),
        )


@dataclass(frozen=True, slots=True)
class PreparedKinaseTrainingData:
    positive_positions: np.ndarray
    negative_pool_positions: np.ndarray
    positive_sites: pd.Index
    negative_pool_sites: pd.Index
    positive_values: np.ndarray
    labels: np.ndarray


class NegativePoolSampler:
    """Resolve initial negative pools for one ensemble iteration."""

    @staticmethod
    def sample_initial_negative_positions(
        *,
        negative_pool_sites: pd.Index,
        negative_pool_positions: np.ndarray,
        site_positions: Mapping[object, int],
        positive_size: int,
        negative_sampling_rng: np.random.Generator,
        ensemble_override: object | None,
        kinase: str,
        ensemble_index: int,
    ) -> np.ndarray:
        override_sites = getattr(ensemble_override, "initial_negative_sites", None)
        if override_sites is not None:
            validated_sites = validate_override_sites(
                available_sites=negative_pool_sites,
                sampled_sites=override_sites,
                expected_size=positive_size,
                context=(
                    f"initial negatives for kinase={kinase}, ensemble={ensemble_index}"
                ),
            )
            return np.fromiter(
                (site_positions[site] for site in validated_sites),
                dtype=int,
                count=len(validated_sites),
            )

        return negative_sampling_rng.choice(
            negative_pool_positions,
            size=positive_size,
            replace=len(negative_pool_positions) < positive_size,
        )

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
        negative_indices = NegativePoolSampler.sample_initial_negative_positions(
            negative_pool_sites=negative_pool.index,
            negative_pool_positions=np.arange(len(negative_pool.index), dtype=int),
            site_positions={site: i for i, site in enumerate(negative_pool.index)},
            positive_size=positive_size,
            negative_sampling_rng=negative_sampling_rng,
            ensemble_override=ensemble_override,
            kinase=kinase,
            ensemble_index=ensemble_index,
        )
        return negative_pool.index.take(negative_indices).tolist()


class TraceRecorder:
    """Emit and collect structured prediction trace data."""

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
            if trace_level != PREDICTION_TRACE_LEVEL_NONE and debug_kinases is None
            else set(debug_kinases or [])
        )
        debug_traces: dict[str, KinasePredictionDebugTrace] | None = (
            {} if trace_level != PREDICTION_TRACE_LEVEL_NONE else None
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
            trace_state.trace_level != PREDICTION_TRACE_LEVEL_NONE
            and kinase in trace_state.traced_kinases
        )

    def start_kinase(
        self,
        *,
        trace_state: PredictionTraceState,
        kinase: str,
        substrates: list[str],
        negative_pool_sites: pd.Index,
    ) -> None:
        if not self.is_traced_kinase(trace_state=trace_state, kinase=kinase):
            return
        negative_pool_site_list = negative_pool_sites.tolist()
        if trace_state.debug_traces is not None:
            trace_state.debug_traces[kinase] = KinasePredictionDebugTrace(
                kinase=kinase,
                candidate_substrates=list(substrates),
                negative_pool_sites=negative_pool_site_list,
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
                    for pool_index, site in enumerate(negative_pool_site_list, start=1)
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


class EnsemblePredictor(EnsemblePredictorContract):
    """Default adaptive-ensemble kinase predictor implementation.

    This concrete implementation powers the built-in prediction engine. Advanced
    users may replace it via `EnsemblePredictorContract`.
    """

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
        self._indexed_feature_mat_cache: (
            tuple[
                tuple[int, int, int, tuple[int, int], int, tuple[int, ...]],
                IndexedFeatureMatrix,
            ]
            | None
        ) = None

    def _get_indexed_feature_matrix(
        self,
        feature_mat: pd.DataFrame,
    ) -> IndexedFeatureMatrix:
        feature_values = feature_mat.to_numpy(copy=False)
        cache_signature = IndexedFeatureMatrix.cache_signature_for(
            feature_mat,
            feature_values=feature_values,
        )
        cache = self._indexed_feature_mat_cache
        if cache is not None and cache[0] == cache_signature:
            return cache[1]
        indexed_feature_mat = IndexedFeatureMatrix.from_feature_matrix(
            feature_mat,
            feature_values=feature_values,
        )
        self._indexed_feature_mat_cache = (cache_signature, indexed_feature_mat)
        return indexed_feature_mat

    def clear_cache(self) -> None:
        self._indexed_feature_mat_cache = None

    def _prepare_kinase_training_data(
        self,
        *,
        kinase: str,
        substrates: list[str],
        feature_mat: pd.DataFrame,
    ) -> tuple[IndexedFeatureMatrix, PreparedKinaseTrainingData]:
        indexed_feature_mat = self._get_indexed_feature_matrix(feature_mat)
        prepared_training_data = indexed_feature_mat.prepare_kinase(
            substrates=substrates
        )
        if len(prepared_training_data.negative_pool_positions) == 0:
            msg = f"No negative pool sites available to train predictor for {kinase}"
            raise PredictionConfigurationError(msg)
        return indexed_feature_mat, prepared_training_data

    @staticmethod
    def _resolve_ensemble_override(
        *,
        request: PredictionRequest,
        kinase: str,
        ensemble_index: int,
    ) -> object | None:
        if request.sampling_trace is None:
            return None
        return request.sampling_trace.get_ensemble_override(
            kinase=kinase,
            ensemble_index=ensemble_index,
        )

    def _sample_initial_negatives(
        self,
        *,
        indexed_feature_mat: IndexedFeatureMatrix,
        prepared_training_data: PreparedKinaseTrainingData,
        negative_sampling_rng: np.random.Generator,
        ensemble_override: object | None,
        kinase: str,
        ensemble_index: int,
    ) -> tuple[np.ndarray, list[str]]:
        sampled_negative_positions = (
            self.negative_pool_sampler.sample_initial_negative_positions(
                negative_pool_sites=prepared_training_data.negative_pool_sites,
                negative_pool_positions=prepared_training_data.negative_pool_positions,
                site_positions=indexed_feature_mat.site_positions,
                positive_size=len(prepared_training_data.positive_positions),
                negative_sampling_rng=negative_sampling_rng,
                ensemble_override=ensemble_override,
                kinase=kinase,
                ensemble_index=ensemble_index,
            )
        )
        negative_sites = indexed_feature_mat.feature_index.take(
            sampled_negative_positions
        ).tolist()
        return sampled_negative_positions, negative_sites

    @staticmethod
    def _build_training_inputs(
        *,
        indexed_feature_mat: IndexedFeatureMatrix,
        prepared_training_data: PreparedKinaseTrainingData,
        sampled_negative_positions: np.ndarray,
    ) -> tuple[np.ndarray, pd.Index]:
        train_values = np.concatenate(
            [
                prepared_training_data.positive_values,
                indexed_feature_mat.feature_values[sampled_negative_positions, :],
            ],
            axis=0,
        )
        train_index = prepared_training_data.positive_sites.append(
            indexed_feature_mat.feature_index.take(sampled_negative_positions)
        )
        return train_values, train_index

    def _run_ensemble_sampling(
        self,
        *,
        request: PredictionRequest,
        trace_state: PredictionTraceState,
        sampling_session: PredictionSamplingSession,
        prepared_training_data: PreparedKinaseTrainingData,
        indexed_feature_mat: IndexedFeatureMatrix,
        kinase: str,
        ensemble_index: int,
        negative_sites: list[str],
        ensemble_override: object | None,
        train_values: np.ndarray,
        train_index: pd.Index,
        resampling_rng: np.random.Generator,
        capture_trace: bool,
    ) -> tuple[np.ndarray, object | None]:
        resolved_trace_level: PredictionTraceLevel = (
            request.trace_level if capture_trace else PREDICTION_TRACE_LEVEL_NONE
        )
        resolved_trace_sink = trace_state.trace_sink if capture_trace else None
        return multi_ada_sampling(
            train_mat=None,
            test_mat=None,
            labels=prepared_training_data.labels,
            kernel=self.kernel,
            n_iterations=request.n_iterations,
            resampling_rng=resampling_rng,
            capture_trace=capture_trace,
            trace_level=resolved_trace_level,
            trace_sink=resolved_trace_sink,
            kinase=kinase,
            ensemble_index=ensemble_index,
            initial_negative_sites=negative_sites,
            debug_top_n=request.debug_top_n,
            svm_mode=request.svm_mode,
            sampling_policy=sampling_session.policy,
            sampling_override=ensemble_override,
            train_values=train_values,
            train_index=train_index,
            test_values=indexed_feature_mat.feature_values,
            test_index=indexed_feature_mat.feature_index,
            return_values=True,
        )

    def _prepare_kinase_prediction(
        self,
        *,
        kinase: str,
        substrates: list[str],
        feature_mat: pd.DataFrame,
        trace_state: PredictionTraceState,
        sampling_session: PredictionSamplingSession,
    ) -> tuple[
        np.random.Generator,
        np.random.Generator,
        IndexedFeatureMatrix,
        PreparedKinaseTrainingData,
        bool,
        np.ndarray,
    ]:
        negative_sampling_rng, resampling_rng = (
            sampling_session.random_source.generators_for_kinase(kinase=kinase)
        )
        indexed_feature_mat, prepared_training_data = (
            self._prepare_kinase_training_data(
                kinase=kinase,
                substrates=substrates,
                feature_mat=feature_mat,
            )
        )
        self.trace_recorder.start_kinase(
            trace_state=trace_state,
            kinase=kinase,
            substrates=substrates,
            negative_pool_sites=prepared_training_data.negative_pool_sites,
        )
        is_traced_kinase = self.trace_recorder.is_traced_kinase(
            trace_state=trace_state,
            kinase=kinase,
        )
        capture_ensemble_trace = (
            is_traced_kinase and trace_state.debug_traces is not None
        )
        kinase_scores = np.zeros(len(feature_mat.index), dtype=float)
        return (
            negative_sampling_rng,
            resampling_rng,
            indexed_feature_mat,
            prepared_training_data,
            capture_ensemble_trace,
            kinase_scores,
        )

    def _run_ensemble_iteration(
        self,
        *,
        request: PredictionRequest,
        trace_state: PredictionTraceState,
        sampling_session: PredictionSamplingSession,
        prepared_training_data: PreparedKinaseTrainingData,
        indexed_feature_mat: IndexedFeatureMatrix,
        negative_sampling_rng: np.random.Generator,
        resampling_rng: np.random.Generator,
        kinase: str,
        ensemble_index: int,
        capture_ensemble_trace: bool,
    ) -> tuple[np.ndarray, object | None]:
        ensemble_override = self._resolve_ensemble_override(
            request=request,
            kinase=kinase,
            ensemble_index=ensemble_index,
        )
        sampled_negative_positions, negative_sites = self._sample_initial_negatives(
            indexed_feature_mat=indexed_feature_mat,
            prepared_training_data=prepared_training_data,
            negative_sampling_rng=negative_sampling_rng,
            ensemble_override=ensemble_override,
            kinase=kinase,
            ensemble_index=ensemble_index,
        )
        train_values, train_index = self._build_training_inputs(
            indexed_feature_mat=indexed_feature_mat,
            prepared_training_data=prepared_training_data,
            sampled_negative_positions=sampled_negative_positions,
        )
        self.trace_recorder.record_initial_negatives(
            trace_state=trace_state,
            kinase=kinase,
            ensemble_index=ensemble_index,
            negative_sites=negative_sites,
        )
        return self._run_ensemble_sampling(
            request=request,
            trace_state=trace_state,
            sampling_session=sampling_session,
            prepared_training_data=prepared_training_data,
            indexed_feature_mat=indexed_feature_mat,
            kinase=kinase,
            ensemble_index=ensemble_index,
            negative_sites=negative_sites,
            ensemble_override=ensemble_override,
            train_values=train_values,
            train_index=train_index,
            resampling_rng=resampling_rng,
            capture_trace=capture_ensemble_trace,
        )

    @staticmethod
    def _aggregate_ensemble_scores(
        *,
        kinase_scores: np.ndarray,
        score_values: np.ndarray,
    ) -> None:
        kinase_scores += score_values

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
        (
            negative_sampling_rng,
            resampling_rng,
            indexed_feature_mat,
            prepared_training_data,
            capture_ensemble_trace,
            kinase_scores,
        ) = self._prepare_kinase_prediction(
            kinase=kinase,
            substrates=substrates,
            feature_mat=feature_mat,
            trace_state=trace_state,
            sampling_session=sampling_session,
        )

        for ensemble_idx in range(request.ensemble_size):
            ensemble_index = ensemble_idx + 1
            score_values, ensemble_trace = self._run_ensemble_iteration(
                request=request,
                trace_state=trace_state,
                sampling_session=sampling_session,
                prepared_training_data=prepared_training_data,
                indexed_feature_mat=indexed_feature_mat,
                negative_sampling_rng=negative_sampling_rng,
                resampling_rng=resampling_rng,
                kinase=kinase,
                ensemble_index=ensemble_index,
                capture_ensemble_trace=capture_ensemble_trace,
            )
            if capture_ensemble_trace:
                self.trace_recorder.record_ensemble_trace(
                    trace_state=trace_state,
                    kinase=kinase,
                    ensemble_trace=ensemble_trace,
                )
            self._aggregate_ensemble_scores(
                kinase_scores=kinase_scores,
                score_values=score_values,
            )

        self.trace_recorder.flush_kinase(trace_state=trace_state, kinase=kinase)
        return KinasePredictionBatch(
            kinase=kinase,
            score_values=kinase_scores,
            score_index=indexed_feature_mat.feature_index,
        )


__all__ = [
    "EnsemblePredictor",
    "KinasePredictionBatch",
    "NegativePoolSampler",
    "PredictionTraceState",
    "PredictionSamplingSession",
    "TraceRecorder",
]
