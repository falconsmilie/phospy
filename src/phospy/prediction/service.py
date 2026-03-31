from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..scoring import KinaseScoringResult
from ..types import PredictionSvmMode, PredictionTraceFormat, PredictionTraceLevel
from ..validation.errors import (
    InputCompatibilityError,
    PredictionConfigurationError,
)
from .models import KinasePredictionDebugTrace, KinasePredictionResult
from .sampling import (
    coerce_sampling_trace,
    make_prediction_random_generators,
    multi_ada_sampling,
    validate_override_sites,
)
from .traces import PredictionSamplingTrace, TraceSink, create_trace_sink
from .validation import (
    validate_positive_int,
    validate_probability_threshold,
    validate_svm_mode,
    validate_trace_format,
    validate_trace_level,
)


@dataclass(frozen=True, slots=True)
class PredictionRequest:
    combined_scores: pd.DataFrame
    ensemble_size: int
    top: int
    score_threshold: float
    inclusion: int
    n_iterations: int
    random_state: int | None
    debug_kinases: list[str] | None
    debug_top_n: int
    svm_mode: PredictionSvmMode
    sampling_trace: PredictionSamplingTrace | None
    trace_level: PredictionTraceLevel
    trace_sink: TraceSink | None


@dataclass(slots=True)
class PredictionTraceState:
    trace_level: PredictionTraceLevel
    trace_sink: TraceSink | None
    traced_kinases: set[str]
    debug_traces: dict[str, KinasePredictionDebugTrace] | None


@dataclass(frozen=True, slots=True)
class KinasePredictionBatch:
    kinase: str
    scores: pd.Series


class PredictionRequestFactory:
    def __init__(self, *, default_svm_mode: PredictionSvmMode) -> None:
        self.default_svm_mode = default_svm_mode

    def create(
        self,
        *,
        combined_scores: pd.DataFrame,
        ensemble_size: int,
        top: int,
        score_threshold: float,
        inclusion: int,
        n_iterations: int,
        random_state: int | None,
        capture_debug_trace: bool,
        debug_kinases: list[str] | None,
        debug_top_n: int,
        svm_mode: PredictionSvmMode | None,
        sampling_trace: PredictionSamplingTrace | str | Path | None,
        trace_level: PredictionTraceLevel | None,
        trace_sink: TraceSink | str | Path | None,
        trace_sink_format: PredictionTraceFormat,
    ) -> PredictionRequest:
        validate_positive_int(ensemble_size, name="ensemble_size")
        validate_positive_int(top, name="top")
        validate_probability_threshold(score_threshold, name="score_threshold")
        validate_positive_int(inclusion, name="inclusion")
        validate_positive_int(n_iterations, name="n_iterations")
        validate_positive_int(debug_top_n, name="debug_top_n")

        resolved_svm_mode = (
            self.default_svm_mode if svm_mode is None else validate_svm_mode(svm_mode)
        )
        resolved_trace_level = validate_trace_level(
            "summary"
            if trace_level is None and capture_debug_trace
            else trace_level or "none"
        )
        resolved_trace_format = validate_trace_format(trace_sink_format)

        return PredictionRequest(
            combined_scores=combined_scores,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            debug_kinases=debug_kinases,
            debug_top_n=debug_top_n,
            svm_mode=resolved_svm_mode,
            sampling_trace=coerce_sampling_trace(sampling_trace),
            trace_level=resolved_trace_level,
            trace_sink=(
                create_trace_sink(trace_sink, fmt=resolved_trace_format)
                if resolved_trace_level == "full"
                else None
            ),
        )


class CandidateSelector:
    def select(
        self,
        combined_scores: pd.DataFrame,
        *,
        top: int,
        score_threshold: float,
        inclusion: int,
    ) -> dict[str, list[str]]:
        return build_candidate_substrate_list(
            combined_scores,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
        )


class NegativePoolSampler:
    def sample_initial_negative_sites(
        self,
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
    def create_state(
        self,
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
            trace_state.trace_sink.write_frame(
                "trace_selected_candidates",
                pd.DataFrame(
                    {
                        "kinase": [kinase] * len(substrates),
                        "candidate_rank": np.arange(1, len(substrates) + 1),
                        "site": substrates,
                    }
                ),
            )
            negative_sites = negative_pool.index.tolist()
            trace_state.trace_sink.write_frame(
                "trace_negative_pool",
                pd.DataFrame(
                    {
                        "kinase": [kinase] * len(negative_sites),
                        "pool_index": np.arange(1, len(negative_sites) + 1),
                        "site": negative_sites,
                    }
                ),
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
        trace_state.trace_sink.write_frame(
            "trace_initial_negatives",
            pd.DataFrame(
                {
                    "kinase": [kinase] * len(negative_sites),
                    "ensemble": [ensemble_index] * len(negative_sites),
                    "draw": np.arange(1, len(negative_sites) + 1),
                    "site": negative_sites,
                }
            ),
        )

    def record_ensemble_trace(
        self,
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

    def flush_final(self, *, trace_state: PredictionTraceState) -> None:
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
        master_rng: np.random.Generator,
        trace_state: PredictionTraceState,
    ) -> KinasePredictionBatch:
        negative_sampling_rng, resampling_rng = make_prediction_random_generators(
            master_rng
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
                    sampling_override=ensemble_override,
                )
            kinase_scores += series

        self.trace_recorder.flush_kinase(trace_state=trace_state, kinase=kinase)
        return KinasePredictionBatch(kinase=kinase, scores=kinase_scores)


class PredictionAggregator:
    def empty_result(
        self,
        *,
        request: PredictionRequest,
    ) -> KinasePredictionResult:
        empty = pd.DataFrame(index=request.combined_scores.index.copy(), dtype=float)
        return KinasePredictionResult(
            pred_matrix=empty,
            substrate_list={},
            debug_traces={} if request.trace_level != "none" else None,
            trace_level=request.trace_level,
            trace_sink=request.trace_sink,
        )

    def initialize_prediction_matrix(
        self,
        *,
        feature_mat: pd.DataFrame,
        substrate_list: dict[str, list[str]],
    ) -> pd.DataFrame:
        return pd.DataFrame(
            0.0,
            index=feature_mat.index.copy(),
            columns=list(substrate_list),
        )

    def add_kinase_scores(
        self,
        *,
        pred_matrix: pd.DataFrame,
        batch: KinasePredictionBatch,
    ) -> None:
        pred_matrix.loc[:, batch.kinase] += batch.scores

    def finalize(
        self,
        *,
        pred_matrix: pd.DataFrame,
        substrate_list: dict[str, list[str]],
        request: PredictionRequest,
        trace_state: PredictionTraceState,
    ) -> KinasePredictionResult:
        pred_matrix /= float(request.ensemble_size)
        return KinasePredictionResult(
            pred_matrix=pred_matrix,
            substrate_list=substrate_list,
            debug_traces=trace_state.debug_traces,
            trace_level=request.trace_level,
            trace_sink=request.trace_sink,
        )


class KinasePredictor:
    """Predict kinase-substrate relationships from phosphosite score matrices.

    This is a narrow native Python port of PhosR's score-to-prediction seam.
    It mirrors the broad structure of ``kinaseSubstratePred()``: candidate
    substrates are selected from the combined score matrix, then an ensemble of
    adaptive SVM models is used to produce a kinase prediction matrix.

    Use ``svm_mode='r_parity'`` when you want settings that more closely match
    PhosR's e1071-based learner seam. The default mode preserves the standard
    scikit-learn behaviour.
    """

    def __init__(
        self,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
        *,
        request_factory: PredictionRequestFactory | None = None,
        candidate_selector: CandidateSelector | None = None,
        negative_pool_sampler: NegativePoolSampler | None = None,
        trace_recorder: TraceRecorder | None = None,
        prediction_aggregator: PredictionAggregator | None = None,
        ensemble_predictor: EnsemblePredictor | None = None,
    ) -> None:
        self.kernel = kernel
        self.svm_mode = validate_svm_mode(svm_mode)
        self.request_factory = request_factory or PredictionRequestFactory(
            default_svm_mode=self.svm_mode
        )
        self.candidate_selector = candidate_selector or CandidateSelector()
        self.negative_pool_sampler = negative_pool_sampler or NegativePoolSampler()
        self.trace_recorder = trace_recorder or TraceRecorder()
        self.prediction_aggregator = prediction_aggregator or PredictionAggregator()
        self.ensemble_predictor = ensemble_predictor or EnsemblePredictor(
            kernel=self.kernel,
            negative_pool_sampler=self.negative_pool_sampler,
            trace_recorder=self.trace_recorder,
        )

    def predict(
        self,
        combined_scores: pd.DataFrame,
        ensemble_size: int = 10,
        top: int = 50,
        score_threshold: float = 0.8,
        inclusion: int = 20,
        n_iterations: int = 5,
        random_state: int | None = None,
        capture_debug_trace: bool = False,
        debug_kinases: list[str] | None = None,
        debug_top_n: int = 10,
        svm_mode: PredictionSvmMode | None = None,
        sampling_trace: PredictionSamplingTrace | str | Path | None = None,
        trace_level: PredictionTraceLevel | None = None,
        trace_sink: TraceSink | str | Path | None = None,
        trace_sink_format: PredictionTraceFormat = "csv",
    ) -> KinasePredictionResult:
        request = self.request_factory.create(
            combined_scores=combined_scores,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            capture_debug_trace=capture_debug_trace,
            debug_kinases=debug_kinases,
            debug_top_n=debug_top_n,
            svm_mode=svm_mode,
            sampling_trace=sampling_trace,
            trace_level=trace_level,
            trace_sink=trace_sink,
            trace_sink_format=trace_sink_format,
        )
        substrate_list = self.candidate_selector.select(
            request.combined_scores,
            top=request.top,
            score_threshold=request.score_threshold,
            inclusion=request.inclusion,
        )
        if not substrate_list:
            return self.prediction_aggregator.empty_result(request=request)

        feature_mat = request.combined_scores.astype(float)
        pred_matrix = self.prediction_aggregator.initialize_prediction_matrix(
            feature_mat=feature_mat,
            substrate_list=substrate_list,
        )
        trace_state = self.trace_recorder.create_state(
            substrate_list=substrate_list,
            trace_level=request.trace_level,
            debug_kinases=request.debug_kinases,
            trace_sink=request.trace_sink,
        )
        master_rng = np.random.default_rng(request.random_state)

        for kinase, substrates in substrate_list.items():
            batch = self.ensemble_predictor.predict_kinase(
                kinase=kinase,
                substrates=substrates,
                feature_mat=feature_mat,
                request=request,
                master_rng=master_rng,
                trace_state=trace_state,
            )
            self.prediction_aggregator.add_kinase_scores(
                pred_matrix=pred_matrix,
                batch=batch,
            )

        self.trace_recorder.flush_final(trace_state=trace_state)
        return self.prediction_aggregator.finalize(
            pred_matrix=pred_matrix,
            substrate_list=substrate_list,
            request=request,
            trace_state=trace_state,
        )

    def predict_from_scoring_result(
        self,
        scoring_result: KinaseScoringResult,
        ensemble_size: int = 10,
        top: int = 50,
        score_threshold: float = 0.8,
        inclusion: int = 20,
        n_iterations: int = 5,
        random_state: int | None = None,
        allow_profile_only_fallback: bool = False,
        capture_debug_trace: bool = False,
        debug_kinases: list[str] | None = None,
        debug_top_n: int = 10,
        svm_mode: PredictionSvmMode | None = None,
        sampling_trace: PredictionSamplingTrace | str | Path | None = None,
        trace_level: PredictionTraceLevel | None = None,
        trace_sink: TraceSink | str | Path | None = None,
        trace_sink_format: PredictionTraceFormat = "csv",
    ) -> KinasePredictionResult:
        if scoring_result.combined_scores is not None:
            feature_mat = scoring_result.combined_scores
        elif allow_profile_only_fallback:
            feature_mat = scoring_result.profile_scores
        else:
            msg = (
                "scoring_result does not contain combined_scores; pass "
                "allow_profile_only_fallback=True to use profile_scores instead"
            )
            raise InputCompatibilityError(msg)

        return self.predict(
            combined_scores=feature_mat,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            capture_debug_trace=capture_debug_trace,
            debug_kinases=debug_kinases,
            debug_top_n=debug_top_n,
            svm_mode=svm_mode,
            sampling_trace=sampling_trace,
            trace_level=trace_level,
            trace_sink=trace_sink,
            trace_sink_format=trace_sink_format,
        )


def build_candidate_substrate_list(
    combined_scores: pd.DataFrame,
    top: int = 50,
    score_threshold: float = 0.8,
    inclusion: int = 20,
) -> dict[str, list[str]]:
    """Select candidate kinase substrates from the combined score matrix."""

    validate_positive_int(top, name="top")
    validate_probability_threshold(score_threshold, name="score_threshold")
    validate_positive_int(inclusion, name="inclusion")

    substrate_list: dict[str, list[str]] = {}
    for kinase in combined_scores.columns:
        selected = combined_scores.loc[:, kinase].nlargest(top)
        sites = selected.loc[selected > score_threshold].index.tolist()
        if len(sites) >= inclusion:
            substrate_list[kinase] = sites
    return substrate_list


__all__ = [
    "CandidateSelector",
    "EnsemblePredictor",
    "KinasePredictor",
    "NegativePoolSampler",
    "PredictionAggregator",
    "PredictionRequestFactory",
    "TraceRecorder",
    "build_candidate_substrate_list",
]
