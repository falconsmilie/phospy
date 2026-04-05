from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..scoring import KinaseScoringResult
from ..types import PredictionSvmMode, PredictionTraceFormat, PredictionTraceLevel
from ..validation.errors import InputCompatibilityError
from ..validation.prediction import PredictionRequest
from .aggregation import PredictionAggregator
from .candidates import CandidateSelector, build_candidate_substrate_list
from .execution import EnsemblePredictor, NegativePoolSampler, TraceRecorder
from .models import KinasePredictionResult
from .traces import (
    PredictionSamplingTrace,
    TraceSink,
    build_prediction_runtime_session,
)
from .validation import validate_svm_mode


class PredictionRequestFactory:
    def __init__(self, *, default_svm_mode: PredictionSvmMode) -> None:
        self.default_svm_mode = validate_svm_mode(default_svm_mode)

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
        return PredictionRequest.validate_request(
            combined_scores=combined_scores,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            debug_kinases=debug_kinases,
            debug_top_n=debug_top_n,
            svm_mode=svm_mode,
            sampling_trace=sampling_trace,
            trace_level=trace_level,
            trace_sink=trace_sink,
            default_svm_mode=self.default_svm_mode,
            capture_debug_trace=capture_debug_trace,
            trace_sink_format=trace_sink_format,
        )


class PredictionExecutionRunner:
    """Execute a validated prediction request against resolved runtime resources."""

    def __init__(
        self,
        *,
        candidate_selector: CandidateSelector,
        prediction_aggregator: PredictionAggregator,
        trace_recorder: TraceRecorder,
        ensemble_predictor: EnsemblePredictor,
    ) -> None:
        self.candidate_selector = candidate_selector
        self.prediction_aggregator = prediction_aggregator
        self.trace_recorder = trace_recorder
        self.ensemble_predictor = ensemble_predictor

    def run(self, request: PredictionRequest) -> KinasePredictionResult:
        substrate_list = self.candidate_selector.select(
            request.combined_scores,
            top=request.top,
            score_threshold=request.score_threshold,
            inclusion=request.inclusion,
        )
        if not substrate_list:
            return self.prediction_aggregator.empty_result(request=request)

        feature_mat = request.combined_scores
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
        self.request_factory = request_factory or PredictionRequestFactory(
            default_svm_mode=svm_mode
        )
        self.svm_mode = self.request_factory.default_svm_mode
        self.candidate_selector = candidate_selector or CandidateSelector()
        self.negative_pool_sampler = negative_pool_sampler or NegativePoolSampler()
        self.trace_recorder = trace_recorder or TraceRecorder()
        self.prediction_aggregator = prediction_aggregator or PredictionAggregator()
        self.ensemble_predictor = ensemble_predictor or EnsemblePredictor(
            kernel=self.kernel,
            negative_pool_sampler=self.negative_pool_sampler,
            trace_recorder=self.trace_recorder,
        )

    def _build_execution_runner(self) -> PredictionExecutionRunner:
        return PredictionExecutionRunner(
            candidate_selector=self.candidate_selector,
            prediction_aggregator=self.prediction_aggregator,
            trace_recorder=self.trace_recorder,
            ensemble_predictor=self.ensemble_predictor,
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
        return self.predict_request(request)

    def predict_request(self, request: PredictionRequest) -> KinasePredictionResult:
        with build_prediction_runtime_session(request) as runtime_session:
            result = self._build_execution_runner().run(runtime_session.runtime_request)
            return runtime_session.transfer_trace_sink_ownership(result)

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

        request = self.request_factory.create(
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
        return self.predict_request(request)


__all__ = [
    "CandidateSelector",
    "EnsemblePredictor",
    "KinasePredictor",
    "NegativePoolSampler",
    "PredictionAggregator",
    "PredictionExecutionRunner",
    "PredictionRequestFactory",
    "TraceRecorder",
    "build_candidate_substrate_list",
]
