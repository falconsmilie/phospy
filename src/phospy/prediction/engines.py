from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    format_no_candidate_kinases_message,
)
from ..internal.defaults import (
    DEFAULT_MOTIF_FLANK_SIZE,
    DEFAULT_PREDICTION_DEBUG_TOP_N,
    DEFAULT_PREDICTION_ENSEMBLE_SIZE,
    DEFAULT_PREDICTION_INCLUSION,
    DEFAULT_PREDICTION_N_ITERATIONS,
    DEFAULT_PREDICTION_SCORE_THRESHOLD,
    DEFAULT_PREDICTION_TOP,
)
from ..internal.types import (
    PREDICTION_SVM_MODE_DEFAULT,
    PREDICTION_TRACE_FORMAT_CSV,
    PredictionSvmMode,
    PredictionTraceFormat,
    PredictionTraceLevel,
)
from ..references import ReferenceBundle
from ..validation.domain.prediction import (
    validate_ensemble_predictor,
    validate_kinase_prediction_batch,
)
from ..validation.requests.prediction import PredictionRequest
from ..validation.requests.workflow import (
    WorkflowInputs,
    validate_workflow_request,
)
from .aggregation import PredictionAggregator
from .candidates import CandidateSelector, build_candidate_substrate_list
from .contracts import EnsemblePredictorContract
from .execution import (
    EnsemblePredictor,
    NegativePoolSampler,
    PredictionSamplingSession,
    TraceRecorder,
)
from .motif_scoring import MotifScoringResult
from .profiles import (
    KinaseProfilePolicy,
    KinaseProfileResult,
    build_kinase_substrate_profiles,
)
from .results import KinasePredictionResult
from .scoring import KinaseScorer, KinaseScoringResult
from .traces import (
    PredictionSamplingTrace,
    TraceSink,
    build_prediction_runtime_session,
)
from .validation import validate_svm_mode


class PredictionRequestFactory:
    """Construct validated prediction requests from public predictor inputs."""

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
        ensemble_predictor: EnsemblePredictorContract | None,
    ) -> None:
        self.candidate_selector = candidate_selector
        self.prediction_aggregator = prediction_aggregator
        self.trace_recorder = trace_recorder
        self.ensemble_predictor = validate_ensemble_predictor(ensemble_predictor)

    @staticmethod
    def _merge_cleanup_errors(
        *,
        clear_cache_error: BaseException | None,
        flush_final_error: BaseException | None,
    ) -> BaseException | None:
        if clear_cache_error is None:
            return flush_final_error
        if flush_final_error is None:
            return clear_cache_error
        return RuntimeError(
            "Prediction cleanup failed: clear_cache raised "
            f"{clear_cache_error!r}; flush_final raised {flush_final_error!r}"
        )

    def _cleanup_prediction_resources(
        self,
        *,
        trace_state,
    ) -> BaseException | None:
        clear_cache_error: BaseException | None = None
        flush_final_error: BaseException | None = None

        try:
            self.ensemble_predictor.clear_cache()
        except BaseException as error:
            clear_cache_error = error

        try:
            self.trace_recorder.flush_final(trace_state=trace_state)
        except BaseException as error:
            flush_final_error = error

        return self._merge_cleanup_errors(
            clear_cache_error=clear_cache_error,
            flush_final_error=flush_final_error,
        )

    def run(self, request: PredictionRequest) -> KinasePredictionResult:
        substrate_list = self.candidate_selector.select(
            request.combined_scores,
            top=request.top,
            score_threshold=request.score_threshold,
            inclusion=request.inclusion,
        )
        if not substrate_list:
            self.prediction_aggregator.raise_no_candidate_kinases(request=request)

        if self.ensemble_predictor is None:
            msg = "ensemble_predictor is required when candidate kinases are present"
            raise InputCompatibilityError(msg)

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
        sampling_session = PredictionSamplingSession.from_request(request)
        execution_error: BaseException | None = None
        try:
            for kinase, substrates in substrate_list.items():
                batch = self.ensemble_predictor.predict_kinase(
                    kinase=kinase,
                    substrates=substrates,
                    feature_mat=feature_mat,
                    request=request,
                    trace_state=trace_state,
                    sampling_session=sampling_session,
                )
                batch = validate_kinase_prediction_batch(
                    batch=batch,
                    requested_kinase=kinase,
                    feature_index=feature_mat.index,
                )
                self.prediction_aggregator.add_kinase_scores(
                    pred_matrix=pred_matrix,
                    batch=batch,
                )
        except BaseException as error:
            execution_error = error

        cleanup_error = self._cleanup_prediction_resources(trace_state=trace_state)
        if execution_error is not None:
            if cleanup_error is not None:
                raise execution_error from cleanup_error
            raise execution_error
        if cleanup_error is not None:
            raise cleanup_error

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

    Use ``svm_mode='r_parity'`` when you want the closest supported
    parity-oriented preset, including the R-like learner, sampling, and
    final-scoring seams. The default mode preserves the recommended stable,
    column-order-invariant behaviour.
    """

    def __init__(
        self,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
        *,
        request_factory: PredictionRequestFactory | None = None,
        candidate_selector: CandidateSelector | None = None,
        negative_pool_sampler: NegativePoolSampler | None = None,
        trace_recorder: TraceRecorder | None = None,
        prediction_aggregator: PredictionAggregator | None = None,
        ensemble_predictor: EnsemblePredictorContract | None = None,
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
        ensemble_size: int = DEFAULT_PREDICTION_ENSEMBLE_SIZE,
        top: int = DEFAULT_PREDICTION_TOP,
        score_threshold: float = DEFAULT_PREDICTION_SCORE_THRESHOLD,
        inclusion: int = DEFAULT_PREDICTION_INCLUSION,
        n_iterations: int = DEFAULT_PREDICTION_N_ITERATIONS,
        random_state: int | None = None,
        capture_debug_trace: bool = False,
        debug_kinases: list[str] | None = None,
        debug_top_n: int = DEFAULT_PREDICTION_DEBUG_TOP_N,
        svm_mode: PredictionSvmMode | None = None,
        sampling_trace: PredictionSamplingTrace | str | Path | None = None,
        trace_level: PredictionTraceLevel | None = None,
        trace_sink: TraceSink | str | Path | None = None,
        trace_sink_format: PredictionTraceFormat = PREDICTION_TRACE_FORMAT_CSV,
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
        ensemble_size: int = DEFAULT_PREDICTION_ENSEMBLE_SIZE,
        top: int = DEFAULT_PREDICTION_TOP,
        score_threshold: float = DEFAULT_PREDICTION_SCORE_THRESHOLD,
        inclusion: int = DEFAULT_PREDICTION_INCLUSION,
        n_iterations: int = DEFAULT_PREDICTION_N_ITERATIONS,
        random_state: int | None = None,
        allow_profile_only_fallback: bool = False,
        capture_debug_trace: bool = False,
        debug_kinases: list[str] | None = None,
        debug_top_n: int = DEFAULT_PREDICTION_DEBUG_TOP_N,
        svm_mode: PredictionSvmMode | None = None,
        sampling_trace: PredictionSamplingTrace | str | Path | None = None,
        trace_level: PredictionTraceLevel | None = None,
        trace_sink: TraceSink | str | Path | None = None,
        trace_sink_format: PredictionTraceFormat = PREDICTION_TRACE_FORMAT_CSV,
    ) -> KinasePredictionResult:
        if scoring_result.combined_scores is not None:
            feature_mat = scoring_result.combined_scores
            feature_name = "combined_scores"
        elif allow_profile_only_fallback:
            feature_mat = scoring_result.profile_scores
            feature_name = "profile_scores"
        else:
            msg = (
                "scoring_result does not contain combined_scores; pass "
                "allow_profile_only_fallback=True to use profile_scores instead"
            )
            raise InputCompatibilityError(msg)

        feature_values = feature_mat.to_numpy(dtype=float)
        if not np.isfinite(feature_values).all():
            if not np.isfinite(feature_values).any():
                msg = format_no_candidate_kinases_message(
                    source_name=feature_name,
                    top=top,
                    score_threshold=score_threshold,
                    inclusion=inclusion,
                    kinase_count=int(feature_mat.shape[1]),
                    site_count=int(feature_mat.shape[0]),
                    effective_top=min(int(top), int(feature_mat.shape[0])),
                    qualifying_kinases=0,
                    max_qualifying_sites=0,
                )
                msg = (
                    f"{msg} Candidate scoring could not proceed because all values in "
                    f"{feature_name} are non-finite."
                )
                raise NoCandidateKinasesError(msg)
            msg = (
                "scoring_result contains non-finite values in "
                f"{feature_name}; regenerate scores with finite inputs before "
                "prediction"
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
        try:
            return self.predict_request(request)
        except NoCandidateKinasesError as error:
            if feature_name == "profile_scores":
                msg = (
                    f"{error} Fallback path considered: profile_scores was used "
                    "because combined_scores was unavailable."
                )
                raise NoCandidateKinasesError(msg) from error
            raise


@dataclass(frozen=True, slots=True)
class KinaseWorkflowExecutionResult:
    """Trusted execution outputs for one validated kinase workflow request."""

    profile_result: KinaseProfileResult
    motif_result: MotifScoringResult | None
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult


class KinaseWorkflowExecutor:
    """Prediction-domain executor for the validated kinase workflow path."""

    def __init__(
        self,
        *,
        flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
    ) -> None:
        self.flank_size = flank_size
        self.kernel = kernel
        self.svm_mode = svm_mode

    def validate_request(
        self,
        *,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
        min_substrates: int = 1,
        min_motif_size: int = 1,
        allow_profile_only_fallback: bool = False,
        ensemble_size: int = DEFAULT_PREDICTION_ENSEMBLE_SIZE,
        top: int = DEFAULT_PREDICTION_TOP,
        score_threshold: float = DEFAULT_PREDICTION_SCORE_THRESHOLD,
        inclusion: int = DEFAULT_PREDICTION_INCLUSION,
        n_iterations: int = DEFAULT_PREDICTION_N_ITERATIONS,
        random_state: int | None = None,
        svm_mode: PredictionSvmMode | None = None,
        profile_policy: KinaseProfilePolicy | None = None,
    ) -> WorkflowInputs:
        return validate_workflow_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            reference_bundle=reference_bundle,
            min_substrates=min_substrates,
            min_motif_size=min_motif_size,
            allow_profile_only_fallback=allow_profile_only_fallback,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            svm_mode=svm_mode,
            profile_policy=profile_policy,
            flank_size=self.flank_size,
            default_svm_mode=self.svm_mode,
        )

    def execute_validated_request(
        self,
        request: WorkflowInputs,
    ) -> KinaseWorkflowExecutionResult:
        phospho_matrix = request.phospho_matrix

        profile_result = build_kinase_substrate_profiles(
            substrate_map=request.substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=request.min_substrates,
            policy=request.profile_policy,
        )

        scorer = KinaseScorer(profile_result.profile_matrix)
        motif_result: MotifScoringResult | None = None
        scoring_matrix = phospho_matrix

        if request.motif_scorer is not None:
            scoring_matrix = phospho_matrix.loc[list(request.scoring_site_index)]
            motif_result = request.motif_scorer.score_sequences(
                seqs=request.site_sequences,
                site_index=request.scoring_site_index,
                min_motif_size=request.min_motif_size,
            )
            scoring_result = scorer.score(
                phospho_matrix=scoring_matrix,
                motif_scores=motif_result.motif_scores,
                motif_sizes=motif_result.motif_sizes,
                profile_sizes=profile_result.substrate_counts.astype(float),
                allow_profile_only_fallback=request.allow_profile_only_fallback,
            )
        else:
            scoring_result = scorer.score(phospho_matrix=scoring_matrix)

        predictor = KinasePredictor(
            kernel=self.kernel,
            svm_mode=request.predictor_svm_mode,
        )
        prediction_result = predictor.predict_from_scoring_result(
            scoring_result=scoring_result,
            ensemble_size=request.ensemble_size,
            top=request.top,
            score_threshold=request.score_threshold,
            inclusion=request.inclusion,
            n_iterations=request.n_iterations,
            random_state=request.random_state,
            allow_profile_only_fallback=request.allow_profile_only_fallback,
            svm_mode=request.svm_mode,
        )

        return KinaseWorkflowExecutionResult(
            profile_result=profile_result,
            motif_result=motif_result,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
        )


__all__ = [
    "CandidateSelector",
    "EnsemblePredictor",
    "KinasePredictor",
    "KinaseWorkflowExecutionResult",
    "KinaseWorkflowExecutor",
    "NegativePoolSampler",
    "PredictionAggregator",
    "PredictionExecutionRunner",
    "PredictionRequestFactory",
    "TraceRecorder",
    "build_candidate_substrate_list",
]
