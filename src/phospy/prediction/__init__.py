"""Prediction domain.

This package owns kinase prediction engines, sampling behaviour, scoring
components, prediction execution, and prediction result models."""

from __future__ import annotations

from .contracts import EnsemblePredictorContract
from .engines import (
    KinasePredictor,
    KinaseWorkflowExecutionResult,
    KinaseWorkflowExecutor,
    PredictionExecutionRunner,
    PredictionRequestFactory,
    build_candidate_substrate_list,
)
from .policies import PredictionSamplingPolicy, resolve_prediction_sampling_policy
from .results import (
    AdaptiveSamplingEnsembleTrace,
    AdaptiveSamplingIterationTrace,
    KinasePredictionDebugTrace,
    KinasePredictionResult,
    PredMatResult,
    SamplingTraceOverrideEnsemble,
)
from .scoring import (
    KinaseScorer,
    KinaseScoringResult,
    combine_profile_and_motif_scores,
)
from .traces import (
    DirectoryTraceSink,
    PredictionSamplingTrace,
    TraceSink,
    prediction_debug_trace_tables,
)

__all__ = [
    "EnsemblePredictorContract",
    "AdaptiveSamplingEnsembleTrace",
    "AdaptiveSamplingIterationTrace",
    "KinasePredictionDebugTrace",
    "KinasePredictionResult",
    "PredMatResult",
    "DirectoryTraceSink",
    "KinasePredictor",
    "KinaseScorer",
    "KinaseScoringResult",
    "KinaseWorkflowExecutionResult",
    "KinaseWorkflowExecutor",
    "PredictionExecutionRunner",
    "PredictionRequestFactory",
    "PredictionSamplingPolicy",
    "PredictionSamplingTrace",
    "SamplingTraceOverrideEnsemble",
    "TraceSink",
    "build_candidate_substrate_list",
    "combine_profile_and_motif_scores",
    "prediction_debug_trace_tables",
    "resolve_prediction_sampling_policy",
]
