"""Prediction domain.

This package owns kinase prediction engines, sampling behaviour, scoring
components, prediction execution, and prediction result models.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

__all__ = [
    "EnsemblePredictorContract",
    "AdaptiveSamplingEnsembleTrace",
    "AdaptiveSamplingIterationTrace",
    "DEFAULT_KINASE_PROFILE_POLICY",
    "KinaseMotifScorer",
    "KinasePredictionDebugTrace",
    "KinasePredictionResult",
    "KinaseProfilePolicy",
    "KinaseProfileResult",
    "PredMatResult",
    "MotifScoringResult",
    "ValidatedMotifLibrary",
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
    "build_kinase_substrate_profiles",
    "build_validated_motif_library",
    "build_candidate_substrate_list",
    "combine_profile_and_motif_scores",
    "create_frequency_matrix",
    "frequency_scoring",
    "minmax_scale_columns",
    "prediction_debug_trace_tables",
    "resolve_prediction_sampling_policy",
    "score_phosphosite_motifs",
]

_EXPORTS: dict[str, tuple[str, str]] = {
    "EnsemblePredictorContract": (".contracts", "EnsemblePredictorContract"),
    "AdaptiveSamplingEnsembleTrace": (
        ".results",
        "AdaptiveSamplingEnsembleTrace",
    ),
    "AdaptiveSamplingIterationTrace": (
        ".results",
        "AdaptiveSamplingIterationTrace",
    ),
    "DEFAULT_KINASE_PROFILE_POLICY": (
        ".profiles",
        "DEFAULT_KINASE_PROFILE_POLICY",
    ),
    "KinaseMotifScorer": (".motif_scoring", "KinaseMotifScorer"),
    "KinasePredictionDebugTrace": (".results", "KinasePredictionDebugTrace"),
    "KinasePredictionResult": (".results", "KinasePredictionResult"),
    "KinaseProfilePolicy": (".profiles", "KinaseProfilePolicy"),
    "KinaseProfileResult": (".profiles", "KinaseProfileResult"),
    "PredMatResult": (".results", "PredMatResult"),
    "MotifScoringResult": (".motif_scoring", "MotifScoringResult"),
    "ValidatedMotifLibrary": (".motif_scoring", "ValidatedMotifLibrary"),
    "DirectoryTraceSink": (".traces", "DirectoryTraceSink"),
    "KinasePredictor": (".engines", "KinasePredictor"),
    "KinaseScorer": (".scoring", "KinaseScorer"),
    "KinaseScoringResult": (".scoring", "KinaseScoringResult"),
    "KinaseWorkflowExecutionResult": (".engines", "KinaseWorkflowExecutionResult"),
    "KinaseWorkflowExecutor": (".engines", "KinaseWorkflowExecutor"),
    "PredictionExecutionRunner": (".engines", "PredictionExecutionRunner"),
    "PredictionRequestFactory": (".engines", "PredictionRequestFactory"),
    "PredictionSamplingPolicy": (".policies", "PredictionSamplingPolicy"),
    "PredictionSamplingTrace": (".traces", "PredictionSamplingTrace"),
    "SamplingTraceOverrideEnsemble": (".results", "SamplingTraceOverrideEnsemble"),
    "TraceSink": (".traces", "TraceSink"),
    "build_kinase_substrate_profiles": (
        ".profiles",
        "build_kinase_substrate_profiles",
    ),
    "build_validated_motif_library": (
        ".motif_scoring",
        "build_validated_motif_library",
    ),
    "build_candidate_substrate_list": (".engines", "build_candidate_substrate_list"),
    "combine_profile_and_motif_scores": (
        ".scoring",
        "combine_profile_and_motif_scores",
    ),
    "create_frequency_matrix": (".motif_scoring", "create_frequency_matrix"),
    "frequency_scoring": (".motif_scoring", "frequency_scoring"),
    "minmax_scale_columns": (".motif_scoring", "minmax_scale_columns"),
    "prediction_debug_trace_tables": (".traces", "prediction_debug_trace_tables"),
    "resolve_prediction_sampling_policy": (
        ".policies",
        "resolve_prediction_sampling_policy",
    ),
    "score_phosphosite_motifs": (".motif_scoring", "score_phosphosite_motifs"),
}


def __getattr__(name: str) -> Any:
    if name not in _EXPORTS:
        msg = f"module {__name__!r} has no attribute {name!r}"
        raise AttributeError(msg)
    module_name, attribute_name = _EXPORTS[name]
    module = import_module(module_name, __name__)
    value = getattr(module, attribute_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
