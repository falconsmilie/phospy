from __future__ import annotations

import inspect

from phospy.api.simple_workflow_composition import (
    SimpleKinaseExecutionGraph,
    create_default_simple_kinase_execution_graph,
)
from phospy.prediction.contracts import EnsemblePredictorContract
from phospy.prediction.engines import PredictionRequestFactory
from phospy.prediction.execution import (
    EnsemblePredictor,
    KinasePredictionBatch,
    NegativePoolSampler,
    PredictionSamplingSession,
    PredictionTraceState,
    TraceRecorder,
)
from phospy.prediction.results import (
    AdaptiveSamplingEnsembleTrace,
    AdaptiveSamplingIterationTrace,
    KinasePredictionDebugTrace,
    SamplingTraceOverrideEnsemble,
)
from phospy.prediction.trace_replay import PredictionSamplingTrace
from phospy.prediction.trace_runtime import (
    DirectoryTraceSink,
    PredictionExecutionContext,
    PredictionRuntimeSession,
    TraceSink,
    build_prediction_runtime_session,
    create_trace_sink,
)


def test_advanced_prediction_and_workflow_seams_have_minimal_docstrings() -> None:
    documented_objects = (
        EnsemblePredictorContract,
        PredictionRequestFactory,
        KinasePredictionBatch,
        PredictionTraceState,
        PredictionSamplingSession,
        NegativePoolSampler,
        TraceRecorder,
        EnsemblePredictor,
        AdaptiveSamplingIterationTrace,
        AdaptiveSamplingEnsembleTrace,
        KinasePredictionDebugTrace,
        SamplingTraceOverrideEnsemble,
        PredictionSamplingTrace,
        TraceSink,
        DirectoryTraceSink,
        PredictionExecutionContext,
        PredictionRuntimeSession,
        create_trace_sink,
        build_prediction_runtime_session,
        SimpleKinaseExecutionGraph,
        create_default_simple_kinase_execution_graph,
    )

    undocumented = [
        object_.__qualname__
        for object_ in documented_objects
        if not inspect.getdoc(object_)
    ]
    assert undocumented == []
