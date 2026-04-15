from __future__ import annotations

from dataclasses import dataclass

from ..activities.analysis import KinaseActivityAnalyzer
from ..internal.defaults import DEFAULT_MOTIF_FLANK_SIZE
from ..internal.types import PREDICTION_SVM_MODE_DEFAULT, PredictionSvmMode
from ..prediction.engines import KinaseWorkflowExecutor
from ..preprocessing.modes import AnalysisReadyDatasetBuilder
from ..references import BundledReferenceProvider, ReferenceProvider

__all__ = [
    "SimpleKinaseExecutionGraph",
    "create_default_simple_kinase_execution_graph",
]


@dataclass(frozen=True, slots=True)
class SimpleKinaseExecutionGraph:
    analysis_ready_builder: AnalysisReadyDatasetBuilder
    reference_provider: ReferenceProvider
    activity_analyzer: KinaseActivityAnalyzer
    workflow_executor: KinaseWorkflowExecutor


def create_default_simple_kinase_execution_graph(
    *,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
    analysis_ready_builder: AnalysisReadyDatasetBuilder | None = None,
    reference_provider: ReferenceProvider | None = None,
    activity_analyzer: KinaseActivityAnalyzer | None = None,
    workflow_executor: KinaseWorkflowExecutor | None = None,
) -> SimpleKinaseExecutionGraph:
    return SimpleKinaseExecutionGraph(
        analysis_ready_builder=(
            AnalysisReadyDatasetBuilder()
            if analysis_ready_builder is None
            else analysis_ready_builder
        ),
        reference_provider=(
            BundledReferenceProvider()
            if reference_provider is None
            else reference_provider
        ),
        activity_analyzer=(
            KinaseActivityAnalyzer() if activity_analyzer is None else activity_analyzer
        ),
        workflow_executor=(
            KinaseWorkflowExecutor(
                flank_size=flank_size,
                kernel=kernel,
                svm_mode=svm_mode,
            )
            if workflow_executor is None
            else workflow_executor
        ),
    )
