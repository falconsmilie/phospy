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
    """Advanced collaborator graph for `SimpleKinaseWorkflow`.

    This is an extension seam for advanced users who need to replace default
    workflow collaborators as a unit. The graph is treated as a stable public
    composition contract, while concrete collaborator implementations may evolve
    independently.
    """

    analysis_ready_builder: AnalysisReadyDatasetBuilder
    reference_provider: ReferenceProvider
    activity_analyzer: KinaseActivityAnalyzer
    workflow_executor: KinaseWorkflowExecutor

    def __post_init__(self) -> None:
        required_collaborators = {
            "analysis_ready_builder": (self.analysis_ready_builder, ("build",)),
            "reference_provider": (self.reference_provider, ("resolve",)),
            "activity_analyzer": (self.activity_analyzer, ("run",)),
            "workflow_executor": (
                self.workflow_executor,
                ("validate_request", "execute_validated_request"),
            ),
        }
        for collaborator_name, (
            collaborator,
            methods,
        ) in required_collaborators.items():
            if collaborator is None:
                msg = (
                    f"SimpleKinaseExecutionGraph collaborator "
                    f"'{collaborator_name}' cannot be None."
                )
                raise ValueError(msg)
            for method_name in methods:
                if not hasattr(collaborator, method_name):
                    msg = (
                        f"SimpleKinaseExecutionGraph collaborator "
                        f"'{collaborator_name}' must define '{method_name}()'."
                    )
                    raise TypeError(msg)


def create_default_simple_kinase_execution_graph(
    *,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
    kernel: str = "rbf",
    svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
) -> SimpleKinaseExecutionGraph:
    """Build the default `SimpleKinaseExecutionGraph` used by public workflows."""

    return SimpleKinaseExecutionGraph(
        analysis_ready_builder=AnalysisReadyDatasetBuilder(),
        reference_provider=BundledReferenceProvider(),
        activity_analyzer=KinaseActivityAnalyzer(),
        workflow_executor=KinaseWorkflowExecutor(
            flank_size=flank_size,
            kernel=kernel,
            svm_mode=svm_mode,
        ),
    )
