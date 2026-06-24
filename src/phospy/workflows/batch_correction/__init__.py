"""Batch-correction workflow collaborators."""

from typing import TYPE_CHECKING

from phospy.workflows.batch_correction.contracts import (
    BatchCorrectionWorkflowRequest,
    BatchCorrectionWorkflowResult,
)
from phospy.workflows.batch_correction.interpreter import (
    BatchCorrectionDiagnosticRequirements,
    BatchCorrectionPlanInterpreter,
    EligibleControlSiteRow,
    ReplicateStructure,
    ResolvedBatchCorrectionPlan,
)
from phospy.workflows.batch_correction.provenance import (
    BatchCorrectionProvenanceRecorder,
)

if TYPE_CHECKING:
    from phospy.workflows.batch_correction.workflow import BatchCorrectionWorkflow


def __getattr__(name: str) -> object:
    if name == "BatchCorrectionWorkflow":
        from phospy.workflows.batch_correction.workflow import BatchCorrectionWorkflow

        return BatchCorrectionWorkflow
    raise AttributeError(name)


__all__ = [
    "BatchCorrectionDiagnosticRequirements",
    "BatchCorrectionPlanInterpreter",
    "BatchCorrectionProvenanceRecorder",
    "BatchCorrectionWorkflow",
    "BatchCorrectionWorkflowRequest",
    "BatchCorrectionWorkflowResult",
    "EligibleControlSiteRow",
    "ReplicateStructure",
    "ResolvedBatchCorrectionPlan",
]
