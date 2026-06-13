"""Public workflow entrypoints."""

from phospy.workflows.differential.public import DifferentialAnalysisWorkflow
from phospy.workflows.enrichment.public import EnrichmentWorkflow
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.signalome.public import SignalomeWorkflow

__all__ = [
    "DifferentialAnalysisWorkflow",
    "EnrichmentWorkflow",
    "SignalomeWorkflow",
    "KinaseWorkflow",
]
