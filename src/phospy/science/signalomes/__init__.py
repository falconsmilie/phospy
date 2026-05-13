"""Signalome domain package."""

from phospy.science.signalomes.assignments import (
    build_module_assignments,
    select_kinase_substrates,
)
from phospy.science.signalomes.expanded import build_expanded_signalome_table
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAlignmentDiagnostics,
    SignalomeAlignmentInputDiagnostics,
    SignalomeAssignments,
    SignalomeClusterCandidateScore,
    SignalomeModules,
    SignalomeModuleSelectionDiagnostics,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
)
from phospy.science.signalomes.modules import build_signalome_module_table
from phospy.science.signalomes.network import (
    build_kinase_network,
    build_kinase_network_with_diagnostics,
)

__all__ = [
    "KinaseNetwork",
    "SignalomeAlignmentDiagnostics",
    "SignalomeAlignmentInputDiagnostics",
    "SignalomeAssignments",
    "SignalomeClusterCandidateScore",
    "SignalomeModuleSelectionDiagnostics",
    "SignalomeNetworkCorrelationDiagnostics",
    "SignalomeModules",
    "build_expanded_signalome_table",
    "build_kinase_network",
    "build_kinase_network_with_diagnostics",
    "build_module_assignments",
    "build_signalome_module_table",
    "select_kinase_substrates",
    "SignalomeScorePreconditioningDiagnostics",
]
