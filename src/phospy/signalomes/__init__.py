"""Signalome domain package."""

from phospy.signalomes.assignments import (
    build_module_assignments,
    select_kinase_substrates,
)
from phospy.signalomes.expanded import build_expanded_signalome_table
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeClusterCandidateScore,
    SignalomeModules,
    SignalomeModuleSelectionDiagnostics,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
)
from phospy.signalomes.modules import build_signalome_module_table
from phospy.signalomes.network import (
    build_kinase_network,
    build_kinase_network_with_diagnostics,
)

__all__ = [
    "KinaseNetwork",
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
