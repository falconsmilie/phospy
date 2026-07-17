"""Signalome sidecar table schema wrappers."""

from phospy.science.tables.signalome.assignments import SignalomeAssignmentsTable
from phospy.science.tables.signalome.context import (
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
)
from phospy.science.tables.signalome.modules import SignalomeModulesTable
from phospy.science.tables.signalome.network import (
    KinaseNetworkCandidateCorrelationsTable,
    KinaseNetworkEdgesTable,
    KinaseNetworkNodesTable,
)

__all__ = [
    "KinaseNetworkCandidateCorrelationsTable",
    "KinaseNetworkEdgesTable",
    "KinaseNetworkNodesTable",
    "SignalomeAssignmentsTable",
    "SignalomeModulesTable",
    "SignalomeProteinSiteContext",
    "SignalomeSiteContext",
]
