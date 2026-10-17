"""Signalome sidecar table schema wrappers."""

from phospy.signalome.tables.assignments import SignalomeAssignmentsTable
from phospy.signalome.tables.context import (
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
)
from phospy.signalome.tables.modules import SignalomeModulesTable
from phospy.signalome.tables.network import (
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
