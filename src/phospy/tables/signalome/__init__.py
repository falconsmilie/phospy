"""Signalome sidecar table schema wrappers."""

from phospy.tables.signalome.assignments import SignalomeAssignmentsTable
from phospy.tables.signalome.context import (
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
)
from phospy.tables.signalome.modules import SignalomeModulesTable
from phospy.tables.signalome.network import (
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
