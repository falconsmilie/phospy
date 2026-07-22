"""Compatibility re-exports for signalome table schemas."""

from phospy.science.tables.signalome import (
    KinaseNetworkCandidateCorrelationsTable,
    KinaseNetworkEdgesTable,
    KinaseNetworkNodesTable,
    SignalomeAssignmentsTable,
    SignalomeModulesTable,
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
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
