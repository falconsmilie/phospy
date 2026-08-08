"""Signalome sidecar table schema wrappers."""

__phospy_contracts_facade_role__ = "science_owned_public_table_contract"

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
