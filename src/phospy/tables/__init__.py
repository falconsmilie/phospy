"""Internal scientific table schema wrappers."""

from phospy.tables.activity import (
    ActivityCountSeries,
    ActivityMatrix,
    ActivityTargetTable,
)
from phospy.tables.base import TableSchema
from phospy.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.tables.differential import (
    filter_differential_results,
    rank_differential_results,
)
from phospy.tables.kinase import (
    KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS,
    KinasePredictionMatrix,
    KinaseScoreMatrix,
    KinaseSubstrateContributionTable,
)
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference
from phospy.tables.signalome import (
    KinaseNetworkCandidateCorrelationsTable,
    KinaseNetworkEdgesTable,
    KinaseNetworkNodesTable,
    SignalomeAssignmentsTable,
    SignalomeModulesTable,
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
)

__all__ = [
    "ActivityCountSeries",
    "ActivityMatrix",
    "ActivityTargetTable",
    "KinasePredictionMatrix",
    "KinaseNetworkCandidateCorrelationsTable",
    "KinaseNetworkEdgesTable",
    "KinaseNetworkNodesTable",
    "KinaseScoreMatrix",
    "KinaseSubstrateContributionTable",
    "KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS",
    "KinaseSubstrateReference",
    "PhosphoIntensityMatrix",
    "filter_differential_results",
    "rank_differential_results",
    "SampleMetadataTable",
    "SignalomeAssignmentsTable",
    "SignalomeModulesTable",
    "SignalomeProteinSiteContext",
    "SignalomeSiteContext",
    "SiteMetadataTable",
    "SiteSequenceReference",
    "TableSchema",
    "TotalProteinMatrix",
]
