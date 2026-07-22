"""Compatibility aggregate for scientific table schema imports."""

from phospy.frames.table_schema import TableSchema
from phospy.science.tables.activity import (
    ActivityCountSeries,
    ActivityMatrix,
    ActivityTargetTable,
)
from phospy.science.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.science.tables.differential import (
    filter_differential_results,
    rank_differential_results,
)
from phospy.science.tables.kinase import (
    KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS,
    KinasePredictionMatrix,
    KinaseScoreMatrix,
    KinaseSubstrateContributionTable,
)
from phospy.science.tables.references import (
    KinaseSubstrateReference,
    SiteSequenceReference,
)
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
