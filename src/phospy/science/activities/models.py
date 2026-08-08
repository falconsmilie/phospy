"""Compatibility import route for activity-domain result models.

Owned implementations are split by responsibility under this package. This
module preserves the historical ``phospy.science.activities.models`` import path
by re-exporting the same class and constant objects.
"""

__phospy_contracts_facade_role__ = "science_owned_public_model"

from phospy.science.activities.diagnostics import (
    ActivityMethodDiagnostics,
    KseaZScoreActivityDiagnostics,
    SsgseaSubstrateEnrichmentActivityDiagnostics,
    WeightedSubstrateActivityDiagnostics,
)
from phospy.science.activities.inputs import KinaseActivityInputs, PredMatOverlapSummary
from phospy.science.activities.method_models import (
    KSEA_ZSCORE_ACTIVITY_METHOD,
    SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD,
    SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD,
    ActivityMethodMetadata,
    ActivityMethodSummary,
)
from phospy.science.activities.results import KinaseActivityResult

__all__ = [
    "ActivityMethodDiagnostics",
    "ActivityMethodMetadata",
    "ActivityMethodSummary",
    "KSEA_ZSCORE_ACTIVITY_METHOD",
    "KinaseActivityInputs",
    "KinaseActivityResult",
    "KseaZScoreActivityDiagnostics",
    "PredMatOverlapSummary",
    "SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD",
    "SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD",
    "SsgseaSubstrateEnrichmentActivityDiagnostics",
    "WeightedSubstrateActivityDiagnostics",
]
