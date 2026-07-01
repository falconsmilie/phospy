"""Public compatibility wrapper for internal contract ownership."""

from phospy.contracts import results as _result_contracts
from phospy.contracts.results import (
    DifferentialAnalysisResult,
    DifferentialModelDiagnostics,
    EnrichmentResultRecord,
    EnrichmentWorkflowResult,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    PhosphositeImportResult,
    SignalomeWorkflowResult,
)

# Compatibility aliases intentionally re-exported at module scope.
KinaseEligibilityReport = _result_contracts.KinaseEligibilityReport
KinaseWorkflowPreprocessingAttritionSummary = (
    _result_contracts.KinaseWorkflowPreprocessingAttritionSummary
)
KinaseWorkflowScoringAttritionSummary = (
    _result_contracts.KinaseWorkflowScoringAttritionSummary
)
KinaseWorkflowSiteAttritionSummary = (
    _result_contracts.KinaseWorkflowSiteAttritionSummary
)

__all__ = [
    "DifferentialAnalysisResult",
    "DifferentialModelDiagnostics",
    "EnrichmentResultRecord",
    "EnrichmentWorkflowResult",
    "KinaseActivityResult",
    "KinaseEligibilityReport",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowPreprocessingAttritionSummary",
    "KinaseWorkflowResult",
    "KinaseWorkflowScoringAttritionSummary",
    "KinaseWorkflowSiteAttritionSummary",
    "PhosphositeImportResult",
    "SignalomeWorkflowResult",
]
