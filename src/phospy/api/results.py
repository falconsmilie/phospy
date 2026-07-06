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
    ResultCaveat,
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
KinaseWorkflowAttritionProvenance = _result_contracts.KinaseWorkflowAttritionProvenance
KinaseWorkflowCaveat = _result_contracts.KinaseWorkflowCaveat

__all__ = [
    "DifferentialAnalysisResult",
    "DifferentialModelDiagnostics",
    "EnrichmentResultRecord",
    "EnrichmentWorkflowResult",
    "KinaseActivityResult",
    "KinaseEligibilityReport",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowAttritionProvenance",
    "KinaseWorkflowCaveat",
    "KinaseWorkflowPreprocessingAttritionSummary",
    "KinaseWorkflowResult",
    "KinaseWorkflowScoringAttritionSummary",
    "KinaseWorkflowSiteAttritionSummary",
    "PhosphositeImportResult",
    "ResultCaveat",
    "SignalomeWorkflowResult",
]
