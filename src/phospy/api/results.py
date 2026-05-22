"""Public compatibility wrapper for internal contract ownership."""

from phospy.contracts import results as _result_contracts
from phospy.contracts.results import (
    DifferentialAnalysisResult,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
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
    "KinaseActivityResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowResult",
    "SignalomeWorkflowResult",
]
