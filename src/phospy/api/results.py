"""Public compatibility wrapper for internal contract ownership."""

from phospy.contracts import results as _result_contracts
from phospy.contracts.results import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_NOT_REPORTED,
    IMPORTER_QUALITY_STATUS_REPORTED,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
    DifferentialAnalysisResult,
    DifferentialContrastDefinition,
    DifferentialDesignMatrixSummary,
    DifferentialEmpiricalBayesProvenance,
    DifferentialFixedEffectCovariateProvenance,
    DifferentialMissingValuePolicyProvenance,
    DifferentialPolicyProvenance,
    DifferentialReplicatePolicyProvenance,
    DifferentialStatisticalTestingProvenance,
    DifferentialTechnicalReplicateGroup,
    DifferentialUnsupportedDesignPolicyProvenance,
    EnrichmentResultRecord,
    EnrichmentWorkflowResult,
    ImporterDetectedIntensityColumn,
    ImporterDuplicateKeySummary,
    ImporterFlaggedRowSummary,
    ImporterLocalisationConfidenceSummary,
    ImporterMissingIntensitySummary,
    ImporterQualityCount,
    ImporterQualityReport,
    ImporterQualityStatus,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    PhosphositeImportResult,
    ProteinAwareMappingDiagnostics,
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
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
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "DifferentialAnalysisResult",
    "DifferentialContrastDefinition",
    "DifferentialDesignMatrixSummary",
    "DifferentialEmpiricalBayesProvenance",
    "DifferentialFixedEffectCovariateProvenance",
    "DifferentialMissingValuePolicyProvenance",
    "DifferentialPolicyProvenance",
    "DifferentialReplicatePolicyProvenance",
    "DifferentialStatisticalTestingProvenance",
    "DifferentialTechnicalReplicateGroup",
    "DifferentialUnsupportedDesignPolicyProvenance",
    "EnrichmentResultRecord",
    "EnrichmentWorkflowResult",
    "IMPORTER_QUALITY_STATUS_NOT_APPLICABLE",
    "IMPORTER_QUALITY_STATUS_NOT_REPORTED",
    "IMPORTER_QUALITY_STATUS_REPORTED",
    "ImporterDetectedIntensityColumn",
    "ImporterDuplicateKeySummary",
    "ImporterFlaggedRowSummary",
    "ImporterLocalisationConfidenceSummary",
    "ImporterMissingIntensitySummary",
    "ImporterQualityCount",
    "ImporterQualityReport",
    "ImporterQualityStatus",
    "KinaseActivityResult",
    "KinaseEligibilityReport",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowPreprocessingAttritionSummary",
    "KinaseWorkflowResult",
    "KinaseWorkflowScoringAttritionSummary",
    "KinaseWorkflowSiteAttritionSummary",
    "PhosphositeImportResult",
    "ProteinAwareMappingDiagnostics",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
    "ProteinAwareSiteEligibility",
    "SignalomeWorkflowResult",
]
