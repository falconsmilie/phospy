"""Public result models."""

from __future__ import annotations

from phospy.contracts.results.base import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_NOT_REPORTED,
    IMPORTER_QUALITY_STATUS_REPORTED,
    ImporterDetectedIntensityColumn,
    ImporterDuplicateKeySummary,
    ImporterFlaggedRowSummary,
    ImporterLocalisationConfidenceSummary,
    ImporterMissingIntensitySummary,
    ImporterQualityCount,
    ImporterQualityReport,
    ImporterQualityStatus,
    PhosphositeImportResult,
)
from phospy.contracts.results.differential import (
    DifferentialAnalysisResult,
    DifferentialContrastDefinition,
    DifferentialDesignMatrixSummary,
    DifferentialEmpiricalBayesProvenance,
    DifferentialFixedEffectCovariateProvenance,
    DifferentialMissingValuePolicyProvenance,
    DifferentialModelDiagnostics,
    DifferentialPolicyProvenance,
    DifferentialReplicatePolicyProvenance,
    DifferentialStatisticalTestingProvenance,
    DifferentialTechnicalReplicateGroup,
    DifferentialUnsupportedDesignPolicyProvenance,
)
from phospy.contracts.results.enrichment import (
    EnrichmentResultRecord,
    EnrichmentWorkflowResult,
)
from phospy.contracts.results.kinase import (
    ActivityMethodDiagnostics,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    KseaZScoreActivityDiagnostics,
    SsgseaSubstrateEnrichmentActivityDiagnostics,
    WeightedSubstrateActivityDiagnostics,
)
from phospy.contracts.results.kinase import (
    KinaseEligibilityReport as KinaseEligibilityReport,
)
from phospy.contracts.results.kinase import (
    KinaseWorkflowAttritionProvenance as KinaseWorkflowAttritionProvenance,
)
from phospy.contracts.results.kinase import (
    KinaseWorkflowCaveat as KinaseWorkflowCaveat,
)
from phospy.contracts.results.kinase import (
    KinaseWorkflowPreprocessingAttritionSummary as KinaseWorkflowPreprocessingAttritionSummary,
)
from phospy.contracts.results.kinase import (
    KinaseWorkflowScoringAttritionSummary as KinaseWorkflowScoringAttritionSummary,
)
from phospy.contracts.results.kinase import (
    KinaseWorkflowSiteAttritionSummary as KinaseWorkflowSiteAttritionSummary,
)
from phospy.contracts.results.preprocessing import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
    ProteinAwareMappingDiagnostics,
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)
from phospy.contracts.results.signalome import SignalomeWorkflowResult

__all__ = [
    "ActivityMethodDiagnostics",
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "DifferentialAnalysisResult",
    "DifferentialContrastDefinition",
    "DifferentialDesignMatrixSummary",
    "DifferentialEmpiricalBayesProvenance",
    "DifferentialFixedEffectCovariateProvenance",
    "DifferentialMissingValuePolicyProvenance",
    "DifferentialModelDiagnostics",
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
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowAttritionProvenance",
    "KinaseWorkflowCaveat",
    "KinaseWorkflowResult",
    "KseaZScoreActivityDiagnostics",
    "PhosphositeImportResult",
    "ProteinAwareMappingDiagnostics",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
    "ProteinAwareSiteEligibility",
    "SsgseaSubstrateEnrichmentActivityDiagnostics",
    "SignalomeWorkflowResult",
    "WeightedSubstrateActivityDiagnostics",
]
