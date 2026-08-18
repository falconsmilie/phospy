# pyright: reportUnsupportedDunderAll=false
# ruff: noqa: F401
"""Advanced supported PhosPy API.

Use this namespace for specialized configuration, diagnostic result models,
reference-resource helpers, and table-inspection helpers that are supported but
outside the stable ``phospy.api`` compatibility tier.
"""

from __future__ import annotations

from phospy._api_inventory import ADVANCED_PUBLIC_API
from phospy.advanced.configs import (
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
    CorrectionMaskPolicy,
    CorrectionMissingnessPolicy,
    DatasetBatchCorrectionConfig,
    DatasetBatchCorrectionMethod,
    DatasetComparisonBuildingConfig,
    DatasetComparisonBuildingPolicy,
    DatasetComparisonPair,
    DatasetGroupCoverageFilterConfig,
    DatasetIntensityTransformConfig,
    DatasetIntensityTransformPolicy,
    DatasetLocalisationMode,
    DatasetMissingDataConfig,
    DatasetMissingDataInputScale,
    DatasetMissingDataKnnNoOverlapPolicy,
    DatasetMissingDataPolicy,
    DatasetNormalisationConfig,
    DatasetNormalisationPolicy,
    DatasetPreprocessingBatchCorrectionConfig,
    DatasetProteinAwarePreparationConfig,
    DatasetProteinAwarePreparationMappingPolicy,
    DatasetProteinAwarePreparationPolicy,
    DatasetRuvReadinessConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteMatrixDuplicateSitePolicy,
    DatasetSiteMatrixMissingDataPolicy,
    DatasetSiteMatrixPolicy,
    DatasetSiteSequenceConflictPolicy,
    DatasetSiteSequenceResolutionConfig,
    DatasetSiteSequenceResolutionMode,
    DatasetTotalProteinCorrectionConfig,
    DatasetTotalProteinCorrectionDuplicatePolicy,
    DatasetTotalProteinCorrectionIdentityConfig,
    DatasetTotalProteinCorrectionIdentityMatchingPolicy,
    DatasetTotalProteinCorrectionIdentityMode,
    DatasetTotalProteinCorrectionPolicy,
    DatasetTotalProteinCorrectionUnmatchedPolicy,
    DifferentialAnalysisConfig,
    DifferentialImputedValuePolicy,
    DifferentialReliabilityProfile,
    EmpiricalBayesConfig,
    EnrichmentIdentifierKind,
    EnrichmentMethod,
    EnrichmentOutsideBackgroundPolicy,
    KinaseActivityConfig,
    KinaseActivityMethod,
    KinaseActivityPValueMethod,
    KinaseActivitySsgseaRankingDirection,
    KinaseAdaptivePolicy,
    KinaseAttritionPolicy,
    KinaseAttritionViolationMode,
    KinasePredictionConfig,
    KinasePredictionMode,
    KinaseProfileMissingValueStrategy,
    KinaseReferenceDisplayAmbiguityPolicy,
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    KinaseScoringMode,
    KinaseSiteSequenceConflictPolicy,
    LocalisationPolicy,
    LocalisationRequirement,
    MultipleTestingConfig,
    MultipleTestingCorrection,
    MultipleTestingMethod,
    ObservationMask,
    OriginallyMissingCellTracking,
    PairedDesignPolicy,
    ProfileSelfInclusionPolicy,
    ReferenceContextCompatibilityPolicy,
    SignalomeAssignmentPolicy,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeKinaseNetworkPolicy,
    SignalomeMode,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeScorePreconditioningPolicy,
    SignalomeValidationConfig,
    SpsRuvBatchCorrectionConfig,
    SpsRuvBatchCorrectionMethod,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.advanced.results import (
    DifferentialModelDiagnostics,
    KinaseEligibilityReport,
    KinaseWorkflowAttritionProvenance,
    KinaseWorkflowCaveat,
    KinaseWorkflowPreprocessingAttritionSummary,
    KinaseWorkflowScoringAttritionSummary,
    KinaseWorkflowSiteAttritionSummary,
)
from phospy.io.bundles.kinase_library import (
    KinaseLibraryResourceLoader,
    load_kinase_library_resource,
)
from phospy.io.publishers.workflows import (
    publish_dataset,
    publish_kinase_workflow,
    publish_signalome_workflow,
)
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteAnnotation,
    ControlSiteSet,
    ControlSiteSourceMetadata,
    ControlSiteStatus,
)
from phospy.science.differential.policy_models import TechnicalReplicatePolicy
from phospy.science.references.kinase_library_models import (
    KinaseLibraryResource,
    KinaseLibraryResourceLoadRequest,
)
from phospy.science.tables.differential import (
    filter_differential_results,
    rank_differential_results,
)

__all__ = ADVANCED_PUBLIC_API
