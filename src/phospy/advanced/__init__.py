# ruff: noqa: F401
"""The advanced supported PhosPy API defined by ADR-0054."""

from __future__ import annotations

from phospy._api_inventory import ADVANCED_PUBLIC_API as __all__
from phospy.differential.config import (
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
    MultipleTestingConfig,
)
from phospy.differential.diagnostics import (
    DifferentialModelDiagnostics,
    QuantificationDepthTrendDiagnostics,
)
from phospy.differential.moderation.config import EmpiricalBayesConfig
from phospy.differential.policies import (
    DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1,
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
    DifferentialImputedValuePolicy,
    DifferentialReliabilityProfile,
    MultipleTestingMethod,
    PairedDesignPolicy,
    TechnicalReplicatePolicy,
)
from phospy.differential.tables import (
    filter_differential_results,
    rank_differential_results,
)
from phospy.enrichment.contracts.vocabulary import (
    EnrichmentIdentifierKind,
    EnrichmentMethod,
    EnrichmentOutsideBackgroundPolicy,
)
from phospy.evidence.config import LocalisationPolicy, LocalisationRequirement
from phospy.io.bundles.kinase_library import (
    KinaseLibraryResourceLoader,
    load_kinase_library_resource,
)
from phospy.io.publishers.workflows import (
    publish_dataset,
    publish_kinase_workflow,
    publish_signalome_workflow,
)
from phospy.kinase.activity.config import (
    KinaseActivityConfig,
    KinaseActivityMethod,
    KinaseActivityPValueMethod,
    KinaseActivitySsgseaRankingDirection,
)
from phospy.kinase.prediction.config import (
    KinaseAdaptivePolicy,
    KinasePredictionConfig,
    KinasePredictionMode,
)
from phospy.kinase.references.models import (
    KinaseLibraryResource,
    KinaseLibraryResourceLoadRequest,
)
from phospy.kinase.results import (
    KinaseEligibilityReport,
    KinaseWorkflowAttritionProvenance,
    KinaseWorkflowCaveat,
    KinaseWorkflowPreprocessingAttritionSummary,
    KinaseWorkflowScoringAttritionSummary,
    KinaseWorkflowSiteAttritionSummary,
)
from phospy.kinase.scoring.config import (
    KinaseAttritionPolicy,
    KinaseAttritionViolationMode,
    KinaseProfileMissingValueStrategy,
    KinaseReferenceDisplayAmbiguityPolicy,
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    KinaseScoringMode,
    KinaseSiteSequenceConflictPolicy,
)
from phospy.kinase.scoring.policies import ProfileSelfInclusionPolicy
from phospy.preprocessing.batch_correction.contracts.configuration import (
    DatasetBatchCorrectionConfig,
    DatasetBatchCorrectionMethod,
    DatasetPreprocessingBatchCorrectionConfig,
    SpsRuvBatchCorrectionConfig,
    SpsRuvBatchCorrectionMethod,
)
from phospy.preprocessing.batch_correction.contracts.missingness import (
    CorrectionMaskPolicy,
    CorrectionMissingnessPolicy,
    ObservationMask,
    OriginallyMissingCellTracking,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.preprocessing.contracts.configuration import (
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
)
from phospy.preprocessing.control_sites.contracts import (
    ControlSiteAnnotation,
    ControlSiteSet,
    ControlSiteSourceMetadata,
    ControlSiteStatus,
)
from phospy.preprocessing.control_sites.discovery_contracts import (
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryValidationError,
    SpsDiscoveryValidationIssue,
    SpsDiscoveryValidationResult,
    SpsReferenceDataset,
)
from phospy.preprocessing.control_sites.discovery_provenance import (
    SpsDiscoveryProvenance,
    SpsReferenceDatasetProvenance,
    SpsSampleConditionAssignment,
    SpsSelectionBoundaryCounts,
)
from phospy.preprocessing.control_sites.discovery_results import (
    SpsDatasetSiteStatistic,
    SpsDiscoveryResult,
    SpsSiteStabilityRecord,
)
from phospy.preprocessing.control_sites.workflow import (
    SpsDiscoveryWorkflow,
)
from phospy.references.config import ReferenceContextCompatibilityPolicy
from phospy.signalome.config import (
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
)
from phospy.statistics.multiple_testing import MultipleTestingCorrection
