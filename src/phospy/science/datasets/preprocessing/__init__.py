"""Internal dataset preprocessing subsystem."""

from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionDiagnostics,
    BatchCorrectionEngine,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
    BatchCorrectionResult,
    LinearResidualizeBatchCorrectionEngine,
)
from phospy.science.datasets.preprocessing.batch_correction_metadata import (
    BatchCorrectionMetadataResolver,
    ResolvedBatchCorrectionMetadata,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStage,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    IntensityTransformPolicy,
    LocalisationEligibilityMode,
    MissingDataPolicy,
    NormalisationPolicy,
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
    TotalProteinCorrectionIdentityMatchingPolicy,
    TotalProteinCorrectionPolicy,
)

__all__ = [
    "DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION",
    "DATASET_PREPROCESSING_STAGE_COMPARISONS",
    "DATASET_PREPROCESSING_STAGE_MISSING_DATA",
    "DATASET_PREPROCESSING_STAGE_SITE_MATRIX",
    "DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION",
    "BatchCorrectionDiagnostics",
    "BatchCorrectionEngine",
    "BatchCorrectionMetadataResolver",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "BatchCorrectionResult",
    "ComparisonBuildingPolicy",
    "IntensityTransformPolicy",
    "LocalisationEligibilityMode",
    "LinearResidualizeBatchCorrectionEngine",
    "MissingDataPolicy",
    "NormalisationPolicy",
    "PreprocessingPipeline",
    "PreprocessingPlan",
    "ResolvedBatchCorrectionMetadata",
    "PreprocessingStageResult",
    "PreprocessingStage",
    "PreprocessingState",
    "SiteMatrixDuplicateSitePolicy",
    "SiteMatrixMissingDataPolicy",
    "SiteMatrixPolicy",
    "SiteSequenceConflictPolicy",
    "SiteSequenceResolutionMode",
    "TotalProteinCorrectionIdentityMatchingPolicy",
    "TotalProteinCorrectionPolicy",
]
