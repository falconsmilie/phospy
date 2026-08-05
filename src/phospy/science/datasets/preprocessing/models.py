"""Compatibility facade for preprocessing planning, trace, and result models."""

from __future__ import annotations

from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION as DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS as DATASET_PREPROCESSING_STAGE_COMPARISONS,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER as DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM as DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_LOCALISATION as DATASET_PREPROCESSING_STAGE_LOCALISATION,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA as DATASET_PREPROCESSING_STAGE_MISSING_DATA,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_NORMALISATION as DATASET_PREPROCESSING_STAGE_NORMALISATION,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT as DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_PROTEIN_AWARE_PREPARATION as DATASET_PREPROCESSING_STAGE_PROTEIN_AWARE_PREPARATION,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX as DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION as DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
)
from phospy.science.datasets.preprocessing.plan import (
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION as DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
)
from phospy.science.datasets.preprocessing.plan import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_BATCH_CORRECTION as PREPROCESSING_STAGE_ORDER_RATIONALE_BATCH_CORRECTION,
)
from phospy.science.datasets.preprocessing.plan import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE as PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
)
from phospy.science.datasets.preprocessing.plan import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER as PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER,
)
from phospy.science.datasets.preprocessing.plan import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM as PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM,
)
from phospy.science.datasets.preprocessing.plan import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA as PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
)
from phospy.science.datasets.preprocessing.plan import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM as PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM,
)
from phospy.science.datasets.preprocessing.plan import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA as PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
)
from phospy.science.datasets.preprocessing.plan import (
    PreprocessingPlan as PreprocessingPlan,
)
from phospy.science.datasets.preprocessing.plan import (
    PreprocessingStageOrderResolution as PreprocessingStageOrderResolution,
)
from phospy.science.datasets.preprocessing.plan import (
    PreprocessingStageOrderValidator as PreprocessingStageOrderValidator,
)
from phospy.science.datasets.preprocessing.plan import (
    TotalProteinCorrectionIdentityPolicy as TotalProteinCorrectionIdentityPolicy,
)
from phospy.science.datasets.preprocessing.plan import (
    reject_external_corrected_output_after_downstream_preprocessing as reject_external_corrected_output_after_downstream_preprocessing,
)
from phospy.science.datasets.preprocessing.plan_interpreter import (
    PreprocessingPlanInterpreter as PreprocessingPlanInterpreter,
)
from phospy.science.datasets.preprocessing.results import (
    ComparisonBuildResult as ComparisonBuildResult,
)
from phospy.science.datasets.preprocessing.results import (
    DuplicateSiteResolutionResult as DuplicateSiteResolutionResult,
)
from phospy.science.datasets.preprocessing.results import (
    PreprocessingReportRow as PreprocessingReportRow,
)
from phospy.science.datasets.preprocessing.results import (
    PreprocessingStage as PreprocessingStage,
)
from phospy.science.datasets.preprocessing.results import (
    PreprocessingStageExecution as PreprocessingStageExecution,
)
from phospy.science.datasets.preprocessing.results import (
    PreprocessingStageResult as PreprocessingStageResult,
)
from phospy.science.datasets.preprocessing.results import (
    StageOwnedPreprocessingReportValue as StageOwnedPreprocessingReportValue,
)
from phospy.science.datasets.preprocessing.trace import (
    PREPROCESSING_STATE_TABLE_KEYS as PREPROCESSING_STATE_TABLE_KEYS,
)
from phospy.science.datasets.preprocessing.trace import (
    PreprocessingState as PreprocessingState,
)
from phospy.science.datasets.preprocessing.trace import (
    PreprocessingStateTableKey as PreprocessingStateTableKey,
)
from phospy.science.datasets.preprocessing.trace import (
    append_row_audit_records as append_row_audit_records,
)
from phospy.science.datasets.preprocessing.trace import (
    empty_preprocessing_row_audit as empty_preprocessing_row_audit,
)

__all__ = [
    "DATASET_PREPROCESSING_STAGE_COMPARISONS",
    "DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION",
    "DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER",
    "DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM",
    "DATASET_PREPROCESSING_STAGE_LOCALISATION",
    "DATASET_PREPROCESSING_STAGE_MISSING_DATA",
    "DATASET_PREPROCESSING_STAGE_NORMALISATION",
    "DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT",
    "DATASET_PREPROCESSING_STAGE_PROTEIN_AWARE_PREPARATION",
    "DATASET_PREPROCESSING_STAGE_SITE_MATRIX",
    "DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION",
    "DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_BATCH_CORRECTION",
    "PREPROCESSING_STATE_TABLE_KEYS",
    "ComparisonBuildResult",
    "DuplicateSiteResolutionResult",
    "append_row_audit_records",
    "empty_preprocessing_row_audit",
    "PreprocessingStateTableKey",
    "PreprocessingPlan",
    "PreprocessingPlanInterpreter",
    "PreprocessingStageOrderValidator",
    "PreprocessingStageOrderResolution",
    "PreprocessingReportRow",
    "PreprocessingStageResult",
    "PreprocessingStageExecution",
    "PreprocessingStage",
    "PreprocessingState",
    "StageOwnedPreprocessingReportValue",
    "TotalProteinCorrectionIdentityPolicy",
    "reject_external_corrected_output_after_downstream_preprocessing",
]
