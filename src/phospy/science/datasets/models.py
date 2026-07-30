"""Compatibility facade for dataset domain model imports."""

from __future__ import annotations

from phospy.science.datasets.construction.analysis_ready import (
    AnalysisReadyPhosphoDataset as AnalysisReadyPhosphoDataset,
)
from phospy.science.datasets.construction.validation import (
    analysis_ready_matrix_missing_value_count as _analysis_ready_matrix_missing_value_count,
)
from phospy.science.datasets.imputation_metadata import (
    IMPUTATION_FEATURE_METADATA_COLUMNS as IMPUTATION_FEATURE_METADATA_COLUMNS,
)
from phospy.science.datasets.imputation_metadata import (
    IMPUTATION_OBSERVATION_SUMMARY_COLUMNS as IMPUTATION_OBSERVATION_SUMMARY_COLUMNS,
)
from phospy.science.datasets.imputation_metadata import (
    ImputationObservationMetadata as ImputationObservationMetadata,
)
from phospy.science.datasets.processing_state import (
    DatasetPreprocessingReport as DatasetPreprocessingReport,
)
from phospy.science.datasets.processing_state import (
    DatasetProcessingState as DatasetProcessingState,
)
from phospy.science.datasets.processing_state import (
    PreprocessingSiteAttritionSummary as PreprocessingSiteAttritionSummary,
)
from phospy.science.datasets.processing_state import (
    RuvReadinessState as RuvReadinessState,
)
from phospy.science.datasets.processing_state import (
    SiteSequenceResolutionReport as SiteSequenceResolutionReport,
)
from phospy.science.datasets.processing_state import (
    require_boolean_observation_mask as _require_boolean_observation_mask,
)

__all__ = [
    "AnalysisReadyPhosphoDataset",
    "DatasetPreprocessingReport",
    "DatasetProcessingState",
    "IMPUTATION_FEATURE_METADATA_COLUMNS",
    "IMPUTATION_OBSERVATION_SUMMARY_COLUMNS",
    "ImputationObservationMetadata",
    "PreprocessingSiteAttritionSummary",
    "RuvReadinessState",
    "SiteSequenceResolutionReport",
    "_analysis_ready_matrix_missing_value_count",
    "_require_boolean_observation_mask",
]
