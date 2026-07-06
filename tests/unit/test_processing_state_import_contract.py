from __future__ import annotations

from phospy.science.datasets import imputation_metadata, processing_state
from phospy.science.datasets.imputation_metadata import (
    ImputationObservationMetadata as DedicatedImputationObservationMetadata,
)
from phospy.science.datasets.models import (
    DatasetPreprocessingReport as ModelsDatasetPreprocessingReport,
)
from phospy.science.datasets.models import (
    ImputationObservationMetadata as ModelsImputationObservationMetadata,
)
from phospy.science.datasets.models import (
    PreprocessingSiteAttritionSummary as ModelsPreprocessingSiteAttritionSummary,
)
from phospy.science.datasets.models import (
    SiteSequenceResolutionReport as ModelsSiteSequenceResolutionReport,
)
from phospy.science.datasets.models import (
    _require_boolean_observation_mask as models_require_boolean_observation_mask,
)
from phospy.science.datasets.processing_state import (
    DatasetPreprocessingReport,
    DatasetProcessingState,
    ImputationObservationMetadata,
    JsonValue,
    MissingDataDiagnosticsV1,
    MissingDataState,
    PreprocessingSiteAttritionSummary,
    SiteSequenceResolutionReport,
    TotalProteinCorrectionDiagnosticsV1,
    require_boolean_observation_mask,
)


def test_processing_state_public_import_contract() -> None:
    assert processing_state.DatasetProcessingState is DatasetProcessingState
    assert processing_state.MissingDataState is MissingDataState
    assert processing_state.MissingDataDiagnosticsV1 is MissingDataDiagnosticsV1
    assert (
        processing_state.TotalProteinCorrectionDiagnosticsV1
        is TotalProteinCorrectionDiagnosticsV1
    )
    assert processing_state.DatasetPreprocessingReport is DatasetPreprocessingReport
    assert (
        processing_state.ImputationObservationMetadata is ImputationObservationMetadata
    )
    assert (
        processing_state.PreprocessingSiteAttritionSummary
        is PreprocessingSiteAttritionSummary
    )
    assert processing_state.SiteSequenceResolutionReport is SiteSequenceResolutionReport
    assert isinstance(processing_state.MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1, int)
    assert isinstance(
        processing_state.TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
        int,
    )
    assert JsonValue is not None
    assert (
        imputation_metadata.ImputationObservationMetadata
        is DedicatedImputationObservationMetadata
    )
    assert DedicatedImputationObservationMetadata is ImputationObservationMetadata


def test_processing_state_models_import_compatibility() -> None:
    assert ModelsDatasetPreprocessingReport is DatasetPreprocessingReport
    assert ModelsImputationObservationMetadata is ImputationObservationMetadata
    assert ModelsPreprocessingSiteAttritionSummary is PreprocessingSiteAttritionSummary
    assert ModelsSiteSequenceResolutionReport is SiteSequenceResolutionReport
    assert models_require_boolean_observation_mask is require_boolean_observation_mask
