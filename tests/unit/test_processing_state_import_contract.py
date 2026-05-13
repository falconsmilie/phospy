from __future__ import annotations

from phospy.science.datasets import processing_state
from phospy.science.datasets.processing_state import (
    DatasetProcessingState,
    JsonValue,
    MissingDataDiagnosticsV1,
    MissingDataState,
    TotalProteinCorrectionDiagnosticsV1,
)


def test_processing_state_public_import_contract() -> None:
    assert processing_state.DatasetProcessingState is DatasetProcessingState
    assert processing_state.MissingDataState is MissingDataState
    assert processing_state.MissingDataDiagnosticsV1 is MissingDataDiagnosticsV1
    assert (
        processing_state.TotalProteinCorrectionDiagnosticsV1
        is TotalProteinCorrectionDiagnosticsV1
    )
    assert isinstance(processing_state.MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1, int)
    assert isinstance(
        processing_state.TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
        int,
    )
    assert JsonValue is not None
