from __future__ import annotations

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    ReferenceBundleBuilder,
)
from phospy.datasets.construction.executor import DatasetBuildExecutor
from phospy.datasets.construction.preprocessing import DatasetPreprocessor
from phospy.datasets.validation.inputs import DatasetInputSourceValidator
from phospy.datasets.validation.request import DatasetBuildRequestValidator
from phospy.preprocessing.batch_correction.validation.design import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
)
from phospy.references.builder_validation import ReferenceBundleBuildRequestValidator


def test_public_dataset_builder_wires_private_validation_adapters() -> None:
    builder = AnalysisReadyDatasetBuilder()

    assert isinstance(builder._validator, DatasetBuildRequestValidator)
    assert isinstance(
        builder._interpreter._source_resolver._reader._source_validator,
        DatasetInputSourceValidator,
    )
    assert isinstance(builder._executor, DatasetBuildExecutor)
    preprocessor = builder._executor._preprocessor
    assert isinstance(preprocessor, DatasetPreprocessor)
    pipeline = preprocessor._runner._pipeline
    batch_stage = pipeline._stages_by_key["batch_correction"]
    assert isinstance(
        batch_stage._metadata_validator,
        BatchDesignMetadataValidator,
    )
    assert isinstance(
        batch_stage._adequacy_validator,
        BatchCorrectionAdequacyValidator,
    )


def test_public_reference_builder_wires_private_request_validator() -> None:
    builder = ReferenceBundleBuilder()

    assert isinstance(
        builder._request_validator,
        ReferenceBundleBuildRequestValidator,
    )
