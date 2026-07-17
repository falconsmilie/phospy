from __future__ import annotations

from phospy.api.builders import AnalysisReadyDatasetBuilder, ReferenceBundleBuilder
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.validation.datasets.batch_correction import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
)
from phospy.validation.datasets.builder_request import DatasetBuildRequestValidator
from phospy.validation.datasets.inputs import DatasetInputSourceValidator
from phospy.validation.references.builder import ReferenceBundleBuildRequestValidator


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
    pipeline = preprocessor._pipeline
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
