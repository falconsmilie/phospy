"""Public builder entrypoints."""

from __future__ import annotations

from phospy.io.bundles.reference_sources import ReferenceSourceTableReader
from phospy.io.readers.dataset_inputs import DatasetPathTableReader
from phospy.science.datasets.builders.contracts import (
    DatasetBuildExecutorContract,
    DatasetBuildInterpreterContract,
    DatasetBuildRequestProtocol,
    DatasetBuildValidatorContract,
)
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.science.datasets.builders.public import (
    AnalysisReadyDatasetBuilder as _AnalysisReadyDatasetBuilder,
)
from phospy.science.datasets.builders.reader import DatasetInputReader
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.preprocessing.stages import (
    SpsRuvStyleBatchCorrectionRunner,
)
from phospy.science.references.builder import (
    ReferenceBundleBuilder as _ReferenceBundleBuilder,
)
from phospy.validation.datasets.batch_correction import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
)
from phospy.validation.datasets.builder_request import DatasetBuildRequestValidator
from phospy.validation.datasets.inputs import DatasetInputSourceValidator
from phospy.validation.datasets.preprocessing import DatasetPreprocessingConfigValidator
from phospy.validation.references.builder import ReferenceBundleBuildRequestValidator
from phospy.workflows.batch_correction.preprocessing_adapter import (
    SpsRuvStyleBatchCorrectionWorkflowRunner,
)


class AnalysisReadyDatasetBuilder(_AnalysisReadyDatasetBuilder):
    """Supported public path with default local table reader wiring.

    The builder validates and interprets user inputs, establishes processing
    state, records construction provenance, and then constructs the strict
    analysis-ready dataset boundary.
    """

    def __init__(
        self,
        *,
        validator: DatasetBuildValidatorContract | None = None,
        interpreter: DatasetBuildInterpreterContract | None = None,
        executor: DatasetBuildExecutorContract | None = None,
        path_reader: DatasetPathTableReader | None = None,
        batch_correction_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
    ) -> None:
        if interpreter is None:
            interpreter = DatasetBuildRequestInterpreter(
                reader=DatasetInputReader(
                    source_validator=DatasetInputSourceValidator(),
                    path_reader=path_reader or DatasetPathTableReader(),
                )
            )
        if executor is None:
            executor = DatasetBuildExecutor(
                preprocessor=DatasetPreprocessor(
                    batch_correction_runner=(
                        batch_correction_runner
                        or SpsRuvStyleBatchCorrectionWorkflowRunner()
                    ),
                    batch_correction_metadata_validator=(
                        BatchDesignMetadataValidator()
                    ),
                    batch_correction_adequacy_validator=(
                        BatchCorrectionAdequacyValidator()
                    ),
                )
            )
        validator = validator or DatasetBuildRequestValidator(
            source_validator=DatasetInputSourceValidator(),
            preprocessing_validator=DatasetPreprocessingConfigValidator(),
        )
        super().__init__(
            validator=validator,
            interpreter=interpreter,
            executor=executor,
        )

    def run(
        self,
        request: DatasetBuildRequestProtocol,
    ) -> AnalysisReadyPhosphoDataset:
        """Validate, interpret, execute, and provenance-stamp a build request."""

        return super().run(request)


class ReferenceBundleBuilder(_ReferenceBundleBuilder):
    """Public reference-bundle builder with default local table reader wiring."""

    def __init__(
        self,
        *,
        source_reader: ReferenceSourceTableReader | None = None,
    ) -> None:
        super().__init__(
            source_reader=source_reader or ReferenceSourceTableReader(),
            request_validator=ReferenceBundleBuildRequestValidator(),
        )


__all__ = ["AnalysisReadyDatasetBuilder", "ReferenceBundleBuilder"]
