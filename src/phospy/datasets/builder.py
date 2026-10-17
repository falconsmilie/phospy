"""Concrete Dataset-owned analysis-ready builder composition root."""

from __future__ import annotations

from typing import cast

from phospy.datasets.construction.contracts import (
    DatasetBuildExecutorContract,
    DatasetBuildInterpreterContract,
    DatasetBuildValidatorContract,
)
from phospy.datasets.construction.executor import DatasetBuildExecutor
from phospy.datasets.construction.input_resolution import DatasetInputReader
from phospy.datasets.construction.interpreter import DatasetBuildRequestInterpreter
from phospy.datasets.construction.preprocessing import DatasetPreprocessor
from phospy.datasets.contracts import DatasetBuildRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.datasets.validation.inputs import DatasetInputSourceValidator
from phospy.datasets.validation.request import DatasetBuildRequestValidator
from phospy.io.readers.dataset_inputs import DatasetPathTableReader
from phospy.preprocessing.dataset_construction import (
    SpsRuvStyleBatchCorrectionRunner,
    create_default_dataset_preprocessing_runner,
)


class AnalysisReadyDatasetBuilder:
    """Supported public path for Dataset construction and construction provenance."""

    def __init__(self) -> None:
        self._init_components(
            validator=None,
            interpreter=None,
            executor=None,
            path_reader=None,
            batch_correction_runner=None,
        )

    @classmethod
    def _with_components(
        cls,
        *,
        validator: DatasetBuildValidatorContract | None = None,
        interpreter: DatasetBuildInterpreterContract | None = None,
        executor: DatasetBuildExecutorContract | None = None,
        path_reader: DatasetPathTableReader | None = None,
        batch_correction_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
    ) -> AnalysisReadyDatasetBuilder:
        builder = cls.__new__(cls)
        builder._init_components(
            validator=validator,
            interpreter=interpreter,
            executor=executor,
            path_reader=path_reader,
            batch_correction_runner=batch_correction_runner,
        )
        return builder

    def _init_components(
        self,
        *,
        validator: DatasetBuildValidatorContract | None,
        interpreter: DatasetBuildInterpreterContract | None,
        executor: DatasetBuildExecutorContract | None,
        path_reader: DatasetPathTableReader | None,
        batch_correction_runner: SpsRuvStyleBatchCorrectionRunner | None,
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
                    runner=create_default_dataset_preprocessing_runner(
                        batch_correction_runner=batch_correction_runner
                    )
                )
            )
        self._validator = validator or DatasetBuildRequestValidator(
            source_validator=DatasetInputSourceValidator()
        )
        self._interpreter = interpreter
        self._executor = executor

    def run(self, request: DatasetBuildRequest) -> AnalysisReadyPhosphoDataset:
        """Validate, interpret, execute, and provenance-stamp a build request."""

        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return cast(AnalysisReadyPhosphoDataset, self._executor.run(interpreted))


__all__ = ["AnalysisReadyDatasetBuilder"]
