"""Public dataset builder shell."""

from __future__ import annotations

from phospy.contracts.requests import DatasetBuildRequest
from phospy.science.datasets.builders.contracts import (
    DatasetBuildExecutorContract,
    DatasetBuildInterpreterContract,
    DatasetBuildValidatorContract,
)
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.science.datasets.builders.validator import DatasetBuildRequestValidator
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset


class AnalysisReadyDatasetBuilder:
    """Build ``AnalysisReadyPhosphoDataset`` from a structured request."""

    def __init__(
        self,
        *,
        validator: DatasetBuildValidatorContract | None = None,
        interpreter: DatasetBuildInterpreterContract | None = None,
        executor: DatasetBuildExecutorContract | None = None,
    ) -> None:
        self._validator = validator or DatasetBuildRequestValidator()
        self._interpreter = interpreter or DatasetBuildRequestInterpreter()
        self._executor = executor or DatasetBuildExecutor()

    def run(self, request: DatasetBuildRequest) -> AnalysisReadyPhosphoDataset:
        """Validate, interpret, and execute the dataset build request."""
        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return self._executor.run(interpreted)
