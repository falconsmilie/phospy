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
    """Supported public path for building ``AnalysisReadyPhosphoDataset``.

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
    ) -> None:
        self._validator = validator or DatasetBuildRequestValidator()
        self._interpreter = interpreter or DatasetBuildRequestInterpreter()
        self._executor = executor or DatasetBuildExecutor()

    def run(self, request: DatasetBuildRequest) -> AnalysisReadyPhosphoDataset:
        """Validate, interpret, execute, and provenance-stamp a build request."""
        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return self._executor.run(interpreted)
