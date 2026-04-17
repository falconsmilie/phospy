"""Public dataset builder shell."""

from __future__ import annotations

from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.builders.contracts import (
    DatasetBuildExecutorContract,
    DatasetBuildInterpreterContract,
    DatasetBuildValidatorContract,
)
from phospy.datasets.builders.executor import DatasetBuildExecutor
from phospy.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.datasets.builders.validator import DatasetBuildRequestValidator
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import PhosPyTransformationError
from phospy.errors.validation import PhosPyValidationError


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

        try:
            validated = self._validator.run(request)
            interpreted = self._interpreter.run(validated)
            return self._executor.run(interpreted)
        except (
            DatasetBuildError,
            PhosPyInputError,
            PhosPyTransformationError,
            PhosPyValidationError,
        ):
            raise
        except Exception as exc:  # pragma: no cover - defensive boundary translation
            raise DatasetBuildError(
                "failed to construct AnalysisReadyPhosphoDataset from build request"
            ) from exc
