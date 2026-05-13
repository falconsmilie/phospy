"""Dataset validators."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from phospy.validation.datasets.analysis_ready import (
        AnalysisReadyDatasetModelBoundaryValidator,
    )
    from phospy.validation.datasets.inputs import DatasetInputSourceValidator
    from phospy.validation.datasets.preprocessing import (
        DatasetPreprocessingConfigValidator,
    )

__all__ = [
    "AnalysisReadyDatasetModelBoundaryValidator",
    "DatasetInputSourceValidator",
    "DatasetPreprocessingConfigValidator",
]


def __getattr__(name: str) -> object:
    if name == "AnalysisReadyDatasetModelBoundaryValidator":
        from phospy.validation.datasets.analysis_ready import (
            AnalysisReadyDatasetModelBoundaryValidator,
        )

        return AnalysisReadyDatasetModelBoundaryValidator
    if name == "DatasetInputSourceValidator":
        from phospy.validation.datasets.inputs import DatasetInputSourceValidator

        return DatasetInputSourceValidator
    if name == "DatasetPreprocessingConfigValidator":
        from phospy.validation.datasets.preprocessing import (
            DatasetPreprocessingConfigValidator,
        )

        return DatasetPreprocessingConfigValidator
    raise AttributeError(name)
