"""Dataset models, dataset construction, and dataset result containers.

This package owns dataset-shaped models, loaders, builders, and dataset result
containers. It does not own preprocessing strategy.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

from .builders import DatasetSiteMatrix
from .models import (
    AnalysisReadyPhosphoDataset,
    AnalysisReadyPreprocessingProvenance,
    AnalysisReadyRowCounts,
    AnalysisReadySiteMatrixStats,
    CoreInputs,
    PhosphoDataset,
)
from .schema import DatasetSchema

__all__ = [
    "AnalysisReadyPhosphoDataset",
    "AnalysisReadyPreprocessingProvenance",
    "AnalysisReadyRowCounts",
    "AnalysisReadySiteMatrixStats",
    "CoreInputs",
    "DatasetLoader",
    "DatasetSchema",
    "DatasetSiteMatrix",
    "LoadedDatasetInputs",
    "PhosphoDataset",
]


def __getattr__(name: str) -> Any:
    if name in {"DatasetLoader", "LoadedDatasetInputs"}:
        module = import_module(".loaders", __name__)
        return getattr(module, name)
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
