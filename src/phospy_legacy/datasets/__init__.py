"""Dataset models, dataset construction, and dataset result containers.

This package owns dataset-shaped models, loaders, builders, and dataset result
containers. It does not own preprocessing strategy.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

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
    "SiteToProteinResolutionDiagnostics",
    "SiteToProteinResolutionResult",
]


_EXPORT_MODULES: dict[str, tuple[str, str]] = {
    "AnalysisReadyPhosphoDataset": (".models", "AnalysisReadyPhosphoDataset"),
    "AnalysisReadyPreprocessingProvenance": (
        ".models",
        "AnalysisReadyPreprocessingProvenance",
    ),
    "AnalysisReadyRowCounts": (".models", "AnalysisReadyRowCounts"),
    "AnalysisReadySiteMatrixStats": (".models", "AnalysisReadySiteMatrixStats"),
    "CoreInputs": (".models", "CoreInputs"),
    "DatasetLoader": (".loaders", "DatasetLoader"),
    "DatasetSchema": (".schema", "DatasetSchema"),
    "DatasetSiteMatrix": (".builders", "DatasetSiteMatrix"),
    "LoadedDatasetInputs": (".loaders", "LoadedDatasetInputs"),
    "PhosphoDataset": (".models", "PhosphoDataset"),
    "SiteToProteinResolutionDiagnostics": (
        ".models",
        "SiteToProteinResolutionDiagnostics",
    ),
    "SiteToProteinResolutionResult": (".models", "SiteToProteinResolutionResult"),
}


def __getattr__(name: str) -> Any:
    module_export = _EXPORT_MODULES.get(name)
    if module_export is not None:
        module_name, export_name = module_export
        module = import_module(module_name, __name__)
        value = getattr(module, export_name)
        globals()[name] = value
        return value
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
