"""Curated convenience entrypoints for the PhosPy package root."""

from importlib import metadata as _metadata

from phospy.api import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.workflows.differential.public import DifferentialAnalysisWorkflow
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.signalome.public import SignalomeWorkflow

try:
    __version__ = _metadata.version("phospy")
except _metadata.PackageNotFoundError:  # pragma: no cover - source tree fallback
    __version__ = "0+unknown"

__all__ = [
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DifferentialAnalysisWorkflow",
    "KinaseWorkflow",
    "SignalomeWorkflow",
]
