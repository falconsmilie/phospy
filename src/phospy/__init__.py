"""Curated convenience entrypoints for the PhosPy package root."""

from phospy.api import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.workflows.differential.public import DifferentialAnalysisWorkflow
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.signalome.public import SignalomeWorkflow

__all__ = [
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DifferentialAnalysisWorkflow",
    "KinaseWorkflow",
    "SignalomeWorkflow",
]
