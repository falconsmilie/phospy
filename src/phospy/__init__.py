"""Curated convenience entrypoints for the PhosPy package root.

`phospy.api` is the authoritative public-contract namespace.
Top-level `phospy` intentionally exposes only the primary product entrypoints.
"""

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
    SignalomeWorkflow,
)

__all__ = [
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "KinaseWorkflow",
    "SignalomeWorkflow",
]
