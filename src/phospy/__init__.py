"""PhosPy public package root.

The root package exposes the supported public surface while implementation code
is organised around domain capability packages such as ``api``, ``datasets``,
``prediction``, ``preprocessing``, ``activities``, ``signalomes``, and
``references``.
"""

from .activities import KinaseActivityAnalyzer
from .api import (
    KinaseWorkflow,
    PredMatWorkflow,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from .datasets import AnalysisReadyPhosphoDataset, PhosphoDataset
from .pipeline import PhosRPipeline
from .prediction import PredMatResult
from .references import (
    BundledReferenceProvider,
    ReferenceBundle,
    ReferenceBundleProvenance,
    ReferenceBundleSourceMetadata,
    ReferenceProvider,
)
from .signalome_maps import SignalomeMapData
from .signalome_networks import SignalomeNetworkData
from .signalomes import SignalomeResult

__all__ = [
    "AnalysisReadyPhosphoDataset",
    "BundledReferenceProvider",
    "KinaseActivityAnalyzer",
    "ReferenceBundle",
    "ReferenceBundleProvenance",
    "ReferenceBundleSourceMetadata",
    "ReferenceProvider",
    "KinaseWorkflow",
    "PhosphoDataset",
    "PhosRPipeline",
    "PredMatResult",
    "PredMatWorkflow",
    "SignalomeMapData",
    "SignalomeNetworkData",
    "SignalomeResult",
    "SignalomeWorkflow",
    "SimpleKinaseWorkflow",
]
