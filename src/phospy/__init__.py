"""PhosPy public package root.

The root package keeps a deliberately small convenience surface for supported
high-level types while implementation code is organised around domain capability
packages such as ``api``, ``datasets``, ``prediction``, ``preprocessing``,
``activities``, ``signalomes``, and ``references``. New code should prefer
those domain packages directly.
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
from .signalomes import SignalomeResult
from .signalomes.maps import SignalomeMapData
from .signalomes.networks import SignalomeNetworkData

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
