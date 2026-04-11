"""PhosPy public package root.

The root package keeps a deliberately small convenience surface for supported
high-level workflows, selected dataset and result types, and the main bundled
reference entry points. Implementation code is organised around domain
capability packages such as ``api``, ``datasets``, ``prediction``,
``preprocessing``, ``activities``, ``signalomes``, and ``references``. New
code should prefer those domain packages directly.
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

# Intentionally retained convenience exports for simple public usage.
# Keep this list small and prefer domain-package imports in new code.
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
