"""Signalome analysis domain.

This package owns signalome construction and related downstream outputs such as
map- and network-oriented result models. Public orchestration should stay out
of this package.
"""

from __future__ import annotations

from .clustering import (
    SignalomeModuleSelectionDiagnostics,
    SignalomeModuleSelectionPolicy,
)
from .maps import SignalomeMapData, build_signalome_map_data
from .networks import (
    SignalomeNetworkData,
    SignalomeNetworkEdge,
    SignalomeNetworkNode,
    build_signalome_network_data,
)
from .results import (
    ExpandedSignalome,
    SignalomeAssignments,
    SignalomeCompatibilityView,
    SignalomeCoreResult,
    SignalomeExportAdapter,
    SignalomeFrameBundle,
    SignalomeKinaseNetwork,
    SignalomeModules,
    SignalomeResult,
    SignalomeVisualizationAdapter,
)


def build_signalome_result(*args: object, **kwargs: object) -> SignalomeResult:
    """Build a structured signalome result from validated aligned inputs."""

    from .analysis import build_signalome_result as _build_signalome_result

    return _build_signalome_result(*args, **kwargs)


__all__ = [
    "ExpandedSignalome",
    "SignalomeAssignments",
    "SignalomeCompatibilityView",
    "SignalomeCoreResult",
    "SignalomeExportAdapter",
    "SignalomeFrameBundle",
    "SignalomeKinaseNetwork",
    "SignalomeMapData",
    "SignalomeModuleSelectionDiagnostics",
    "SignalomeModuleSelectionPolicy",
    "SignalomeModules",
    "SignalomeNetworkData",
    "SignalomeNetworkEdge",
    "SignalomeNetworkNode",
    "SignalomeResult",
    "SignalomeVisualizationAdapter",
    "build_signalome_map_data",
    "build_signalome_network_data",
    "build_signalome_result",
]
