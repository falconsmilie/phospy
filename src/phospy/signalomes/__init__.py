"""Signalome analysis domain.

This package owns signalome construction and related downstream outputs such as
map- and network-oriented result models. Public orchestration should stay out
of this package.
"""

from __future__ import annotations

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
    SignalomeKinaseNetwork,
    SignalomeModules,
    SignalomeResult,
)


def build_signalome_result(*args: object, **kwargs: object) -> SignalomeResult:
    """Build a structured signalome result from validated aligned inputs."""

    from .analysis import build_signalome_result as _build_signalome_result

    return _build_signalome_result(*args, **kwargs)


def execute_signalome_inputs(*args: object, **kwargs: object) -> SignalomeResult:
    """Build a signalome result from trusted signalome inputs."""

    from .analysis import execute_signalome_inputs as _execute_signalome_inputs

    return _execute_signalome_inputs(*args, **kwargs)


__all__ = [
    "ExpandedSignalome",
    "SignalomeAssignments",
    "SignalomeKinaseNetwork",
    "SignalomeMapData",
    "SignalomeModules",
    "SignalomeNetworkData",
    "SignalomeNetworkEdge",
    "SignalomeNetworkNode",
    "SignalomeResult",
    "build_signalome_map_data",
    "build_signalome_network_data",
    "build_signalome_result",
    "execute_signalome_inputs",
]
