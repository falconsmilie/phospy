"""Loaded signalome bundle DTOs."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.api.results import SignalomeWorkflowResult
from phospy.io.bundles._signalome.snapshots import SignalomeWorkflowConfigSnapshot


@dataclass(frozen=True, slots=True)
class LoadedSignalomeWorkflowBundle:
    """Loaded signalome output bundle contents."""

    result: SignalomeWorkflowResult
    config_snapshot: SignalomeWorkflowConfigSnapshot
    manifest_version: int
