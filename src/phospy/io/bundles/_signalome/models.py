"""Loaded signalome bundle DTOs."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.io.bundles._signalome.snapshots import SignalomeWorkflowConfigSnapshot
from phospy.signalome.results import SignalomeWorkflowResult


@dataclass(frozen=True, slots=True)
class LoadedSignalomeWorkflowBundle:
    """Loaded signalome output bundle contents."""

    result: SignalomeWorkflowResult
    config_snapshot: SignalomeWorkflowConfigSnapshot
    manifest_version: int
