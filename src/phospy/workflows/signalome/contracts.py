"""Signalome workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.api.configs import SignalomeConfig
from phospy.api.results import SimpleKinaseWorkflowResult


@dataclass(frozen=True, slots=True)
class ResolvedSignalomeWorkflowRequest:
    """Interpreter output for signalome workflow execution."""

    kinase_result: SimpleKinaseWorkflowResult
    config: SignalomeConfig
