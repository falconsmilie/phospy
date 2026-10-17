"""Loaded kinase bundle DTOs."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.io.bundles._kinase.snapshots import KinaseWorkflowConfigSnapshot
from phospy.kinase.results import KinaseWorkflowResult


@dataclass(frozen=True, slots=True)
class LoadedKinaseWorkflowBundle:
    """Loaded kinase output bundle contents."""

    result: KinaseWorkflowResult
    config_snapshot: KinaseWorkflowConfigSnapshot
    manifest_version: int
