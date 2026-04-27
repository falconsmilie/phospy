"""Compatibility exports for signalome workflow components.

This module remains as a thin import shim while focused component implementations
live in dedicated modules.
"""

from __future__ import annotations

from phospy.workflows.signalome.clustering_runner import SignalomeClusteringRunner
from phospy.workflows.signalome.component_models import (
    SignalomeClusteringRunResult,
    SignalomeContextTableBuildResult,
    SignalomeExecutionMetadata,
    SignalomeModuleTableBuildResult,
    SignalomeNetworkBuildResult,
    SignalomeScaleGuardDecision,
    SignalomeSupportSummary,
)
from phospy.workflows.signalome.context_tables import SignalomeContextTableBuilder
from phospy.workflows.signalome.module_tables import SignalomeModuleTableBuilder
from phospy.workflows.signalome.network_builder import SignalomeNetworkBuilder
from phospy.workflows.signalome.provenance import SignalomeProvenanceBuilder
from phospy.workflows.signalome.result_assembly import SignalomeResultAssembler

__all__ = [
    "SignalomeClusteringRunResult",
    "SignalomeClusteringRunner",
    "SignalomeContextTableBuildResult",
    "SignalomeContextTableBuilder",
    "SignalomeExecutionMetadata",
    "SignalomeModuleTableBuildResult",
    "SignalomeModuleTableBuilder",
    "SignalomeNetworkBuildResult",
    "SignalomeNetworkBuilder",
    "SignalomeProvenanceBuilder",
    "SignalomeResultAssembler",
    "SignalomeScaleGuardDecision",
    "SignalomeSupportSummary",
]
