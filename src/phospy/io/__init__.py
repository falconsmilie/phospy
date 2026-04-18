"""I/O adapters and CLI plumbing for the supported rewrite lane."""

from phospy.io.cli import main
from phospy.io.signalome_bundle import (
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
    LoadedSignalomeWorkflowBundle,
    SignalomeWorkflowConfigSnapshot,
    load_signalome_workflow_bundle,
    save_signalome_workflow_bundle,
)
from phospy.io.simple_kinase_bundle import (
    SIMPLE_KINASE_BUNDLE_MANIFEST_VERSION,
    LoadedSimpleKinaseWorkflowBundle,
    SimpleKinaseWorkflowConfigSnapshot,
    load_simple_kinase_workflow_bundle,
    save_simple_kinase_workflow_bundle,
)

__all__ = [
    "LoadedSimpleKinaseWorkflowBundle",
    "LoadedSignalomeWorkflowBundle",
    "SIGNALOME_BUNDLE_MANIFEST_VERSION",
    "SIMPLE_KINASE_BUNDLE_MANIFEST_VERSION",
    "SignalomeWorkflowConfigSnapshot",
    "SimpleKinaseWorkflowConfigSnapshot",
    "load_signalome_workflow_bundle",
    "load_simple_kinase_workflow_bundle",
    "main",
    "save_signalome_workflow_bundle",
    "save_simple_kinase_workflow_bundle",
]
