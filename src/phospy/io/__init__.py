"""I/O adapters and CLI plumbing for the supported rewrite lane."""

from phospy.io.cli import main
from phospy.io.kinase_bundle import (
    KINASE_BUNDLE_MANIFEST_VERSION,
    KinaseWorkflowConfigSnapshot,
    LoadedKinaseWorkflowBundle,
    load_kinase_workflow_bundle,
    save_kinase_workflow_bundle,
)
from phospy.io.signalome_bundle import (
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
    LoadedSignalomeWorkflowBundle,
    SignalomeWorkflowConfigSnapshot,
    load_signalome_workflow_bundle,
    save_signalome_workflow_bundle,
)

__all__ = [
    "KINASE_BUNDLE_MANIFEST_VERSION",
    "KinaseWorkflowConfigSnapshot",
    "LoadedKinaseWorkflowBundle",
    "LoadedSignalomeWorkflowBundle",
    "SIGNALOME_BUNDLE_MANIFEST_VERSION",
    "SignalomeWorkflowConfigSnapshot",
    "load_kinase_workflow_bundle",
    "load_signalome_workflow_bundle",
    "main",
    "save_kinase_workflow_bundle",
    "save_signalome_workflow_bundle",
]
