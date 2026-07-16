"""Bundle ownership for reproducible workflow serialization."""
# pyright: reportUnsupportedDunderAll=false

from __future__ import annotations

__all__ = [
    "KINASE_BUNDLE_MANIFEST_VERSION",
    "KinaseWorkflowConfigSnapshot",
    "LoadedKinaseWorkflowBundle",
    "LoadedSignalomeWorkflowBundle",
    "ReferenceSourceTableReader",
    "SIGNALOME_BUNDLE_MANIFEST_VERSION",
    "SignalomeWorkflowConfigSnapshot",
    "load_kinase_workflow_bundle",
    "load_signalome_workflow_bundle",
    "save_kinase_workflow_bundle",
    "save_signalome_workflow_bundle",
]


def __getattr__(name: str) -> object:
    if name == "ReferenceSourceTableReader":
        from phospy.io.bundles.reference_sources import ReferenceSourceTableReader

        return ReferenceSourceTableReader
    if name in {
        "KINASE_BUNDLE_MANIFEST_VERSION",
        "KinaseWorkflowConfigSnapshot",
        "LoadedKinaseWorkflowBundle",
        "load_kinase_workflow_bundle",
        "save_kinase_workflow_bundle",
    }:
        from phospy.io.bundles import kinase as _kinase

        return getattr(_kinase, name)
    if name in {
        "SIGNALOME_BUNDLE_MANIFEST_VERSION",
        "LoadedSignalomeWorkflowBundle",
        "SignalomeWorkflowConfigSnapshot",
        "load_signalome_workflow_bundle",
        "save_signalome_workflow_bundle",
    }:
        from phospy.io.bundles import signalome as _signalome

        return getattr(_signalome, name)
    raise AttributeError(name)
