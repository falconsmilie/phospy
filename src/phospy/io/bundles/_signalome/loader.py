"""Signalome bundle loading orchestration."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from phospy.io.bundles._shared.integrity import file_entry_path, verify_bundle_integrity
from phospy.io.bundles._shared.json_files import read_json
from phospy.io.bundles._shared.primitives import require_mapping
from phospy.io.bundles._signalome.constants import MANIFEST_FILENAME
from phospy.io.bundles._signalome.manifest import parse_manifest
from phospy.io.bundles._signalome.models import LoadedSignalomeWorkflowBundle
from phospy.io.bundles._signalome.reconstruction import reconstruct_signalome_result
from phospy.io.bundles._signalome.snapshots import SignalomeWorkflowConfigSnapshot


def load_signalome_workflow_bundle(bundle_root: Path) -> LoadedSignalomeWorkflowBundle:
    """Load a signalome output bundle from disk."""

    root = Path(bundle_root)
    manifest_payload = require_mapping(
        read_json(root / MANIFEST_FILENAME, label="bundle manifest"),
        field_name="bundle manifest",
    )
    sections = parse_manifest(manifest_payload)
    verify_bundle_integrity(
        bundle_root=root,
        manifest_payload=manifest_payload,
        manifest_filename=MANIFEST_FILENAME,
    )
    result = reconstruct_signalome_result(bundle_root=root, sections=sections)

    config_snapshot_path = file_entry_path(
        sections.config_snapshot_entry,
        bundle_root=root,
        field_name="bundle manifest.config_snapshot",
    )
    config_snapshot_payload = require_mapping(
        read_json(config_snapshot_path, label="config snapshot"),
        field_name="config snapshot",
    )
    config_snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        config_snapshot_payload,
        effective_network_min_paired_finite_observations=(
            _effective_network_minimum_from_result_provenance(result.provenance)
        ),
    )

    return LoadedSignalomeWorkflowBundle(
        result=result,
        config_snapshot=config_snapshot,
        manifest_version=sections.manifest_version,
    )


def _effective_network_minimum_from_result_provenance(
    provenance: object,
) -> int | None:
    workflow_parameters = getattr(provenance, "workflow_parameters", None)
    if not isinstance(workflow_parameters, Mapping):
        return None
    signalome_config = workflow_parameters.get("signalome_config")
    if not isinstance(signalome_config, Mapping):
        return None
    output = signalome_config.get("output")
    if not isinstance(output, Mapping):
        return None
    value = output.get("network_min_paired_finite_observations")
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return int(value)
    return None
