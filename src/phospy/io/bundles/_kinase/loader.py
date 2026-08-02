"""Kinase bundle loading orchestration."""

from __future__ import annotations

from pathlib import Path

from phospy.io.bundles._kinase.constants import MANIFEST_FILENAME
from phospy.io.bundles._kinase.manifest import parse_manifest
from phospy.io.bundles._kinase.models import LoadedKinaseWorkflowBundle
from phospy.io.bundles._kinase.reconstruction import reconstruct_kinase_result
from phospy.io.bundles._kinase.snapshots import KinaseWorkflowConfigSnapshot
from phospy.io.bundles._shared.integrity import file_entry_path, verify_bundle_integrity
from phospy.io.bundles._shared.json_files import read_json
from phospy.io.bundles._shared.primitives import require_mapping


def load_kinase_workflow_bundle(bundle_root: Path) -> LoadedKinaseWorkflowBundle:
    """Load a kinase output bundle from disk."""

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
    result = reconstruct_kinase_result(bundle_root=root, sections=sections)

    config_snapshot_path = file_entry_path(
        sections.config_snapshot_entry,
        bundle_root=root,
        field_name="bundle manifest.config_snapshot",
    )
    config_snapshot_payload = require_mapping(
        read_json(config_snapshot_path, label="config snapshot"),
        field_name="config snapshot",
    )
    config_snapshot = KinaseWorkflowConfigSnapshot.from_payload(config_snapshot_payload)

    return LoadedKinaseWorkflowBundle(
        result=result,
        config_snapshot=config_snapshot,
        manifest_version=sections.manifest_version,
    )
