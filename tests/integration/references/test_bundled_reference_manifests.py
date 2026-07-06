from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import pytest

from phospy.science.references import resources as reference_resources
from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.models import (
    Organism,
    RedistributionStatus,
    ReferencePreset,
)
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.references.resources import load_bundled_reference_manifest
from phospy.science.references.validation import validate_bundled_reference_manifests

pytestmark = pytest.mark.integration


def test_bundled_reference_manifests_are_structurally_valid() -> None:
    manifests = validate_bundled_reference_manifests(
        _reference_bundles_root(),
        require_redistribution_allowed=False,
    )

    assert {manifest.reference_id for manifest in manifests} == {"l6_native"}
    rat_manifest = manifests[0]
    assert rat_manifest.files
    assert {item.relative_path for item in rat_manifest.files} == {
        "motif_scores.csv",
        "motif_sizes.csv",
        "site_sequences.csv",
        "substrate_map.csv",
    }
    assert rat_manifest.redistribution_status is RedistributionStatus.APPROVED
    assert rat_manifest.redistribution_allowed is True


def test_runtime_bundled_manifest_validates_hashes_before_table_loading() -> None:
    reference_resources.clear_bundled_reference_manifest_cache()

    manifest = load_bundled_reference_manifest(Organism.RAT)

    assert manifest.reference_id == "l6_native"
    assert manifest.reference_version == "bundled-snapshot-2026-04-16"
    assert all(file.sha256 for file in manifest.files)


def test_runtime_reference_loader_refuses_invalid_manifest_before_workflow_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _write_invalid_bundle(tmp_path)
    reference_resources.clear_bundled_reference_manifest_cache()

    def _patched_bundle_resource(**_: object) -> Path:
        return tmp_path

    monkeypatch.setattr(
        reference_resources,
        "_bundled_reference_bundle_resource",
        _patched_bundle_resource,
    )

    with pytest.raises(ReferenceManifestError, match="hash mismatch"):
        ReferenceResolver().run(ReferencePreset.RAT, dataset_organism=Organism.RAT)


def _reference_bundles_root() -> Path:
    return (
        Path(__file__).resolve().parents[3]
        / "src"
        / "phospy"
        / "data"
        / "reference_bundles"
    )


def _write_invalid_bundle(tmp_path: Path) -> None:
    data = "kinase,site_id\nAKT1,MAPK1;S123;\n"
    (tmp_path / "substrate_map.csv").write_text(data, encoding="utf-8")
    payload = {
        "reference_id": "invalid_runtime_reference",
        "display_name": "Invalid runtime reference",
        "organism": "Rattus norvegicus",
        "organism_common_name": "rat",
        "taxonomy_id": 10116,
        "protein_namespace": "display_site_id",
        "reference_version": "v1",
        "source_name": "unit source",
        "source_version": "source-v1",
        "source_url": "https://example.test/reference",
        "retrieved_at": "2026-06-29",
        "table_sha256": sha256(b"different").hexdigest(),
        "source_publication": None,
        "license_name": "CC0 synthetic",
        "license_url": "https://example.test/license",
        "redistribution_status": "approved",
        "redistribution_notes": "synthetic test fixture",
        "derived_from": ["unit test"],
        "generated_by": "unit test",
        "generated_at_utc": "2026-06-29T00:00:00Z",
        "manifest_schema_version": "1.1",
        "files": [
            {
                "relative_path": "substrate_map.csv",
                "role": "kinase_substrate",
                "format": "csv",
                "sha256": sha256(b"different").hexdigest(),
                "row_count": 1,
                "column_names": ["kinase", "site_id"],
            }
        ],
        "sequence_context_policy": "centered phosphosite sequence window",
        "sequence_window_length": 3,
        "sequence_center_index": 1,
        "allowed_sequence_alphabet": "ACDEFGHIKLMNPQRSTVWY",
    }
    (tmp_path / "manifest.json").write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
