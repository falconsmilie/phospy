from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import pytest

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.validation import validate_bundled_reference_manifests

pytestmark = pytest.mark.release_gate


def test_bundled_reference_manifests_are_release_eligible() -> None:
    validate_bundled_reference_manifests(_reference_bundles_root())


def test_reference_manifest_release_gate_fails_for_invalid_fixture(
    tmp_path: Path,
) -> None:
    bundle_root = tmp_path / "reference_bundles" / "rat" / "invalid"
    bundle_root.mkdir(parents=True)
    data = "kinase,site_id\nAKT1,MAPK1;S123;\n"
    (bundle_root / "substrate_map.csv").write_text(data, encoding="utf-8")
    payload = _valid_manifest_payload(sha256(data.encode("utf-8")).hexdigest())
    payload["license_name"] = ""
    (bundle_root / "manifest.json").write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )

    with pytest.raises(ReferenceManifestError, match="license_name"):
        validate_bundled_reference_manifests(tmp_path / "reference_bundles")


def _reference_bundles_root() -> Path:
    return (
        Path(__file__).resolve().parents[2]
        / "src"
        / "phospy"
        / "data"
        / "reference_bundles"
    )


def _valid_manifest_payload(file_hash: str) -> dict[str, object]:
    return {
        "reference_id": "invalid",
        "display_name": "Invalid",
        "organism": "Rattus norvegicus",
        "taxonomy_id": 10116,
        "protein_namespace": "display_site_id",
        "reference_version": "v1",
        "source_name": "unit source",
        "source_version": None,
        "source_url": "https://example.test/reference",
        "retrieved_at": "2026-06-29",
        "table_sha256": file_hash,
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
                "sha256": file_hash,
                "row_count": 1,
                "column_names": ["kinase", "site_id"],
            }
        ],
        "sequence_context_policy": "centered phosphosite sequence window",
        "sequence_window_length": 3,
        "sequence_center_index": 1,
        "allowed_sequence_alphabet": "ACDEFGHIKLMNPQRSTVWY",
    }
