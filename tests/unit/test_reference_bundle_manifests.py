from __future__ import annotations

import json
from pathlib import Path

import pytest

from phospy.errors.references import (
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    UnsupportedOrganismError,
)
from phospy.science.references import resources as reference_resources
from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.models import Organism, ReferencePreset
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.references.resources import (
    available_bundled_reference_lanes,
    bundled_reference_name_for_organism,
    load_bundled_kinase_substrate_map,
    load_bundled_reference_manifest,
    load_bundled_site_sequences,
    supported_bundled_organisms,
)

_REQUIRED_RUNTIME_FILES = (
    "manifest.json",
    "substrate_map.csv",
    "site_sequences.csv",
)
_APPROVAL_REQUIRED_ORGANISMS = (Organism.HUMAN, Organism.MOUSE)


def _unsupported_approval_required_organism() -> Organism:
    supported = set(supported_bundled_organisms())
    for organism in _APPROVAL_REQUIRED_ORGANISMS:
        if organism not in supported:
            return organism
    pytest.skip("all human/mouse bundled lanes are committed in this test run")


def test_valid_bundled_manifest_loads_for_supported_runtime_lane() -> None:
    manifest = load_bundled_reference_manifest(Organism.RAT)

    assert manifest.reference_id == "l6_native"
    assert manifest.bundle_id == "l6_native"
    assert manifest.display_name
    assert manifest.organism == "Rattus norvegicus"
    assert manifest.organism_common_name == Organism.RAT.value
    assert manifest.taxonomy_id == 10116
    assert manifest.protein_namespace
    assert manifest.identifier_namespace == manifest.protein_namespace
    assert manifest.reference_version == "bundled-snapshot-2026-04-16"
    assert manifest.source_name
    assert manifest.source_url == "https://github.com/PYangLab/PhosR"
    assert manifest.source_license_url
    assert manifest.source_license
    assert "gpl-3" in manifest.source_license.lower()
    assert "phosphositeplus" in manifest.source_license.lower()
    assert "pride" in manifest.source_license.lower()
    assert manifest.redistribution_allowed is True
    assert "not independently verified" in manifest.redistribution_notes.lower()
    assert {item.relative_path for item in manifest.files} == {
        "motif_scores.csv",
        "motif_sizes.csv",
        "site_sequences.csv",
        "substrate_map.csv",
    }
    assert manifest.sequence_window.upstream_residues == 15
    assert manifest.sequence_window.downstream_residues == 15
    assert "site_sequence_derivation" in manifest.supports
    assert manifest.limitations
    assert any("not independently" in item.lower() for item in manifest.limitations)


def test_available_bundled_reference_lanes_reports_manifest_metadata() -> None:
    lanes = available_bundled_reference_lanes()

    rat_lanes = [lane for lane in lanes if lane.organism is Organism.RAT]
    assert len(rat_lanes) == 1
    lane = rat_lanes[0]
    assert lane.organism is Organism.RAT
    assert lane.bundle_id == "l6_native"
    assert lane.source_name
    assert lane.source_version == "bundled-snapshot-2026-04-16"
    assert lane.retrieved_at.isoformat() == "2026-04-16"
    assert lane.redistribution_status
    assert "not independently verified" in lane.redistribution_status.lower()
    assert "site_sequence_derivation" in lane.supports
    assert lane.limitations
    assert lane.to_payload()["organism"] == Organism.RAT.value


def test_supported_bundled_lanes_have_required_runtime_files_and_manifest_metadata() -> (
    None
):
    package_root = Path(__file__).resolve().parents[2] / "src" / "phospy"
    for organism in supported_bundled_organisms():
        bundle_id = bundled_reference_name_for_organism(organism)
        lane_root = (
            package_root / "data" / "reference_bundles" / organism.value / bundle_id
        )
        for filename in _REQUIRED_RUNTIME_FILES:
            assert lane_root.joinpath(filename).is_file(), (
                f"missing bundled reference runtime file: "
                f"{organism.value}/{bundle_id}/{filename}"
            )

        manifest = load_bundled_reference_manifest(organism)
        assert manifest.reference_id == bundle_id
        assert manifest.organism
        assert manifest.organism_common_name == organism.value
        assert manifest.protein_namespace
        assert manifest.source_name
        assert manifest.reference_version
        assert manifest.retrieved_at.isoformat()
        assert manifest.source_license
        assert isinstance(manifest.redistribution_allowed, bool)
        assert manifest.files
        assert manifest.sequence_window.upstream_residues >= 0
        assert manifest.sequence_window.downstream_residues >= 0


@pytest.mark.parametrize(
    ("organism", "preset"),
    [
        (Organism.HUMAN, ReferencePreset.HUMAN),
        (Organism.MOUSE, ReferencePreset.MOUSE),
    ],
)
def test_human_mouse_lane_loads_only_if_committed_and_approved(
    organism: Organism,
    preset: ReferencePreset,
) -> None:
    supported = supported_bundled_organisms()
    if organism not in supported:
        with pytest.raises(UnsupportedOrganismError):
            ReferenceResolver().run(preset, dataset_organism=organism)
        return

    bundle = ReferenceResolver().run(preset, dataset_organism=organism)

    assert bundle.organism is organism
    assert bundle.manifest is not None
    assert bundle.manifest.redistribution_allowed is True
    assert not bundle.kinase_substrate_map.empty
    assert not bundle.site_sequences.empty


@pytest.mark.parametrize("organism", _APPROVAL_REQUIRED_ORGANISMS)
def test_human_mouse_bundle_directories_are_not_stray_or_unapproved(
    organism: Organism,
) -> None:
    organism_root = (
        Path(__file__).resolve().parents[2]
        / "src"
        / "phospy"
        / "data"
        / "reference_bundles"
        / organism.value
    )
    if not organism_root.is_dir():
        assert organism not in supported_bundled_organisms()
        return

    default_bundle = reference_resources._BUNDLED_DEFAULTS.get(organism)
    assert default_bundle is not None, (
        f"{organism.value} bundle data is packaged but not declared in "
        "_BUNDLED_DEFAULTS"
    )
    committed_lanes = sorted(
        lane.name for lane in organism_root.iterdir() if lane.is_dir()
    )
    assert committed_lanes == [default_bundle]

    manifest = load_bundled_reference_manifest(organism)
    assert manifest.redistribution_allowed is True


def test_bundled_manifest_missing_resource_fails_clearly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()

    def _missing_bundle_resource(**_: object) -> object:
        return Path("missing-reference-bundle")

    monkeypatch.setattr(
        reference_resources,
        "_bundled_reference_bundle_resource",
        _missing_bundle_resource,
    )
    with pytest.raises(ReferenceResolutionError, match="directory is missing"):
        load_bundled_reference_manifest(Organism.RAT)


def test_bundled_manifest_hash_mismatch_fails_clearly(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _write_temporary_manifest_bundle(tmp_path, file_hash="0" * 64)
    reference_resources.clear_bundled_reference_manifest_cache()

    def _patched_bundle_resource(**_: object) -> Path:
        return tmp_path

    monkeypatch.setattr(
        reference_resources,
        "_bundled_reference_bundle_resource",
        _patched_bundle_resource,
    )

    with pytest.raises(ReferenceManifestError, match="hash mismatch"):
        load_bundled_reference_manifest(Organism.RAT)


def test_reference_resolution_succeeds_for_supported_runtime_organism() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )

    assert bundle.organism is Organism.RAT
    assert bundle.manifest is not None
    assert bundle.manifest.reference_id == "l6_native"


def test_reference_resolution_fails_for_unsupported_runtime_organism() -> None:
    organism = _unsupported_approval_required_organism()
    preset = (
        ReferencePreset.HUMAN if organism is Organism.HUMAN else ReferencePreset.MOUSE
    )
    with pytest.raises(
        UnsupportedOrganismError,
        match="reference-bundle docs: https://phospy.com/docs/api/guide/#references",
    ) as exc_info:
        ReferenceResolver().run(
            preset,
            dataset_organism=organism,
        )
    message = str(exc_info.value)
    assert "supported bundled organisms:" in message
    assert "ReferenceBundle(organism=..." in message


def test_auto_reference_resolution_succeeds_for_supported_runtime_organism() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.AUTO,
        dataset_organism=Organism.RAT,
    )

    assert bundle.manifest is not None
    assert bundle.manifest.organism_common_name == Organism.RAT.value


@pytest.mark.parametrize("organism", _APPROVAL_REQUIRED_ORGANISMS)
def test_auto_reference_resolution_for_human_mouse_matches_committed_lane_presence(
    organism: Organism,
) -> None:
    if organism not in supported_bundled_organisms():
        with pytest.raises(UnsupportedOrganismError):
            ReferenceResolver().run(ReferencePreset.AUTO, dataset_organism=organism)
        return

    bundle = ReferenceResolver().run(ReferencePreset.AUTO, dataset_organism=organism)

    assert bundle.organism is organism
    assert bundle.manifest is not None


def test_auto_reference_resolution_fails_for_unsupported_runtime_organism() -> None:
    organism = _unsupported_approval_required_organism()
    with pytest.raises(
        UnsupportedOrganismError,
        match="reference-bundle docs: https://phospy.com/docs/api/guide/#references",
    ) as exc_info:
        ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=organism,
        )
    assert "supported bundled organisms:" in str(exc_info.value)


def test_wrong_organism_fails_compatibility_validation_before_bundle_loading() -> None:
    with pytest.raises(
        ReferenceCompatibilityError,
        match="dataset.organism and requested reference preset must match",
    ):
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.MOUSE,
        )


def test_bundled_human_mouse_tables_load_if_present() -> None:
    for organism in _APPROVAL_REQUIRED_ORGANISMS:
        if organism not in supported_bundled_organisms():
            continue

        substrate_map = load_bundled_kinase_substrate_map(organism)
        site_sequences = load_bundled_site_sequences(organism)

        assert not substrate_map.empty
        assert not site_sequences.empty


def test_reference_provenance_includes_manifest_fields_for_bundled_reference() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )

    assert bundle.provenance is not None
    assert bundle.provenance.bundle_id == "l6_native"
    assert bundle.provenance.organism == Organism.RAT.value
    assert bundle.provenance.source_name
    assert bundle.provenance.identifier_namespace
    assert bundle.provenance.sequence_window is not None
    assert bundle.provenance.manifest is not None
    assert bundle.provenance.manifest.get("files") is not None


def test_runtime_reference_loader_does_not_treat_test_fixtures_as_bundled_runtime_data(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()
    patched_defaults = dict(reference_resources._BUNDLED_DEFAULTS)
    patched_defaults[Organism.RAT] = "fragile_support_reference"
    monkeypatch.setattr(reference_resources, "_BUNDLED_DEFAULTS", patched_defaults)

    with pytest.raises(ReferenceResolutionError, match="directory is missing"):
        load_bundled_reference_manifest(Organism.RAT)


def _write_temporary_manifest_bundle(tmp_path: Path, *, file_hash: str) -> None:
    data = "kinase,site_id\nAKT1,MAPK1;S123;\n"
    (tmp_path / "reference.csv").write_text(data, encoding="utf-8")
    payload = {
        "reference_id": "unit_reference",
        "display_name": "Unit reference",
        "organism": "Rattus norvegicus",
        "organism_common_name": "rat",
        "taxonomy_id": 10116,
        "protein_namespace": "display_site_id",
        "reference_version": "v1",
        "source_name": "unit source",
        "source_version": "source-v1",
        "source_url": None,
        "source_publication": None,
        "source_license": "CC0 synthetic",
        "source_license_url": None,
        "redistribution_allowed": True,
        "redistribution_notes": "synthetic test fixture",
        "derived_from": ["unit test"],
        "generated_by": "unit test",
        "generated_at_utc": "2026-06-29T00:00:00Z",
        "manifest_schema_version": "1.0",
        "files": [
            {
                "relative_path": "reference.csv",
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
    (tmp_path / "manifest.json").write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
