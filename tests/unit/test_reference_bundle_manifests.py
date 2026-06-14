from __future__ import annotations

import json
from importlib import resources

import pytest

from phospy.errors.references import (
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    UnsupportedOrganismError,
)
from phospy.science.references import resources as reference_resources
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
_REDISTRIBUTION_APPROVAL_FIELDS = (
    "source_url",
    "license_url",
    "retrieval_method",
    "redistribution_basis",
)
_APPROVAL_REQUIRED_ORGANISMS = (Organism.HUMAN, Organism.MOUSE)


def _unsupported_approval_required_organism() -> Organism:
    supported = set(supported_bundled_organisms())
    for organism in _APPROVAL_REQUIRED_ORGANISMS:
        if organism not in supported:
            return organism
    pytest.skip("all human/mouse bundled lanes are committed in this test run")


def _valid_manifest_payload(
    organism: Organism = Organism.RAT,
    *,
    approved_redistribution: bool = False,
) -> dict[str, object]:
    organism_names = {
        Organism.HUMAN: "Homo sapiens",
        Organism.MOUSE: "Mus musculus",
        Organism.RAT: "Rattus norvegicus",
    }
    payload: dict[str, object] = {
        "bundle_id": "l6_native",
        "organism": organism_names[organism],
        "organism_common_name": organism.value,
        "identifier_namespace": "site_id (GENE_SYMBOL;RESIDUE;)",
        "source_name": "unit test source",
        "source_version": "unit-test-snapshot",
        "retrieved_at": "2026-04-16",
        "license": "unit test license",
        "redistribution_status": (
            "redistribution approved for unit tests"
            if approved_redistribution
            else "unit test redistribution status"
        ),
        "sequence_window": {
            "upstream_residues": 15,
            "downstream_residues": 15,
            "central_residue_required": True,
        },
        "supports": ["site_sequence_derivation"],
        "limitations": ["unit test limitation"],
    }
    if approved_redistribution:
        payload.update(
            {
                "source_url": "https://example.invalid/unit-test-source",
                "license_url": "https://example.invalid/unit-test-license",
                "retrieval_method": "unit test manifest fixture",
                "redistribution_basis": (
                    "unit test source declares redistribution allowed"
                ),
            }
        )
    return payload


def test_valid_bundled_manifest_loads_for_supported_runtime_lane() -> None:
    manifest = load_bundled_reference_manifest(Organism.RAT)

    assert manifest.bundle_id == "l6_native"
    assert manifest.organism == "Rattus norvegicus"
    assert manifest.organism_common_name == Organism.RAT.value
    assert manifest.identifier_namespace
    assert manifest.source_name
    assert manifest.source_version == "bundled-snapshot-2026-04-16"
    assert manifest.retrieved_at.isoformat() == "2026-04-16"
    assert manifest.source_url == "https://github.com/PYangLab/PhosR"
    assert manifest.license_url
    assert manifest.retrieval_method
    assert manifest.redistribution_basis
    assert manifest.license
    assert "gpl-3" in manifest.license.lower()
    assert "phosphositeplus" in manifest.license.lower()
    assert "pride" in manifest.license.lower()
    assert manifest.redistribution_status
    assert "not independently verified" in manifest.redistribution_status.lower()
    assert manifest.source_files is not None
    assert "substrate_map" in manifest.source_files
    assert manifest.provenance_notes is not None
    assert any(
        "not independent" in note.lower()
        or "not captured" in note.lower()
        or "not independently" in note.lower()
        for note in manifest.provenance_notes
    )
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
    for organism in supported_bundled_organisms():
        bundle_id = bundled_reference_name_for_organism(organism)
        lane_root = (
            resources.files("phospy")
            .joinpath("data")
            .joinpath("reference_bundles")
            .joinpath(organism.value)
            .joinpath(bundle_id)
        )
        for filename in _REQUIRED_RUNTIME_FILES:
            assert lane_root.joinpath(filename).is_file(), (
                f"missing bundled reference runtime file: "
                f"{organism.value}/{bundle_id}/{filename}"
            )

        manifest = load_bundled_reference_manifest(organism)
        assert manifest.bundle_id == bundle_id
        assert manifest.organism
        assert manifest.organism_common_name == organism.value
        assert manifest.identifier_namespace
        assert manifest.source_name
        assert manifest.source_version
        assert manifest.retrieved_at.isoformat()
        assert manifest.license
        assert manifest.redistribution_status
        assert manifest.supports
        assert manifest.limitations
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
    assert bundle.manifest.redistribution_status
    assert not bundle.kinase_substrate_map.empty
    assert not bundle.site_sequences.empty


@pytest.mark.parametrize("organism", _APPROVAL_REQUIRED_ORGANISMS)
def test_human_mouse_bundle_directories_are_not_stray_or_unapproved(
    organism: Organism,
) -> None:
    organism_root = (
        resources.files("phospy")
        .joinpath("data")
        .joinpath("reference_bundles")
        .joinpath(organism.value)
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

    manifest_resource = organism_root.joinpath(default_bundle).joinpath("manifest.json")
    raw_manifest = json.loads(manifest_resource.read_text(encoding="utf-8"))
    for field_name in _REDISTRIBUTION_APPROVAL_FIELDS:
        assert raw_manifest.get(field_name), (
            f"{organism.value}/{default_bundle}/manifest.json must include {field_name}"
        )
    load_bundled_reference_manifest(organism)


def test_bundled_manifest_missing_resource_fails_clearly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()

    def _missing_resource(**_: object) -> object:
        raise ReferenceResolutionError("bundled reference resource is missing")

    monkeypatch.setattr(reference_resources, "_read_json_resource", _missing_resource)
    with pytest.raises(ReferenceResolutionError, match="resource is missing"):
        load_bundled_reference_manifest(Organism.RAT)


def test_bundled_manifest_malformed_payload_fails_clearly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()

    def _malformed_payload(**_: object) -> object:
        return {
            "organism": Organism.RAT.value,
            "source_name": "missing required fields",
        }

    monkeypatch.setattr(reference_resources, "_read_json_resource", _malformed_payload)
    with pytest.raises(ReferenceResolutionError, match="missing required field"):
        load_bundled_reference_manifest(Organism.RAT)


@pytest.mark.parametrize(
    "missing_field",
    [
        "organism",
        "bundle_id",
        "identifier_namespace",
        "source_name",
        "source_version",
        "retrieved_at",
        "license",
        "redistribution_status",
        "sequence_window",
        "supports",
        "limitations",
    ],
)
def test_bundled_manifest_missing_required_metadata_fails_clearly(
    monkeypatch: pytest.MonkeyPatch,
    missing_field: str,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()
    payload = _valid_manifest_payload()
    payload.pop(missing_field)

    def _missing_metadata_payload(**_: object) -> object:
        return payload

    monkeypatch.setattr(
        reference_resources,
        "_read_json_resource",
        _missing_metadata_payload,
    )
    with pytest.raises(
        ReferenceResolutionError,
        match=rf"missing required field\(s\).*{missing_field}",
    ):
        load_bundled_reference_manifest(Organism.RAT)


def test_bundled_manifest_preserves_optional_provenance_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()
    payload = _valid_manifest_payload()
    payload.update(
        {
            "source_files": {
                "substrate_map": {
                    "path": "reference/substrate_map.csv",
                    "sha256": "a" * 64,
                }
            },
            "provenance_notes": ["unit test provenance note"],
        }
    )

    def _manifest_with_optional_metadata(**_: object) -> object:
        return payload

    monkeypatch.setattr(
        reference_resources,
        "_read_json_resource",
        _manifest_with_optional_metadata,
    )

    manifest = load_bundled_reference_manifest(Organism.RAT)

    assert manifest.source_files == payload["source_files"]
    assert manifest.provenance_notes == ("unit test provenance note",)
    manifest_payload = manifest.to_payload()
    assert manifest_payload["source_files"] == payload["source_files"]
    assert manifest_payload["provenance_notes"] == ("unit test provenance note",)


@pytest.mark.parametrize(
    ("field_name", "invalid_value", "message"),
    [
        ("source_files", [], "source_files must be an object"),
        ("source_files", {}, "source_files must not be empty"),
        ("provenance_notes", ["valid", ""], r"provenance_notes\[1\]"),
    ],
)
def test_bundled_manifest_rejects_invalid_optional_provenance_metadata(
    monkeypatch: pytest.MonkeyPatch,
    field_name: str,
    invalid_value: object,
    message: str,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()
    payload = _valid_manifest_payload()
    payload[field_name] = invalid_value

    def _manifest_with_invalid_optional_metadata(**_: object) -> object:
        return payload

    monkeypatch.setattr(
        reference_resources,
        "_read_json_resource",
        _manifest_with_invalid_optional_metadata,
    )

    with pytest.raises(ReferenceResolutionError, match=message):
        load_bundled_reference_manifest(Organism.RAT)


def test_human_mouse_bundled_manifest_requires_redistribution_approval_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()
    patched_defaults = dict(reference_resources._BUNDLED_DEFAULTS)
    patched_defaults[Organism.HUMAN] = "test_human"
    monkeypatch.setattr(reference_resources, "_BUNDLED_DEFAULTS", patched_defaults)

    def _human_manifest_without_approval_metadata(**_: object) -> object:
        return _valid_manifest_payload(Organism.HUMAN)

    monkeypatch.setattr(
        reference_resources,
        "_read_json_resource",
        _human_manifest_without_approval_metadata,
    )

    with pytest.raises(
        ReferenceResolutionError,
        match="requires redistribution approval metadata",
    ):
        load_bundled_reference_manifest(Organism.HUMAN)


def test_human_mouse_bundled_manifest_rejects_unapproved_redistribution_status(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()
    patched_defaults = dict(reference_resources._BUNDLED_DEFAULTS)
    patched_defaults[Organism.MOUSE] = "test_mouse"
    monkeypatch.setattr(reference_resources, "_BUNDLED_DEFAULTS", patched_defaults)
    payload = _valid_manifest_payload(Organism.MOUSE, approved_redistribution=True)
    payload["redistribution_status"] = "redistribution unclear pending review"

    def _mouse_manifest_with_unapproved_status(**_: object) -> object:
        return payload

    monkeypatch.setattr(
        reference_resources,
        "_read_json_resource",
        _mouse_manifest_with_unapproved_status,
    )

    with pytest.raises(
        ReferenceResolutionError,
        match="redistribution_status must explicitly say",
    ):
        load_bundled_reference_manifest(Organism.MOUSE)


def test_human_mouse_bundled_manifest_preserves_approval_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()
    patched_defaults = dict(reference_resources._BUNDLED_DEFAULTS)
    patched_defaults[Organism.HUMAN] = "test_human"
    monkeypatch.setattr(reference_resources, "_BUNDLED_DEFAULTS", patched_defaults)
    payload = _valid_manifest_payload(Organism.HUMAN, approved_redistribution=True)

    def _approved_human_manifest(**_: object) -> object:
        return payload

    monkeypatch.setattr(
        reference_resources,
        "_read_json_resource",
        _approved_human_manifest,
    )

    try:
        manifest = load_bundled_reference_manifest(Organism.HUMAN)

        assert manifest.source_url == payload["source_url"]
        assert manifest.license_url == payload["license_url"]
        assert manifest.retrieval_method == payload["retrieval_method"]
        assert manifest.redistribution_basis == payload["redistribution_basis"]
        manifest_payload = manifest.to_payload()
        for field_name in _REDISTRIBUTION_APPROVAL_FIELDS:
            assert manifest_payload[field_name] == payload[field_name]
    finally:
        reference_resources.clear_bundled_reference_manifest_cache()


def test_reference_resolution_succeeds_for_supported_runtime_organism() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )

    assert bundle.organism is Organism.RAT
    assert bundle.manifest is not None
    assert bundle.manifest.bundle_id == "l6_native"


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
    assert bundle.provenance.manifest.get("supports") is not None


def test_runtime_reference_loader_does_not_treat_test_fixtures_as_bundled_runtime_data(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference_resources.clear_bundled_reference_manifest_cache()
    patched_defaults = dict(reference_resources._BUNDLED_DEFAULTS)
    patched_defaults[Organism.RAT] = "fragile_support_reference"
    monkeypatch.setattr(reference_resources, "_BUNDLED_DEFAULTS", patched_defaults)

    with pytest.raises(ReferenceResolutionError, match="resource is missing"):
        load_bundled_reference_manifest(Organism.RAT)
