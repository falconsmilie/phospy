from __future__ import annotations

import pytest

from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.science.references import resources as reference_resources
from phospy.science.references.models import Organism, ReferencePreset
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.references.resources import (
    available_bundled_reference_lanes,
    load_bundled_reference_manifest,
)


def _valid_manifest_payload() -> dict[str, object]:
    return {
        "bundle_id": "l6_native",
        "organism": "Rattus norvegicus",
        "organism_common_name": Organism.RAT.value,
        "identifier_namespace": "site_id (GENE_SYMBOL;RESIDUE;)",
        "source_name": "unit test source",
        "source_version": "unit-test-snapshot",
        "retrieved_at": "2026-04-16",
        "license": "unit test license",
        "redistribution_status": "unit test redistribution status",
        "sequence_window": {
            "upstream_residues": 15,
            "downstream_residues": 15,
            "central_residue_required": True,
        },
        "supports": ["site_sequence_derivation"],
        "limitations": ["unit test limitation"],
    }


def test_valid_bundled_manifest_loads_for_supported_runtime_lane() -> None:
    manifest = load_bundled_reference_manifest(Organism.RAT)

    assert manifest.bundle_id == "l6_native"
    assert manifest.organism == "Rattus norvegicus"
    assert manifest.organism_common_name == Organism.RAT.value
    assert manifest.identifier_namespace
    assert manifest.source_name
    assert manifest.source_version == "bundled-snapshot-2026-04-16"
    assert manifest.retrieved_at.isoformat() == "2026-04-16"
    assert manifest.license
    assert manifest.redistribution_status
    assert manifest.sequence_window.upstream_residues == 15
    assert manifest.sequence_window.downstream_residues == 15
    assert "site_sequence_derivation" in manifest.supports
    assert manifest.limitations


def test_available_bundled_reference_lanes_reports_manifest_metadata() -> None:
    lanes = available_bundled_reference_lanes()

    assert len(lanes) == 1
    lane = lanes[0]
    assert lane.organism is Organism.RAT
    assert lane.bundle_id == "l6_native"
    assert lane.source_name
    assert lane.source_version == "bundled-snapshot-2026-04-16"
    assert lane.retrieved_at.isoformat() == "2026-04-16"
    assert lane.redistribution_status
    assert "site_sequence_derivation" in lane.supports
    assert lane.limitations
    assert lane.to_payload()["organism"] == Organism.RAT.value


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


def test_reference_resolution_succeeds_for_supported_runtime_organism() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )

    assert bundle.organism is Organism.RAT
    assert bundle.manifest is not None
    assert bundle.manifest.bundle_id == "l6_native"


def test_reference_resolution_fails_for_unsupported_runtime_organism() -> None:
    with pytest.raises(
        UnsupportedOrganismError,
        match="reference-bundle docs: https://phospy.com/docs/api/guide/#references",
    ) as exc_info:
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.HUMAN,
        )
    message = str(exc_info.value)
    assert "supported bundled organisms: rat" in message
    assert "ReferenceBundle(organism=..." in message


def test_auto_reference_resolution_succeeds_for_supported_runtime_organism() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.AUTO,
        dataset_organism=Organism.RAT,
    )

    assert bundle.manifest is not None
    assert bundle.manifest.organism_common_name == Organism.RAT.value


def test_auto_reference_resolution_fails_for_unsupported_runtime_organism() -> None:
    with pytest.raises(
        UnsupportedOrganismError,
        match="reference-bundle docs: https://phospy.com/docs/api/guide/#references",
    ) as exc_info:
        ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=Organism.HUMAN,
        )
    assert "supported bundled organisms: rat" in str(exc_info.value)


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
