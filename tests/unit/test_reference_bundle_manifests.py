from __future__ import annotations

import pytest

from phospy.errors.references import ReferenceResolutionError, UnsupportedOrganismError
from phospy.references import resources as reference_resources
from phospy.references.models import Organism, ReferencePreset
from phospy.references.resolution import ReferenceResolver
from phospy.references.resources import load_bundled_reference_manifest


def test_valid_bundled_manifest_loads_for_supported_runtime_lane() -> None:
    manifest = load_bundled_reference_manifest(Organism.RAT)

    assert manifest.bundle_id == "l6_native"
    assert manifest.organism == "Rattus norvegicus"
    assert manifest.organism_common_name == Organism.RAT.value
    assert manifest.identifier_namespace
    assert manifest.source_name
    assert manifest.sequence_window is not None
    assert manifest.sequence_window.upstream_residues == 15
    assert manifest.sequence_window.downstream_residues == 15
    assert "site_sequence_derivation" in manifest.supports
    assert manifest.limitations


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
    with pytest.raises(ReferenceResolutionError, match="manifest bundle_id"):
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
        UnsupportedOrganismError, match="supported bundled organisms: rat"
    ):
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.HUMAN,
        )


def test_auto_reference_resolution_succeeds_for_supported_runtime_organism() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.AUTO,
        dataset_organism=Organism.RAT,
    )

    assert bundle.manifest is not None
    assert bundle.manifest.organism_common_name == Organism.RAT.value


def test_auto_reference_resolution_fails_for_unsupported_runtime_organism() -> None:
    with pytest.raises(
        UnsupportedOrganismError, match="supported bundled organisms: rat"
    ):
        ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=Organism.HUMAN,
        )


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
