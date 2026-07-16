from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import pandas as pd
import pytest

from phospy.errors.references import ReferenceResolutionError
from phospy.errors.validation import ReferenceValidationError
from phospy.io.bundles.reference_sources import ReferenceSourceTableReader
from phospy.provenance.models import (
    EnvironmentProvenance,
    ReferenceProvenance,
    RunProvenance,
)
from phospy.provenance.serialization import from_payload, to_payload
from phospy.science.references import resolution as reference_resolution
from phospy.science.references.builder import ReferenceBundleBuilder
from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.models import (
    Organism,
    RedistributionStatus,
    ReferenceBundleBuildRequest,
    ReferenceContext,
    ReferenceFileManifest,
    ReferenceManifest,
    ReferencePreset,
    SequenceWindowDefinition,
    reference_context_from_provenance,
)
from phospy.science.references.resolution import (
    BundledReferenceProvider,
    ReferenceResolver,
)
from phospy.science.references.validation import load_reference_manifest

LOCAL_REFERENCE_VERSION = "local-snapshot-v1"
UPSTREAM_SOURCE_VERSION = "upstream-source-v1"


def test_bundled_provider_propagates_manifest_source_version() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )

    assert bundle.manifest is not None
    assert bundle.provenance is not None
    assert bundle.manifest.source_version == "PhosR 1.20.0"
    assert bundle.manifest.reference_version == "bundled-snapshot-2026-04-16"
    assert bundle.provenance.source_version == "PhosR 1.20.0"
    assert bundle.provenance.reference_context is not None
    assert bundle.provenance.reference_context.source_version == "PhosR 1.20.0"
    assert bundle.provenance.manifest is not None
    assert bundle.provenance.manifest["source_version"] == "PhosR 1.20.0"
    assert (
        bundle.provenance.manifest["reference_version"] == "bundled-snapshot-2026-04-16"
    )
    assert bundle.validation_report.provenance_fields["source_version"] == (
        "PhosR 1.20.0"
    )
    assert bundle.validation_report.bundle_version == "bundled-snapshot-2026-04-16"


def test_bundled_provider_rejects_missing_manifest_source_version_before_tables(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = _manifest(source_version=None)
    table_loads: list[str] = []

    monkeypatch.setattr(
        reference_resolution,
        "bundled_reference_name_for_organism",
        lambda organism: manifest.bundle_id,
    )
    monkeypatch.setattr(
        reference_resolution,
        "load_bundled_reference_manifest",
        lambda organism: manifest,
    )

    def _table_loader(organism: Organism) -> pd.DataFrame:
        table_loads.append(organism.value)
        raise AssertionError("tables must not load without manifest source_version")

    monkeypatch.setattr(
        reference_resolution,
        "load_bundled_kinase_substrate_map",
        _table_loader,
    )
    monkeypatch.setattr(
        reference_resolution,
        "load_bundled_site_sequences",
        _table_loader,
    )

    with pytest.raises(ReferenceResolutionError, match="missing source_version"):
        BundledReferenceProvider().run(Organism.RAT)

    assert table_loads == []


def test_stock_builder_preserves_distinct_explicit_reference_and_source_versions(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)

    bundle = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        _build_request(
            kinase_path,
            sequence_path,
            reference_version=f" {LOCAL_REFERENCE_VERSION} ",
        )
    )

    assert bundle.manifest is not None
    assert bundle.provenance is not None
    assert bundle.manifest.reference_version == LOCAL_REFERENCE_VERSION
    assert bundle.manifest.source_version == UPSTREAM_SOURCE_VERSION
    assert bundle.provenance.source_version == UPSTREAM_SOURCE_VERSION
    assert bundle.provenance.reference_context is not None
    assert bundle.provenance.reference_context.source_version == (
        UPSTREAM_SOURCE_VERSION
    )
    assert bundle.provenance.manifest is not None
    assert bundle.provenance.manifest["reference_version"] == LOCAL_REFERENCE_VERSION
    assert bundle.provenance.manifest["source_version"] == UPSTREAM_SOURCE_VERSION


def test_stock_builder_generates_deterministic_content_reference_version_when_omitted(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)

    first = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        ReferenceBundleBuildRequest(
            organism=Organism.RAT,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="unit upstream source",
            source_version=UPSTREAM_SOURCE_VERSION,
            retrieved_at="2026-06-11",
            license="synthetic test license",
            redistribution_status="redistributable synthetic fixture",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
            bundle_id="unit_local_snapshot",
        )
    )
    second = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        ReferenceBundleBuildRequest(
            organism=Organism.RAT,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="unit upstream source",
            source_version=UPSTREAM_SOURCE_VERSION,
            retrieved_at="2026-06-11",
            license="synthetic test license",
            redistribution_status="redistributable synthetic fixture",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
            bundle_id="unit_local_snapshot",
        )
    )

    expected = _expected_generated_reference_version(kinase_path, sequence_path)

    assert first.manifest is not None
    assert second.manifest is not None
    assert first.manifest.reference_version == expected
    assert second.manifest.reference_version == expected
    assert first.manifest.source_version == UPSTREAM_SOURCE_VERSION


def test_generated_reference_version_changes_when_kinase_source_bytes_change(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)
    first = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        _build_request(kinase_path, sequence_path, reference_version=None)
    )

    pd.DataFrame(
        {
            "kinase": [" akt1 ", "Map3k1"],
            "site_id": [" mapk1 ; s123 ", "Mapk1;Y185;"],
            "organism": ["rat", "rat"],
        }
    ).to_csv(kinase_path, index=False)

    second = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        _build_request(kinase_path, sequence_path, reference_version=None)
    )

    assert first.manifest is not None
    assert second.manifest is not None
    assert first.manifest.reference_version != second.manifest.reference_version
    assert second.manifest.reference_version == _expected_generated_reference_version(
        kinase_path,
        sequence_path,
    )


def test_generated_reference_version_changes_when_sequence_source_bytes_change(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)
    first = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        _build_request(kinase_path, sequence_path, reference_version=None)
    )

    pd.DataFrame(
        {
            "site_id": [" mapk1 ; s123 ", "Mapk1;Y185;"],
            "site_sequence": [_window("S"), _window("Y")],
            "display_id": ["mapk1;s123", "MAPK1;Y185;"],
            "organism": ["Rattus norvegicus", "Rattus norvegicus"],
        }
    ).to_csv(sequence_path, index=False)

    second = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        _build_request(kinase_path, sequence_path, reference_version=None)
    )

    assert first.manifest is not None
    assert second.manifest is not None
    assert first.manifest.reference_version != second.manifest.reference_version
    assert second.manifest.reference_version == _expected_generated_reference_version(
        kinase_path,
        sequence_path,
    )


def test_builder_provenance_and_context_use_upstream_source_version(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)

    bundle = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        _build_request(
            kinase_path,
            sequence_path,
            reference_version=LOCAL_REFERENCE_VERSION,
        )
    )

    assert bundle.provenance is not None
    assert bundle.provenance.source_version == UPSTREAM_SOURCE_VERSION
    assert bundle.provenance.reference_context is not None
    assert bundle.provenance.reference_context.source_version == (
        UPSTREAM_SOURCE_VERSION
    )
    context = reference_context_from_provenance(bundle.provenance)
    assert context is not None
    assert context.source_version == UPSTREAM_SOURCE_VERSION
    assert bundle.validation_report.provenance_fields["source_version"] == (
        UPSTREAM_SOURCE_VERSION
    )


def test_builder_serialization_preserves_local_and_upstream_versions(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)
    bundle = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        _build_request(
            kinase_path,
            sequence_path,
            reference_version=LOCAL_REFERENCE_VERSION,
        )
    )

    assert bundle.provenance is not None
    restored = from_payload(to_payload(_run_provenance(bundle.provenance)))

    assert restored.reference is not None
    assert restored.reference.source_version == UPSTREAM_SOURCE_VERSION
    assert restored.reference.reference_context is not None
    assert restored.reference.reference_context.source_version == (
        UPSTREAM_SOURCE_VERSION
    )
    assert restored.reference.manifest is not None
    assert restored.reference.manifest["reference_version"] == LOCAL_REFERENCE_VERSION
    assert restored.reference.manifest["source_version"] == UPSTREAM_SOURCE_VERSION


def test_existing_request_without_reference_version_remains_valid(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)

    bundle = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        ReferenceBundleBuildRequest(
            organism=Organism.RAT,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="unit upstream source",
            source_version=UPSTREAM_SOURCE_VERSION,
            retrieved_at="2026-06-11",
            license="synthetic test license",
            redistribution_status="redistributable synthetic fixture",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
            bundle_id="unit_local_snapshot",
        )
    )

    assert bundle.manifest is not None
    assert bundle.manifest.reference_version.startswith("local-snapshot-sha256-")
    assert bundle.manifest.source_version == UPSTREAM_SOURCE_VERSION


def test_existing_optional_positional_arguments_are_not_reinterpreted(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)
    sequence_window = SequenceWindowDefinition(
        upstream_residues=15,
        downstream_residues=15,
        central_residue_required=True,
    )

    request = ReferenceBundleBuildRequest(
        Organism.RAT,
        kinase_path,
        sequence_path,
        "unit upstream source",
        UPSTREAM_SOURCE_VERSION,
        "2026-06-11",
        "synthetic test license",
        "redistributable synthetic fixture",
        "display_id (GENE_SYMBOL;RESIDUE;)",
        sequence_window,
        "positional_bundle",
        "rat",
        ("kinase_workflow",),
        ("legacy positional limitation",),
    )

    assert request.sequence_window is sequence_window
    assert request.bundle_id == "positional_bundle"
    assert request.organism_common_name == "rat"
    assert request.supports == ("kinase_workflow",)
    assert request.limitations == ("legacy positional limitation",)
    assert request.reference_version is None

    bundle = ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
        request
    )

    assert bundle.manifest is not None
    assert bundle.manifest.bundle_id == "positional_bundle"


def test_blank_reference_version_fails_build_request_validation(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(tmp_path)

    with pytest.raises(ReferenceResolutionError, match="reference_version"):
        ReferenceBundleBuilder(source_reader=ReferenceSourceTableReader()).run(
            _build_request(
                kinase_path,
                sequence_path,
                reference_version=" ",
            )
        )


def test_manifest_context_and_provenance_source_version_agreement_succeeds() -> None:
    provenance = ReferenceProvenance(
        source_type="bundled",
        organism="rat",
        bundle_id="unit_local_snapshot",
        source_name="unit upstream source",
        source_version=f" {UPSTREAM_SOURCE_VERSION} ",
        identifier_namespace="display_id",
        manifest={
            "reference_version": LOCAL_REFERENCE_VERSION,
            "source_version": UPSTREAM_SOURCE_VERSION,
        },
        table_fingerprints=(),
        reference_context=_context(source_version=UPSTREAM_SOURCE_VERSION),
    )

    assert provenance.source_version == UPSTREAM_SOURCE_VERSION


def test_provenance_context_source_version_contradiction_fails() -> None:
    with pytest.raises(ReferenceValidationError) as exc_info:
        ReferenceProvenance(
            source_type="bundled",
            organism="rat",
            bundle_id="unit_local_snapshot",
            source_name="unit upstream source",
            source_version=LOCAL_REFERENCE_VERSION,
            identifier_namespace="display_id",
            table_fingerprints=(),
            reference_context=_context(source_version=UPSTREAM_SOURCE_VERSION),
        )

    message = str(exc_info.value)
    assert "Reference provenance source-version mismatch" in message
    assert "provenance.source_version" in message
    assert "reference_context.source_version" in message


def test_provenance_manifest_source_version_contradiction_fails() -> None:
    with pytest.raises(ReferenceValidationError) as exc_info:
        ReferenceProvenance(
            source_type="bundled",
            organism="rat",
            bundle_id="unit_local_snapshot",
            source_name="unit upstream source",
            source_version=LOCAL_REFERENCE_VERSION,
            identifier_namespace="display_id",
            manifest={
                "reference_version": LOCAL_REFERENCE_VERSION,
                "source_version": UPSTREAM_SOURCE_VERSION,
            },
            table_fingerprints=(),
        )

    message = str(exc_info.value)
    assert "Reference provenance source-version mismatch" in message
    assert "provenance.source_version" in message
    assert "manifest.source_version" in message


def test_missing_source_version_does_not_fallback_to_reference_version() -> None:
    provenance = _legacy_unknown_provenance()

    assert provenance.source_version is None
    assert reference_context_from_provenance(provenance) is None


def test_serialized_missing_source_version_does_not_fallback_to_reference_version() -> (
    None
):
    provenance = _legacy_unknown_provenance()

    restored = from_payload(to_payload(_run_provenance(provenance)))

    assert restored.reference is not None
    assert restored.reference.source_version is None
    assert restored.reference.reference_context is None
    assert reference_context_from_provenance(restored.reference) is None
    assert restored.reference.manifest is not None
    assert restored.reference.manifest["reference_version"] == LOCAL_REFERENCE_VERSION


@pytest.mark.parametrize(
    "placeholder",
    [
        "unknown",
        "unspecified",
        "n/a",
        "na",
        "none",
        "null",
        "tbd",
        "not specified",
    ],
)
def test_placeholder_approved_source_version_fails_release_validation(
    tmp_path: Path,
    placeholder: str,
) -> None:
    manifest_path = _write_manifest_bundle(tmp_path, source_version=placeholder)

    with pytest.raises(ReferenceManifestError, match="placeholder"):
        load_reference_manifest(manifest_path, bundled=True)


def _context(*, source_version: str) -> ReferenceContext:
    return ReferenceContext(
        organism="rat",
        protein_namespace="display_id",
        source_name="unit upstream source",
        source_version=source_version,
        proteome_version=None,
        reference_table_sha256="a" * 64,
    )


def _manifest(
    *, source_version: str | None = UPSTREAM_SOURCE_VERSION
) -> ReferenceManifest:
    return ReferenceManifest(
        reference_id="unit_local_snapshot",
        display_name="Unit local snapshot",
        organism="Rattus norvegicus",
        taxonomy_id=10116,
        organism_common_name="rat",
        protein_namespace="display_id",
        reference_version=LOCAL_REFERENCE_VERSION,
        source_name="unit upstream source",
        source_url="https://example.test/reference",
        source_version=source_version,
        retrieved_at="2026-06-29",
        table_sha256="a" * 64,
        license_name="unit license",
        license_url="https://example.test/license",
        redistribution_status=RedistributionStatus.APPROVED,
        redistribution_notes="synthetic fixture",
        derived_from=("unit test",),
        generated_by="unit test",
        generated_at_utc="2026-06-29T00:00:00Z",
        manifest_schema_version="1.1",
        files=(
            ReferenceFileManifest(
                relative_path="reference.csv",
                role="kinase_substrate",
                format="csv",
                sha256="a" * 64,
            ),
        ),
    )


def _legacy_unknown_provenance() -> ReferenceProvenance:
    return ReferenceProvenance(
        source_type="legacy",
        organism="rat",
        bundle_id="unit_local_snapshot",
        source_name="unit upstream source",
        source_version=None,
        identifier_namespace="display_id",
        manifest={
            "organism": "rat",
            "protein_namespace": "display_id",
            "source_name": "unit upstream source",
            "source_version": None,
            "reference_version": LOCAL_REFERENCE_VERSION,
            "table_sha256": "a" * 64,
        },
        table_fingerprints=(),
    )


def _run_provenance(reference: ReferenceProvenance) -> RunProvenance:
    return RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="test",
            python_version="3.12",
            dependency_versions={},
            platform={},
        ),
        input_tables=(),
        preprocessing_stages=(),
        reference=reference,
        workflow_name="unit_test",
        workflow_parameters={},
        random_state=None,
        random_seed_policy=None,
        output_tables=(),
    )


def _build_request(
    kinase_path: Path,
    sequence_path: Path,
    *,
    reference_version: str | None,
) -> ReferenceBundleBuildRequest:
    return ReferenceBundleBuildRequest(
        organism=Organism.RAT,
        kinase_substrate_path=kinase_path,
        site_sequence_path=sequence_path,
        source_name="unit upstream source",
        source_version=UPSTREAM_SOURCE_VERSION,
        retrieved_at="2026-06-11",
        license="synthetic test license",
        redistribution_status="redistributable synthetic fixture",
        identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
        bundle_id="unit_local_snapshot",
        reference_version=reference_version,
    )


def _expected_generated_reference_version(
    kinase_path: Path,
    sequence_path: Path,
) -> str:
    canonical = (
        f"kinase_substrate:{sha256(kinase_path.read_bytes()).hexdigest()}\n"
        f"site_sequences:{sha256(sequence_path.read_bytes()).hexdigest()}\n"
    )
    return f"local-snapshot-sha256-{sha256(canonical.encode('ascii')).hexdigest()}"


def _write_reference_sources(tmp_path: Path) -> tuple[Path, Path]:
    kinase_path = tmp_path / "kinase.csv"
    sequence_path = tmp_path / "sequences.csv"
    pd.DataFrame(
        {
            "kinase": [" akt1 ", "Map2k1"],
            "site_id": [" mapk1 ; s123 ", "Mapk1;Y185;"],
            "organism": ["rat", "rat"],
        }
    ).to_csv(kinase_path, index=False)
    pd.DataFrame(
        {
            "site_id": [" mapk1 ; s123 ", "Mapk1;Y185;"],
            "site_sequence": [_window("S").lower(), _window("Y")],
            "display_id": ["mapk1;s123", "MAPK1;Y185;"],
            "organism": ["rat", "rat"],
        }
    ).to_csv(sequence_path, index=False)
    return kinase_path, sequence_path


def _write_manifest_bundle(tmp_path: Path, *, source_version: str) -> Path:
    data_path = tmp_path / "reference.csv"
    data = "kinase,site_id\nAKT1,MAPK1;S123;\n"
    data_path.write_text(data, encoding="utf-8")
    file_hash = sha256(data_path.read_bytes()).hexdigest()
    manifest_payload: dict[str, object] = {
        "reference_id": "unit_local_snapshot",
        "display_name": "Unit local snapshot",
        "organism": "Rattus norvegicus",
        "taxonomy_id": 10116,
        "protein_namespace": "display_id",
        "reference_version": LOCAL_REFERENCE_VERSION,
        "source_name": "unit upstream source",
        "source_version": source_version,
        "source_url": "https://example.test/reference",
        "retrieved_at": "2026-06-29",
        "table_sha256": file_hash,
        "source_publication": None,
        "license_name": "CC0 synthetic",
        "license_url": "https://example.test/license",
        "redistribution_status": "approved",
        "redistribution_allowed": True,
        "redistribution_notes": "synthetic test fixture",
        "derived_from": ["unit test"],
        "generated_by": "unit test",
        "generated_at_utc": "2026-06-29T00:00:00Z",
        "manifest_schema_version": "1.1",
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
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest_payload, indent=2), encoding="utf-8")
    return manifest_path


def _window(center: str) -> str:
    return f"{'A' * 15}{center}{'A' * 15}"
