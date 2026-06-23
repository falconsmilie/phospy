from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    Organism,
    ReferenceBundle,
    ReferenceBundleValidationReport,
)
from phospy.errors.validation import ReferenceValidationError
from phospy.science.references.models import (
    ReferenceManifest,
    SequenceWindowDefinition,
)
from phospy.validation.references.bundle import ReferenceBundleValidator


def test_valid_reference_bundle_exposes_structured_validation_report() -> None:
    bundle = _bundle(manifest=_manifest())

    report = bundle.validation_report

    assert isinstance(report, ReferenceBundleValidationReport)
    assert report.bundle_name == "unit_reference"
    assert report.bundle_version == "v1"
    assert report.organism == "Rattus norvegicus"
    assert report.organism_common_name == "rat"
    assert report.identifier_namespace == "display_id (GENE_SYMBOL;RESIDUE;)"
    assert report.kinase_substrate_record_count == 2
    assert report.duplicate_record_count == 0
    assert report.duplicate_records == ()
    assert report.compatibility_warnings == ()

    table_status = {item.table_name: item for item in report.required_tables}
    assert set(table_status) == {"kinase_substrate_map", "site_sequences"}
    assert table_status["kinase_substrate_map"].present is True
    assert table_status["kinase_substrate_map"].required_columns == (
        "kinase",
        "substrate_site",
    )
    assert table_status["kinase_substrate_map"].missing_required_columns == ()
    assert table_status["site_sequences"].missing_values == ()

    file_status = {item.role: item for item in report.required_source_files}
    assert file_status["kinase_substrate"].present is True
    assert file_status["site_sequences"].present is True
    assert file_status["kinase_substrate"].path == "kinase.csv"

    payload = report.to_payload()
    assert payload["bundle_name"] == "unit_reference"
    assert payload["kinase_substrate_record_count"] == 2


def test_report_warns_when_required_source_file_metadata_is_missing() -> None:
    bundle = _bundle(
        manifest=_manifest(
            source_files={
                "kinase_substrate": {
                    "path": "kinase.csv",
                    "role": "kinase-substrate source",
                }
            }
        )
    )

    report = bundle.validation_report

    file_status = {item.role: item.present for item in report.required_source_files}
    assert file_status == {"kinase_substrate": True, "site_sequences": False}
    assert any(
        "source-file metadata is incomplete" in warning for warning in report.warnings
    )


def test_validator_rejects_missing_required_table() -> None:
    with pytest.raises(ReferenceValidationError, match="must be a pandas DataFrame"):
        ReferenceBundleValidator().run(
            organism=Organism.RAT,
            kinase_substrate_map=None,
            site_sequences=_site_sequences(),
        )


def test_validator_rejects_missing_required_column() -> None:
    invalid_sequences = pd.DataFrame(
        {"sequence": ["AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA"]},
        index=pd.Index(["MAPK1;S123;"], name="site_id"),
    )

    with pytest.raises(ReferenceValidationError, match="missing required columns"):
        ReferenceBundleValidator().run(
            organism=Organism.RAT,
            kinase_substrate_map=_kinase_substrate_map(),
            site_sequences=invalid_sequences,
        )


def test_duplicate_reference_records_are_rejected() -> None:
    duplicate_map = pd.DataFrame(
        {
            "kinase": ["AKT1", "AKT1"],
            "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
        }
    )

    with pytest.raises(
        ReferenceValidationError,
        match="duplicate \\(kinase, substrate_site\\) pairs",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=duplicate_map,
            site_sequences=_site_sequences(),
        )


def test_report_warns_when_organism_and_namespace_metadata_are_limited() -> None:
    bundle = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=_kinase_substrate_map(),
        site_sequences=_site_sequences(),
    )

    report = bundle.validation_report

    assert report.organism == "rat"
    assert report.organism_common_name is None
    assert report.identifier_namespace is None
    assert any("organism metadata is limited" in item for item in report.warnings)
    assert any("identifier namespace metadata" in item for item in report.warnings)


def test_report_includes_available_provenance_fields() -> None:
    bundle = _bundle(manifest=_manifest())

    fields = bundle.validation_report.provenance_fields

    assert fields["source_name"] == "unit reference source"
    assert fields["source_version"] == "v1"
    assert fields["license"] == "unit test license"
    assert fields["redistribution_status"] == "redistributable synthetic fixture"
    assert fields["source_files_available"] is True


def test_invalid_bundle_still_rejected_when_sequence_entries_are_missing() -> None:
    with pytest.raises(ReferenceValidationError, match="missing sequence entries"):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {"kinase": ["AKT1"], "substrate_site": ["MAPK1;S123;"]}
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA"]},
                index=pd.Index(["MAPK3;T202;"], name="site_id"),
            ),
        )


def _bundle(*, manifest: ReferenceManifest) -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=_kinase_substrate_map(),
        site_sequences=_site_sequences(),
        manifest=manifest,
    )


def _kinase_substrate_map() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kinase": ["AKT1", "MAP2K1"],
            "substrate_site": ["MAPK1;S123;", "MAPK3;T202;"],
        }
    )


def _site_sequences() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_sequence": [
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ]
        },
        index=pd.Index(["MAPK1;S123;", "MAPK3;T202;"], name="site_id"),
    )


def _manifest(
    *,
    source_files: dict[str, object] | None = None,
) -> ReferenceManifest:
    return ReferenceManifest(
        bundle_id="unit_reference",
        organism="Rattus norvegicus",
        organism_common_name="rat",
        identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
        source_name="unit reference source",
        source_version="v1",
        retrieved_at=pd.Timestamp("2026-06-23").date(),
        license="unit test license",
        redistribution_status="redistributable synthetic fixture",
        sequence_window=SequenceWindowDefinition(
            upstream_residues=15,
            downstream_residues=15,
            central_residue_required=True,
        ),
        supports=("kinase_workflow",),
        limitations=("unit test fixture",),
        source_files=(
            {
                "kinase_substrate": {
                    "path": "kinase.csv",
                    "role": "kinase-substrate source",
                },
                "site_sequences": {
                    "path": "sequences.csv",
                    "role": "site-sequence source",
                },
            }
            if source_files is None
            else source_files
        ),
    )
