from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.validation import (
    ReferenceIdentifierNormalisationValidationError,
    ReferenceValidationError,
)
from phospy.provenance.models import (
    EnvironmentProvenance,
    ReferenceProvenance,
    RunProvenance,
)
from phospy.provenance.serialization import from_payload, to_payload
from phospy.science.references.identifiers import (
    normalise_reference_kinase_id,
    normalise_reference_protein_accession,
    normalise_reference_site_id,
)
from phospy.science.references.models import Organism, ReferenceBundle
from phospy.tables.references import (
    KinaseSubstrateReference,
    ProteinAccessionReference,
    SiteSequenceReference,
)


def test_kinase_identifier_normaliser_normalises_mixed_case_values() -> None:
    records = [
        normalise_reference_kinase_id(
            "akt1",
            table_name="references.kinase_substrate_map",
            column_name="kinase",
            row_position=0,
        ),
        normalise_reference_kinase_id(
            "Akt1",
            table_name="references.kinase_substrate_map",
            column_name="kinase",
            row_position=1,
        ),
        normalise_reference_kinase_id(
            "AKT1",
            table_name="references.kinase_substrate_map",
            column_name="kinase",
            row_position=2,
        ),
    ]

    assert [record.normalised_value for record in records] == ["AKT1", "AKT1", "AKT1"]


def test_site_identifier_normaliser_normalises_supported_inputs() -> None:
    records = [
        normalise_reference_site_id(
            "mapk1;s123",
            table_name="references.kinase_substrate_map",
            column_name="substrate_site",
            row_position=0,
        ),
        normalise_reference_site_id(
            "Mapk1 ; s123",
            table_name="references.kinase_substrate_map",
            column_name="substrate_site",
            row_position=1,
        ),
        normalise_reference_site_id(
            "MAPK1;S123;",
            table_name="references.kinase_substrate_map",
            column_name="substrate_site",
            row_position=2,
        ),
    ]

    assert [record.normalised_value for record in records] == [
        "MAPK1;S123;",
        "MAPK1;S123;",
        "MAPK1;S123;",
    ]


def test_reference_bundle_rejects_duplicate_pairs_after_kinase_normalisation() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["akt1", "AKT1"],
                    "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK1;S123;"], name="site_id"),
            ),
        )


def test_reference_bundle_kinase_substrate_conflict_preserves_structured_report() -> (
    None
):
    frame = pd.DataFrame(
        {
            "kinase": ["akt1", "AKT1"],
            "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
            "source_db": ["db_a", "db_b"],
        }
    )
    with pytest.raises(
        ReferenceIdentifierNormalisationValidationError,
        match="contains conflicting payload rows for normalised",
    ) as direct_exc:
        KinaseSubstrateReference(frame=frame)
    with pytest.raises(
        ReferenceValidationError,
        match="contains conflicting payload rows for normalised",
    ) as bundle_exc:
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=frame,
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK1;S123;"], name="site_id"),
            ),
        )

    direct_report = direct_exc.value.identifier_normalisation_report
    report = getattr(bundle_exc.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report == direct_report
    assert report.conflict_count == 4
    assert report.duplicate_identifier_count == 0
    assert any(
        record.table_name == "references.kinase_substrate_map"
        and record.status == "conflict_after_normalisation"
        for record in report.records
    )


def test_reference_bundle_site_sequence_conflict_preserves_structured_report() -> None:
    frame = pd.DataFrame(
        {"site_sequence": ["A" * 31, "B" * 31]},
        index=pd.Index(["mapk1;s123", "MAPK1;S123;"], name="site_id"),
    )
    with pytest.raises(
        ReferenceIdentifierNormalisationValidationError,
        match="conflicting site_sequence values after normalisation",
    ) as direct_exc:
        SiteSequenceReference(frame=frame)
    with pytest.raises(
        ReferenceValidationError,
        match="conflicting site_sequence values after normalisation",
    ) as bundle_exc:
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {"kinase": ["AKT1"], "substrate_site": ["MAPK1;S123;"]}
            ),
            site_sequences=frame,
        )

    direct_report = direct_exc.value.identifier_normalisation_report
    report = getattr(bundle_exc.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report == direct_report
    assert report.conflict_count == 2
    assert report.duplicate_identifier_count == 0
    assert any(
        record.table_name == "references.site_sequences"
        and record.status == "conflict_after_normalisation"
        and record.column_name == "index"
        for record in report.records
    )


def test_reference_bundle_duplicate_after_normalisation_preserves_structured_report() -> (
    None
):
    frame = pd.DataFrame(
        {
            "kinase": ["akt1", "AKT1"],
            "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
        }
    )
    with pytest.raises(
        ReferenceIdentifierNormalisationValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ) as direct_exc:
        KinaseSubstrateReference(frame=frame)
    with pytest.raises(
        ReferenceValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ) as bundle_exc:
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=frame,
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK1;S123;"], name="site_id"),
            ),
        )

    direct_report = direct_exc.value.identifier_normalisation_report
    report = getattr(bundle_exc.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report == direct_report
    assert report.duplicate_identifier_count == 4
    assert report.conflict_count == 0
    assert any(
        record.table_name == "references.kinase_substrate_map"
        and record.status == "duplicate_after_normalisation"
        for record in report.records
    )


def test_invalid_site_identifier_failure_exposes_structured_report() -> None:
    with pytest.raises(
        ReferenceIdentifierNormalisationValidationError,
        match="site identifiers must use 'GENE;SITE;' format",
    ) as exc_info:
        SiteSequenceReference(
            frame=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK1-S123"], name="site_id"),
            )
        )

    report = exc_info.value.identifier_normalisation_report
    assert report.original_row_count == 1
    assert report.normalised_row_count == 0
    assert report.invalid_identifier_count == 1
    assert report.changed_identifier_count == 0
    assert report.duplicate_identifier_count == 0
    assert report.conflict_count == 0
    assert len(report.records) == 1
    assert report.records[0].original_value == "MAPK1-S123"
    assert report.records[0].normalised_value is None
    assert report.records[0].status == "invalid"


def test_duplicate_kinase_identifier_failure_exposes_structured_report() -> None:
    with pytest.raises(
        ReferenceIdentifierNormalisationValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ) as exc_info:
        KinaseSubstrateReference(
            frame=pd.DataFrame(
                {
                    "kinase": ["akt1", "AKT1"],
                    "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
                }
            )
        )

    report = exc_info.value.identifier_normalisation_report
    assert report.invalid_identifier_count == 0
    assert report.changed_identifier_count == 1
    assert report.duplicate_identifier_count == 4
    assert report.conflict_count == 0
    assert any(
        record.original_value == "akt1" and record.normalised_value == "AKT1"
        for record in report.records
    )
    duplicate_records = [
        record
        for record in report.records
        if record.status == "duplicate_after_normalisation"
    ]
    assert duplicate_records


def test_conflict_after_normalisation_failure_exposes_structured_report() -> None:
    with pytest.raises(
        ReferenceIdentifierNormalisationValidationError,
        match="contains conflicting payload rows for normalised",
    ) as exc_info:
        KinaseSubstrateReference(
            frame=pd.DataFrame(
                {
                    "kinase": ["akt1", "AKT1"],
                    "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
                    "source_db": ["db_a", "db_b"],
                }
            )
        )

    report = exc_info.value.identifier_normalisation_report
    assert report.invalid_identifier_count == 0
    assert report.changed_identifier_count == 1
    assert report.duplicate_identifier_count == 0
    assert report.conflict_count == 4
    assert any(
        record.original_value == "akt1" and record.normalised_value == "AKT1"
        for record in report.records
    )
    conflict_records = [
        record
        for record in report.records
        if record.status == "conflict_after_normalisation"
    ]
    assert conflict_records


def test_site_sequence_duplicate_after_normalisation_reports_duplicate_status() -> None:
    with pytest.raises(ReferenceValidationError) as exc_info:
        SiteSequenceReference(
            frame=pd.DataFrame(
                {"site_sequence": ["A" * 31, "A" * 31]},
                index=pd.Index(["mapk1;s123", "MAPK1;S123;"], name="site_id"),
            )
        )

    report = getattr(exc_info.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report.duplicate_identifier_count == 2
    assert report.conflict_count == 0
    classified = [
        record
        for record in report.records
        if record.status.endswith("_after_normalisation")
    ]
    assert {record.status for record in classified} == {"duplicate_after_normalisation"}


def test_site_sequence_conflict_after_normalisation_reports_conflict_status() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="conflicting site_sequence values after normalisation",
    ) as exc_info:
        SiteSequenceReference(
            frame=pd.DataFrame(
                {"site_sequence": ["A" * 31, "B" * 31]},
                index=pd.Index(["mapk1;s123", "MAPK1;S123;"], name="site_id"),
            )
        )

    report = getattr(exc_info.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report.duplicate_identifier_count == 0
    assert report.conflict_count == 2
    conflict_records = [
        record
        for record in report.records
        if record.status == "conflict_after_normalisation"
    ]
    assert len(conflict_records) == 2


def test_kinase_substrate_duplicate_after_normalisation_reports_duplicate_status() -> (
    None
):
    with pytest.raises(
        ReferenceValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ) as exc_info:
        KinaseSubstrateReference(
            frame=pd.DataFrame(
                {
                    "kinase": ["akt1", "AKT1"],
                    "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
                }
            )
        )

    report = getattr(exc_info.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report.duplicate_identifier_count == 4
    assert report.conflict_count == 0
    classified = [
        record
        for record in report.records
        if record.status.endswith("_after_normalisation")
    ]
    assert {record.status for record in classified} == {"duplicate_after_normalisation"}


def test_kinase_substrate_conflict_after_normalisation_with_divergent_payload() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="contains conflicting payload rows for normalised",
    ) as exc_info:
        KinaseSubstrateReference(
            frame=pd.DataFrame(
                {
                    "kinase": ["akt1", "AKT1"],
                    "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
                    "source_db": ["db_a", "db_b"],
                }
            )
        )

    report = getattr(exc_info.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report.duplicate_identifier_count == 0
    assert report.conflict_count == 4
    conflict_records = [
        record
        for record in report.records
        if record.status == "conflict_after_normalisation"
    ]
    assert len(conflict_records) == 4


def test_kinase_substrate_without_payload_columns_cannot_report_conflicts() -> None:
    with pytest.raises(ReferenceValidationError) as exc_info:
        KinaseSubstrateReference(
            frame=pd.DataFrame(
                {
                    "kinase": ["akt1", "AKT1"],
                    "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
                }
            )
        )

    report = getattr(exc_info.value, "identifier_normalisation_report", None)
    assert report is not None
    assert report.conflict_count == 0
    assert report.duplicate_identifier_count > 0


def test_reference_identifier_normalisation_does_not_fuzzy_match_or_add_synonyms() -> (
    None
):
    site_record = normalise_reference_site_id(
        "MAPK1-S123",
        table_name="references.kinase_substrate_map",
        column_name="substrate_site",
        row_position=0,
    )
    kinase_record = normalise_reference_kinase_id(
        "akt-1",
        table_name="references.kinase_substrate_map",
        column_name="kinase",
        row_position=0,
    )

    assert site_record.status == "invalid"
    assert kinase_record.normalised_value == "AKT-1"


def test_protein_accession_normaliser_trims_whitespace_only() -> None:
    trimmed = normalise_reference_protein_accession(
        " P28482 ",
        table_name="tests.reference.protein_accession_fixture",
        column_name="protein_accession",
        row_position=0,
    )
    mixed_case = normalise_reference_protein_accession(
        " p28482 ",
        table_name="tests.reference.protein_accession_fixture",
        column_name="protein_accession",
        row_position=1,
    )

    assert trimmed.original_value == " P28482 "
    assert trimmed.normalised_value == "P28482"
    assert trimmed.status == "normalised"
    assert mixed_case.normalised_value == "p28482"
    assert mixed_case.normalised_value != "P28482"
    assert mixed_case.status == "normalised"


def test_protein_accession_reference_reports_duplicates_after_trim() -> None:
    with pytest.raises(
        ReferenceIdentifierNormalisationValidationError,
        match="duplicate protein_accession after normalisation",
    ) as exc_info:
        ProteinAccessionReference(
            frame=pd.DataFrame(
                {
                    "protein_accession": [" P28482 ", "P28482"],
                }
            )
        )

    report = exc_info.value.identifier_normalisation_report
    assert report.changed_identifier_count == 1
    assert report.duplicate_identifier_count == 2
    assert report.conflict_count == 0
    duplicate_records = [
        record
        for record in report.records
        if record.status == "duplicate_after_normalisation"
    ]
    assert len(duplicate_records) == 2
    assert {record.normalised_value for record in duplicate_records} == {"P28482"}
    provenance = ReferenceProvenance(
        source_type="explicit",
        organism="rat",
        bundle_id=None,
        table_fingerprints=(),
        identifier_normalisation=report,
    )
    payload = to_payload(
        RunProvenance(
            environment=EnvironmentProvenance(
                package_name="phospy",
                package_version="test",
                python_version="3.13",
                dependency_versions={},
                platform={},
            ),
            input_tables=(),
            preprocessing_stages=(),
            reference=provenance,
            workflow_name="unit_test",
            workflow_parameters={},
            random_state=None,
            random_seed_policy=None,
            output_tables=(),
        )
    )
    restored = from_payload(payload)
    assert restored.reference is not None
    assert restored.reference.identifier_normalisation is not None
    assert restored.reference.identifier_normalisation.duplicate_identifier_count == 2


def test_protein_accession_reference_reports_conflicts_after_trim() -> None:
    with pytest.raises(
        ReferenceIdentifierNormalisationValidationError,
        match="conflicting payload for protein_accession after normalisation",
    ) as exc_info:
        ProteinAccessionReference(
            frame=pd.DataFrame(
                {
                    "protein_accession": [" P28482 ", "P28482"],
                    "source_db": ["db_a", "db_b"],
                }
            )
        )

    report = exc_info.value.identifier_normalisation_report
    assert report.changed_identifier_count == 1
    assert report.duplicate_identifier_count == 0
    assert report.conflict_count == 2
    conflict_records = [
        record
        for record in report.records
        if record.status == "conflict_after_normalisation"
    ]
    assert len(conflict_records) == 2


def test_protein_accession_production_ingestion_normalises_values() -> None:
    reference = ProteinAccessionReference(
        frame=pd.DataFrame(
            {
                "protein_accession": [" P28482 "],
                "source_db": ["phosphositeplus"],
            }
        )
    )

    assert reference.frame.loc[0, "protein_accession"] == "P28482"
    report = reference.identifier_normalisation
    assert report is not None
    accession_records = [
        record for record in report.records if record.column_name == "protein_accession"
    ]
    assert len(accession_records) == 1
    assert accession_records[0].status == "normalised"


def test_protein_accession_normalisation_records_are_serialisable() -> None:
    reference = ProteinAccessionReference(
        frame=pd.DataFrame(
            {
                "protein_accession": [" P28482 ", "p28482"],
                "source_db": ["db_a", "db_b"],
            }
        )
    )
    report = reference.identifier_normalisation
    assert report is not None
    provenance = ReferenceProvenance(
        source_type="explicit",
        organism="rat",
        bundle_id=None,
        table_fingerprints=(),
        identifier_normalisation=report,
    )
    payload = to_payload(
        RunProvenance(
            environment=EnvironmentProvenance(
                package_name="phospy",
                package_version="test",
                python_version="3.13",
                dependency_versions={},
                platform={},
            ),
            input_tables=(),
            preprocessing_stages=(),
            reference=provenance,
            workflow_name="unit_test",
            workflow_parameters={},
            random_state=None,
            random_seed_policy=None,
            output_tables=(),
        )
    )
    restored = from_payload(payload)

    assert restored.reference is not None
    assert restored.reference.identifier_normalisation is not None
    restored_records = restored.reference.identifier_normalisation.records
    assert len(restored_records) == 2
    assert restored_records[0].identifier_kind == "protein_accession"
    assert restored_records[0].original_value == " P28482 "
    assert restored_records[0].normalised_value == "P28482"
