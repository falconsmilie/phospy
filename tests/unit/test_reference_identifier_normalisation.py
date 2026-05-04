from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.validation import ReferenceValidationError
from phospy.references.identifiers import (
    normalise_reference_kinase_id,
    normalise_reference_site_id,
)
from phospy.references.models import Organism, ReferenceBundle
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference


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
        match="conflicting site_sequence values after canonicalization",
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
