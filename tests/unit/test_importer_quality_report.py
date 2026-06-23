from __future__ import annotations

from phospy.api import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_NOT_REPORTED,
    IMPORTER_QUALITY_STATUS_REPORTED,
    ImporterDetectedIntensityColumn,
    ImporterDuplicateKeySummary,
    ImporterFlaggedRowSummary,
    ImporterLocalisationConfidenceSummary,
    ImporterMissingIntensitySummary,
    ImporterQualityCount,
    ImporterQualityReport,
)


def test_full_importer_quality_report_construction_and_payload() -> None:
    report = ImporterQualityReport(
        source_name="mapped_fixture",
        row_count_status=IMPORTER_QUALITY_STATUS_REPORTED,
        rows_read=4,
        rows_retained=3,
        rows_dropped=1,
        intensity_column_status=IMPORTER_QUALITY_STATUS_REPORTED,
        detected_intensity_columns=(
            ImporterDetectedIntensityColumn(
                source_column="Intensity A",
                sample_id="sample_a",
            ),
            ImporterDetectedIntensityColumn(
                source_column="Intensity B",
                sample_id="sample_b",
            ),
        ),
        missing_intensity=ImporterMissingIntensitySummary(
            status=IMPORTER_QUALITY_STATUS_REPORTED,
            total_missing_values=2,
            rows_with_any_missing_intensity=1,
            missing_values_by_sample_id={"sample_a": 0, "sample_b": 2},
            missing_values_by_source_column={"Intensity A": 0, "Intensity B": 2},
        ),
        localisation_confidence=ImporterLocalisationConfidenceSummary(
            status=IMPORTER_QUALITY_STATUS_REPORTED,
            source_column="Localization prob",
            output_column="localisation_confidence",
            scale="probability",
            row_count=3,
            missing_count=1,
            invalid_count=1,
            invalid_examples=("row2:'bad':not_numeric",),
        ),
        flagged_rows=ImporterFlaggedRowSummary(
            contaminant=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_REPORTED,
                count=1,
                source_column="Potential contaminant",
                policy="remove",
            ),
            reverse=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_REPORTED,
                count=0,
                source_column="Reverse",
                policy="remove",
            ),
            decoy=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="format does not report decoys",
            ),
        ),
        duplicate_keys=ImporterDuplicateKeySummary(
            site_key=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_REPORTED,
                count=1,
                source_column="site_key",
            ),
            display_key=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_REPORTED,
                count=2,
                source_column="display_id",
            ),
            duplicate_site_candidate_rows=1,
        ),
        format_specific={"maxquant": {"removed_rows": 1}},
        warnings=("duplicate site candidates were retained",),
    )

    payload = report.to_payload()

    assert payload["source_name"] == "mapped_fixture"
    assert payload["rows_read"] == 4
    assert payload["rows_retained"] == 3
    assert payload["rows_dropped"] == 1
    assert payload["detected_intensity_columns"] == [
        {"source_column": "Intensity A", "sample_id": "sample_a"},
        {"source_column": "Intensity B", "sample_id": "sample_b"},
    ]
    assert payload["missing_intensity"]["total_missing_values"] == 2
    assert payload["localisation_confidence"]["invalid_count"] == 1
    assert payload["flagged_rows"]["contaminant"]["count"] == 1
    assert payload["flagged_rows"]["decoy"]["status"] == (
        IMPORTER_QUALITY_STATUS_NOT_APPLICABLE
    )
    assert payload["duplicate_keys"]["display_key"]["count"] == 2
    assert payload["format_specific"] == {"maxquant": {"removed_rows": 1}}
    assert payload["warnings"] == ["duplicate site candidates were retained"]


def test_minimal_importer_quality_report_uses_explicit_not_reported_fields() -> None:
    report = ImporterQualityReport()
    payload = report.to_payload()

    assert report.row_count_status == IMPORTER_QUALITY_STATUS_NOT_REPORTED
    assert report.rows_read is None
    assert payload["intensity_column_status"] == IMPORTER_QUALITY_STATUS_NOT_REPORTED
    assert payload["detected_intensity_columns"] == []
    assert payload["missing_intensity"]["status"] == (
        IMPORTER_QUALITY_STATUS_NOT_REPORTED
    )
    assert payload["localisation_confidence"]["status"] == (
        IMPORTER_QUALITY_STATUS_NOT_REPORTED
    )


def test_importer_quality_report_preserves_warnings() -> None:
    report = ImporterQualityReport(warnings=("check duplicate upstream rows",))

    assert report.warnings == ("check duplicate upstream rows",)
    assert report.to_payload()["warnings"] == ["check duplicate upstream rows"]


def test_importer_quality_report_represents_not_applicable_fields() -> None:
    report = ImporterQualityReport(
        localisation_confidence=ImporterLocalisationConfidenceSummary(
            status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
            reason="localisation column was not mapped",
        ),
        flagged_rows=ImporterFlaggedRowSummary(
            contaminant=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="contaminant flags are not present",
            ),
            reverse=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="reverse flags are not present",
            ),
            decoy=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="decoy flags are not present",
            ),
        ),
        duplicate_keys=ImporterDuplicateKeySummary(
            site_key=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="site_key column was not mapped",
            ),
            display_key=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="display_id column was not mapped",
            ),
        ),
    )

    payload = report.to_payload()

    assert payload["localisation_confidence"]["status"] == (
        IMPORTER_QUALITY_STATUS_NOT_APPLICABLE
    )
    assert payload["flagged_rows"]["contaminant"]["reason"] == (
        "contaminant flags are not present"
    )
    assert payload["duplicate_keys"]["site_key"]["status"] == (
        IMPORTER_QUALITY_STATUS_NOT_APPLICABLE
    )
