from __future__ import annotations

from collections.abc import Iterator, Mapping

import numpy as np
import pytest

from phospy.contracts.results import (
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
from phospy.errors import PhosPyInputError


class _DuplicateKeyMapping(Mapping[object, object]):
    def __iter__(self) -> Iterator[object]:
        return iter(("duplicate", "duplicate"))

    def __len__(self) -> int:
        return 2

    def __getitem__(self, key: object) -> object:
        if key == "duplicate":
            return "value"
        raise KeyError(key)


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


def test_importer_quality_format_specific_is_recursively_immutable() -> None:
    source_metadata = {
        "maxquant": {
            "removed_rows": [1],
            "columns": {"score": "Localization prob"},
        }
    }
    report = ImporterQualityReport(format_specific=source_metadata)

    source_metadata["maxquant"]["removed_rows"].append(2)
    source_metadata["maxquant"]["columns"]["score"] = "mutated"

    maxquant = report.format_specific["maxquant"]
    assert isinstance(maxquant, Mapping)
    assert maxquant["removed_rows"] == (1,)
    assert maxquant["columns"]["score"] == "Localization prob"

    with pytest.raises(TypeError):
        maxquant["columns"]["score"] = "changed"  # type: ignore[index]

    payload = report.to_payload()
    payload["format_specific"]["maxquant"]["removed_rows"].append(3)  # type: ignore[union-attr]
    payload["format_specific"]["maxquant"]["columns"]["score"] = "payload"  # type: ignore[index]

    assert report.to_payload()["format_specific"] == {
        "maxquant": {
            "removed_rows": [1],
            "columns": {"score": "Localization prob"},
        }
    }


@pytest.mark.parametrize(
    ("format_specific", "expected"),
    (
        ({1: "numeric-key"}, "keys must be strings"),
        (_DuplicateKeyMapping(), "duplicate JSON object key"),
        ({"bad": float("nan")}, "finite JSON number"),
        ({"bad": float("inf")}, "finite JSON number"),
        ({"bad": {"unsupported"}}, "got set"),
        ({"bad": np.array([1])}, "got ndarray"),
        ({"bad": bytearray(b"x")}, "got bytearray"),
        ({"bad": object()}, "got object"),
    ),
)
def test_importer_quality_format_specific_rejects_invalid_json_values(
    format_specific: Mapping[object, object],
    expected: str,
) -> None:
    with pytest.raises(PhosPyInputError) as exc_info:
        ImporterQualityReport(format_specific=format_specific)

    message = str(exc_info.value)
    assert "importer_quality.format_specific" in message
    assert expected in message


def test_importer_missing_intensity_count_mappings_are_immutable() -> None:
    by_sample = {"sample_a": 1}
    by_source = {"Intensity A": 1}
    summary = ImporterMissingIntensitySummary(
        status=IMPORTER_QUALITY_STATUS_REPORTED,
        total_missing_values=1,
        rows_with_any_missing_intensity=1,
        missing_values_by_sample_id=by_sample,
        missing_values_by_source_column=by_source,
    )

    by_sample["sample_a"] = 99
    by_source["Intensity A"] = 99

    assert summary.missing_values_by_sample_id == {"sample_a": 1}
    assert summary.missing_values_by_source_column == {"Intensity A": 1}

    with pytest.raises(TypeError):
        summary.missing_values_by_sample_id["sample_a"] = 2  # type: ignore[index]

    payload = summary.to_payload()
    payload["missing_values_by_sample_id"]["sample_a"] = 3  # type: ignore[index]

    assert summary.to_payload()["missing_values_by_sample_id"] == {"sample_a": 1}


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
