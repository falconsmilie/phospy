from __future__ import annotations

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.provenance import RowAttritionRecord, RowAttritionReport


def test_row_attrition_record_accepts_valid_counts() -> None:
    record = RowAttritionRecord(
        stage="missing_data",
        input_rows=5,
        output_rows=3,
        removed_rows=2,
        reason="rows with excessive missing values",
        examples=("SITE_A", "SITE_B"),
    )

    assert record.input_rows == 5
    assert record.output_rows == 3
    assert record.removed_rows == 2
    assert record.examples == ("SITE_A", "SITE_B")
    assert record.to_payload() == {
        "stage": "missing_data",
        "input_rows": 5,
        "output_rows": 3,
        "removed_rows": 2,
        "reason": "rows with excessive missing values",
        "examples": ["SITE_A", "SITE_B"],
    }


@pytest.mark.parametrize(
    "kwargs",
    [
        {"input_rows": -1, "output_rows": 0, "removed_rows": 0},
        {"input_rows": 1, "output_rows": -1, "removed_rows": 2},
        {"input_rows": 1, "output_rows": 1, "removed_rows": -1},
    ],
)
def test_row_attrition_record_rejects_negative_counts(
    kwargs: dict[str, int],
) -> None:
    with pytest.raises(PhosPyInputError, match="must be non-negative"):
        RowAttritionRecord(
            stage="missing_data",
            reason="rows with excessive missing values",
            **kwargs,
        )


@pytest.mark.parametrize(
    "kwargs, message",
    [
        (
            {"input_rows": 3, "output_rows": 4, "removed_rows": 0},
            "output_rows must be less than or equal to input_rows",
        ),
        (
            {"input_rows": 5, "output_rows": 3, "removed_rows": 1},
            "removed_rows must equal input_rows - output_rows",
        ),
    ],
)
def test_row_attrition_record_rejects_inconsistent_counts(
    kwargs: dict[str, int],
    message: str,
) -> None:
    with pytest.raises(PhosPyInputError, match=message):
        RowAttritionRecord(
            stage="missing_data",
            reason="rows with excessive missing values",
            **kwargs,
        )


def test_row_attrition_record_rejects_empty_stage() -> None:
    with pytest.raises(PhosPyInputError, match="stage must be non-empty"):
        RowAttritionRecord(
            stage=" ",
            input_rows=1,
            output_rows=1,
            removed_rows=0,
            reason="no rows removed",
        )


def test_row_attrition_record_rejects_empty_reason() -> None:
    with pytest.raises(PhosPyInputError, match="reason must be non-empty"):
        RowAttritionRecord(
            stage="missing_data",
            input_rows=1,
            output_rows=1,
            removed_rows=0,
            reason="",
        )


def test_row_attrition_report_composes_typed_records() -> None:
    records = (
        RowAttritionRecord(
            stage="site_matrix",
            input_rows=5,
            output_rows=4,
            removed_rows=1,
            reason="invalid site identifiers",
        ),
        RowAttritionRecord(
            stage="missing_data",
            input_rows=4,
            output_rows=3,
            removed_rows=1,
            reason="rows with excessive missing values",
        ),
    )

    report = RowAttritionReport.from_records(records)

    assert report.input_rows == 5
    assert report.final_rows == 3
    assert report.to_payload()["records"] == [record.to_payload() for record in records]


def test_row_attrition_report_rejects_plain_dict_records() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="RowAttritionRecord",
    ):
        RowAttritionReport(
            records=(
                {
                    "stage": "missing_data",
                    "input_rows": 2,
                    "output_rows": 1,
                    "removed_rows": 1,
                    "reason": "rows with excessive missing values",
                },
            ),
            input_rows=2,
            final_rows=1,
        )
