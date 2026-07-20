from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.errors.workflows import WorkflowStageError
from phospy.provenance import RowAttritionRecord, RowAttritionReport
from phospy.workflows._row_attrition import (
    make_row_attrition_record,
    reconcile_row_attrition_report,
)


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


def test_row_attrition_payload_mutation_does_not_change_report() -> None:
    report = RowAttritionReport.from_records(
        (
            RowAttritionRecord(
                stage="site_matrix",
                input_rows=3,
                output_rows=1,
                removed_rows=2,
                reason="duplicate_site_resolution",
                examples=("SITE_A", "SITE_B"),
            ),
        )
    )

    payload = report.to_payload()
    records = payload["records"]
    assert isinstance(records, list)
    first_record = records[0]
    assert isinstance(first_record, dict)
    examples = first_record["examples"]
    assert isinstance(examples, list)
    examples.append("PAYLOAD_ONLY")

    assert report.to_payload()["records"][0]["examples"] == ["SITE_A", "SITE_B"]


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


def test_make_row_attrition_record_samples_removed_ids_deterministically() -> None:
    record = make_row_attrition_record(
        workflow="test_workflow",
        stage="filter_stage",
        reason="test_filter",
        input_site_ids=pd.Index(["S1", "S2", "S3", "S4", "S5", "S6", "S7"]),
        output_site_ids=pd.Index(["S7"]),
    )

    assert record is not None
    assert record.removed_rows == 6
    assert record.examples == ("S1", "S2", "S3", "S4", "S5")


def test_make_row_attrition_record_omits_zero_removal_records() -> None:
    record = make_row_attrition_record(
        workflow="test_workflow",
        stage="filter_stage",
        reason="test_filter",
        input_site_ids=pd.Index(["S1", "S2"]),
        output_site_ids=pd.Index(["S1", "S2"]),
    )

    assert record is None


def test_make_row_attrition_record_rejects_output_site_absent_from_input() -> None:
    with pytest.raises(
        WorkflowStageError,
        match="stage=filter_stage; input_count=2; output_count=2",
    ):
        make_row_attrition_record(
            workflow="test_workflow",
            stage="filter_stage",
            reason="test_filter",
            input_site_ids=pd.Index(["S1", "S2"]),
            output_site_ids=pd.Index(["S1", "S3"]),
        )


def test_row_attrition_reconciliation_rejects_discontinuity() -> None:
    records = (
        RowAttritionRecord(
            stage="first",
            input_rows=4,
            output_rows=3,
            removed_rows=1,
            reason="first_filter",
        ),
        RowAttritionRecord(
            stage="second",
            input_rows=2,
            output_rows=1,
            removed_rows=1,
            reason="second_filter",
        ),
    )

    with pytest.raises(
        WorkflowStageError,
        match="stage=second; input_count=2; output_count=1; expected_preceding_count=3",
    ):
        reconcile_row_attrition_report(
            workflow="test_workflow",
            records=records,
            initial_site_ids=pd.Index(["S1", "S2", "S3", "S4"]),
            final_site_ids=pd.Index(["S1"]),
        )


def test_row_attrition_reconciliation_rejects_wrong_final_count() -> None:
    records = (
        RowAttritionRecord(
            stage="first",
            input_rows=4,
            output_rows=3,
            removed_rows=1,
            reason="first_filter",
        ),
    )

    with pytest.raises(
        WorkflowStageError,
        match="stage=first; input_count=4; output_count=3",
    ):
        reconcile_row_attrition_report(
            workflow="test_workflow",
            records=records,
            initial_site_ids=pd.Index(["S1", "S2", "S3", "S4"]),
            final_site_ids=pd.Index(["S1", "S2"]),
        )
