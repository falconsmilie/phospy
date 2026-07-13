"""Private helpers for causal workflow row-attrition records."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.provenance.models import RowAttritionRecord, RowAttritionReport

_REMOVED_ID_SAMPLE_LIMIT = 5


def make_row_attrition_record(
    *,
    workflow: str,
    stage: str,
    reason: str,
    input_site_ids: pd.Index | Iterable[object],
    output_site_ids: pd.Index | Iterable[object],
) -> RowAttritionRecord | None:
    """Create a typed record from one immediate filter-only stage transition."""

    input_ids = _site_id_tuple(input_site_ids)
    output_ids = _site_id_tuple(output_site_ids)
    input_count = int(len(input_ids))
    output_count = int(len(output_ids))
    if output_count > input_count:
        _raise_row_attrition_error(
            workflow=workflow,
            stage=stage,
            input_count=input_count,
            output_count=output_count,
            detail="output count exceeds input count for filter-only stage",
        )
    input_set = set(input_ids)
    unexpected_output_ids = tuple(
        site_id for site_id in output_ids if site_id not in input_set
    )
    if unexpected_output_ids:
        _raise_row_attrition_error(
            workflow=workflow,
            stage=stage,
            input_count=input_count,
            output_count=output_count,
            detail=(
                "stage output contains site IDs absent from its input; "
                f"unexpected_output_site_examples={list(_examples(unexpected_output_ids))}"
            ),
        )
    output_set = set(output_ids)
    removed_ids = tuple(site_id for site_id in input_ids if site_id not in output_set)
    removed_count = int(len(removed_ids))
    expected_removed_count = int(input_count - output_count)
    if removed_count != expected_removed_count:
        _raise_row_attrition_error(
            workflow=workflow,
            stage=stage,
            input_count=input_count,
            output_count=output_count,
            detail=(
                "removed count does not equal input minus output; "
                f"expected_removed_count={expected_removed_count}; "
                f"observed_removed_count={removed_count}"
            ),
        )
    if removed_count == 0:
        return None
    return RowAttritionRecord(
        stage=str(stage),
        input_rows=input_count,
        output_rows=output_count,
        removed_rows=removed_count,
        reason=str(reason),
        examples=_examples(removed_ids),
    )


def reconcile_row_attrition_report(
    *,
    workflow: str,
    records: Iterable[RowAttritionRecord],
    initial_site_ids: pd.Index | Iterable[object],
    final_site_ids: pd.Index | Iterable[object],
) -> RowAttritionReport | None:
    """Validate stage continuity and reconcile typed records to final rows."""

    record_tuple = tuple(records)
    initial_count = int(len(_site_id_tuple(initial_site_ids)))
    final_count = int(len(_site_id_tuple(final_site_ids)))
    if not record_tuple:
        if initial_count != final_count:
            _raise_row_attrition_error(
                workflow=workflow,
                stage="row_attrition_reconciliation",
                input_count=initial_count,
                output_count=final_count,
                detail=(
                    "final count differs from initial count but no causal "
                    "row-attrition records were emitted"
                ),
            )
        return None
    expected_input_count = initial_count
    for record in record_tuple:
        if not isinstance(record, RowAttritionRecord):
            _raise_row_attrition_error(
                workflow=workflow,
                stage="row_attrition_reconciliation",
                input_count=expected_input_count,
                output_count=final_count,
                detail=("row attrition records must be RowAttritionRecord instances"),
            )
        if record.input_rows != expected_input_count:
            _raise_row_attrition_error(
                workflow=workflow,
                stage=record.stage,
                input_count=int(record.input_rows),
                output_count=int(record.output_rows),
                expected_preceding_count=expected_input_count,
                detail="stage record is discontinuous with preceding output count",
            )
        expected_input_count = int(record.output_rows)
    if expected_input_count != final_count:
        last_record = record_tuple[-1]
        _raise_row_attrition_error(
            workflow=workflow,
            stage=last_record.stage,
            input_count=int(last_record.input_rows),
            output_count=int(last_record.output_rows),
            detail=(
                "accumulated stage records do not reconcile with final result count; "
                f"final_count={final_count}"
            ),
        )
    return RowAttritionReport.from_records(record_tuple)


def _site_id_tuple(values: pd.Index | Iterable[object]) -> tuple[str, ...]:
    if isinstance(values, pd.Index):
        raw_values = values.tolist()
    else:
        raw_values = list(values)
    return tuple(str(value) for value in raw_values)


def _examples(values: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(values[:_REMOVED_ID_SAMPLE_LIMIT])


def _raise_row_attrition_error(
    *,
    workflow: str,
    stage: str,
    input_count: int,
    output_count: int,
    detail: str,
    expected_preceding_count: int | None = None,
) -> None:
    expected_text = (
        ""
        if expected_preceding_count is None
        else f"; expected_preceding_count={int(expected_preceding_count)}"
    )
    raise WorkflowStageError(
        f"{workflow} row attrition internal consistency error; "
        f"stage={stage}; input_count={int(input_count)}; "
        f"output_count={int(output_count)}{expected_text}; {detail}"
    )
