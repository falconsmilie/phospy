"""Builder-facing provenance table adaptation for preprocessing trace records."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.science.datasets.preprocessing.report_schema import (
    PreprocessingOperationRow,
    PreprocessingRowCountRow,
    dataframe_from_operation_rows,
    dataframe_from_row_count_rows,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    resolve_builder_provenance_stage_order,
)

_PREPROCESSING_INPUT_STAGE = "preprocessing_input"
_PREPROCESSING_COMPLETE_STAGE = "preprocessing_complete"


class PreprocessingProvenanceAdapter:
    """Build public-facing preprocessing provenance tables from execution trace."""

    def build_tables(
        self,
        *,
        plan: PreprocessingPlan,
        input_row_count: int,
        output_row_count: int,
        trace: tuple[PreprocessingStageExecution, ...],
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        trace_by_stage = {record.stage: record for record in trace}
        row_count_rows: list[PreprocessingRowCountRow] = [
            PreprocessingRowCountRow(
                stage=_PREPROCESSING_INPUT_STAGE,
                input_rows=input_row_count,
                output_rows=input_row_count,
                dropped_rows=0,
            )
        ]
        operation_rows: list[PreprocessingOperationRow] = []

        row_cursor = input_row_count
        step_order = 1
        emitted_stages: set[str] = set()
        canonical_stage_metadata = resolve_builder_provenance_stage_order(plan)
        for stage_metadata in canonical_stage_metadata:
            stage = stage_metadata.provenance_stage_key
            emitted_stages.add(stage)
            stage_label = stage_metadata.display_label
            record = trace_by_stage.get(stage)
            operation = stage_metadata.operation_name(plan)
            parameters = dict(stage_metadata.serialize_parameters(plan))
            if record is None:
                stage_input_rows = row_cursor
                stage_output_rows = row_cursor
                notes = "stage not scheduled in preprocessing plan"
            else:
                stage_input_rows = int(record.input_rows)
                stage_output_rows = int(record.output_rows)
                notes = (
                    "stage executed"
                    if record.notes is None
                    else str(record.notes).strip()
                )
                parameters = _with_execution_summary(
                    stage=stage,
                    base_parameters=parameters,
                    record=record,
                )
            row_cursor = stage_output_rows
            row_count_rows.append(
                PreprocessingRowCountRow(
                    stage=stage_label,
                    input_rows=stage_input_rows,
                    output_rows=stage_output_rows,
                    dropped_rows=(
                        max(stage_input_rows - stage_output_rows, 0)
                        if record is None
                        else int(max(record.dropped_row_count, 0))
                    ),
                )
            )
            operation_rows.append(
                PreprocessingOperationRow(
                    step_order=step_order,
                    stage=stage_label,
                    operation=operation,
                    parameters=parameters,
                    input_rows=stage_input_rows,
                    output_rows=stage_output_rows,
                    notes=notes,
                )
            )
            step_order += 1

        for record in trace:
            if record.stage in emitted_stages:
                continue
            row_count_rows.append(
                PreprocessingRowCountRow(
                    stage=record.stage,
                    input_rows=int(record.input_rows),
                    output_rows=int(record.output_rows),
                    dropped_rows=int(max(record.dropped_row_count, 0)),
                )
            )
            operation_rows.append(
                PreprocessingOperationRow(
                    step_order=step_order,
                    stage=record.stage,
                    operation=record.operation,
                    parameters=_with_execution_summary(
                        stage=record.stage,
                        base_parameters=dict(record.parameters),
                        record=record,
                    ),
                    input_rows=int(record.input_rows),
                    output_rows=int(record.output_rows),
                    notes=(
                        "stage executed"
                        if record.notes is None
                        else str(record.notes).strip()
                    ),
                )
            )
            step_order += 1

        row_count_rows.append(
            PreprocessingRowCountRow(
                stage=_PREPROCESSING_COMPLETE_STAGE,
                input_rows=output_row_count,
                output_rows=output_row_count,
                dropped_rows=0,
            )
        )

        row_counts = dataframe_from_row_count_rows(row_count_rows)
        operations = dataframe_from_operation_rows(operation_rows)
        return row_counts, operations


def _with_execution_summary(
    *,
    stage: str,
    base_parameters: dict[str, object],
    record: PreprocessingStageExecution,
) -> dict[str, object]:
    diagnostics = (
        {}
        if not isinstance(record.diagnostics, Mapping)
        else {str(key): value for key, value in record.diagnostics.items()}
    )
    input_rows = int(max(record.input_rows, 0))
    output_rows = int(max(record.output_rows, 0))
    dropped_row_count = int(max(record.dropped_row_count, 0))
    imputed_cell_count = int(max(record.imputed_cell_count, 0))
    diagnostics_imputed_cell_count = _resolve_int(
        diagnostics.get("imputed_cell_count"),
        default=0,
    )
    if diagnostics_imputed_cell_count > imputed_cell_count:
        imputed_cell_count = diagnostics_imputed_cell_count
    input_matrix_cells = max(
        int(record.input_shape[0]) * int(record.input_shape[1]),
        0,
    )
    dropped_sample_count = _resolve_int(
        diagnostics.get("dropped_column_count"),
        default=0,
    )
    summary: dict[str, object] = {
        "input_rows": input_rows,
        "output_rows": output_rows,
        "dropped_rows": dropped_row_count,
        "dropped_row_fraction_of_input": _safe_fraction(
            numerator=dropped_row_count,
            denominator=input_rows,
        ),
        "imputed_cell_count": imputed_cell_count,
        "imputed_row_count": int(len(record.imputed_row_ids)),
        "imputed_cell_fraction_of_input_matrix": _safe_fraction(
            numerator=imputed_cell_count,
            denominator=input_matrix_cells,
        ),
        "imputation_scope": _resolve_imputation_scope(
            stage=stage,
            base_parameters=base_parameters,
        ),
        "dropped_sample_count": dropped_sample_count,
        "diagnostic_keys": sorted(diagnostics.keys()),
        "diagnostic_summary": _extract_scalar_diagnostics(diagnostics),
    }
    parameters = dict(base_parameters)
    parameters["execution_summary"] = summary
    return parameters


def _resolve_imputation_scope(
    *,
    stage: str,
    base_parameters: Mapping[str, object],
) -> str | None:
    if stage != "missing_data":
        return None
    policy = str(base_parameters.get("missing_data_policy", "")).strip()
    if policy == "impute_row_median":
        return "per_row"
    if policy in {"impute_knn", "impute_minprob"}:
        return "global_matrix"
    if policy == "forbid":
        return "none"
    return None


def _extract_scalar_diagnostics(values: Mapping[str, object]) -> dict[str, object]:
    summary: dict[str, object] = {}
    for key, value in values.items():
        if value is None or isinstance(value, str | int | float | bool):
            summary[str(key)] = value
            continue
        if isinstance(value, tuple):
            if len(value) <= 10 and all(
                item is None or isinstance(item, str | int | float | bool)
                for item in value
            ):
                summary[str(key)] = list(value)
            continue
        if isinstance(value, list):
            if len(value) <= 10 and all(
                item is None or isinstance(item, str | int | float | bool)
                for item in value
            ):
                summary[str(key)] = list(value)
    return summary


def _resolve_int(value: object, *, default: int) -> int:
    if value is None or isinstance(value, bool):
        return int(default)
    if isinstance(value, int):
        return int(value)
    return int(default)


def _safe_fraction(*, numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return float(numerator) / float(denominator)
