"""Builder-facing provenance table adaptation for preprocessing trace records."""

from __future__ import annotations

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
        canonical_stage_metadata = resolve_builder_provenance_stage_order(plan)
        for stage_metadata in canonical_stage_metadata:
            stage = stage_metadata.provenance_stage_key
            stage_label = stage_metadata.display_label
            record = trace_by_stage.get(stage)
            operation = stage_metadata.operation_name(plan)
            parameters = stage_metadata.serialize_parameters(plan)
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
