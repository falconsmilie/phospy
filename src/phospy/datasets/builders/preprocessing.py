"""Builder-facing adapter for the internal preprocessing subsystem."""

from __future__ import annotations

import pandas as pd

from phospy.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingState,
)
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline

_ROW_COUNT_COLUMNS = ("stage", "input_rows", "output_rows", "dropped_rows")
_OPERATION_COLUMNS = (
    "step_order",
    "stage",
    "operation",
    "parameters",
    "input_rows",
    "output_rows",
    "notes",
)
_PREPROCESSING_INPUT_STAGE = "preprocessing_input"
_PREPROCESSING_COMPLETE_STAGE = "preprocessing_complete"
_STAGE_LABEL_TO_OPERATION = {
    DATASET_PREPROCESSING_STAGE_MISSING_DATA: "apply_missing_data_policy",
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION: "apply_total_protein_correction",
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX: "apply_site_matrix_policy",
    DATASET_PREPROCESSING_STAGE_COMPARISONS: "build_comparisons",
}
_STAGE_LABEL_TO_PARAMETERS: dict[str, tuple[str, ...]] = {
    DATASET_PREPROCESSING_STAGE_MISSING_DATA: (
        "missing_data_policy",
        "missing_data_min_observed_values",
    ),
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION: (
        "total_protein_correction_policy",
    ),
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX: (
        "site_matrix_policy",
        "site_matrix_duplicate_site_strategy",
        "site_matrix_missing_data_policy",
        "site_matrix_minimum_observed_values",
    ),
    DATASET_PREPROCESSING_STAGE_COMPARISONS: (
        "comparison_building_policy",
        "comparison_sample_group_column",
        "comparison_pairs",
    ),
}
_PROVENANCE_CANONICAL_STAGES = (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
)


class DatasetPreprocessor:
    """Translate builder input into internal preprocessing state and run it."""

    def __init__(self, *, pipeline: PreprocessingPipeline | None = None) -> None:
        self._pipeline = pipeline or PreprocessingPipeline()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        plan: PreprocessingPlan,
    ) -> PreprocessedDatasetBuildTables:
        input_row_count = int(len(phospho.index))
        preprocessed_state, trace = self._pipeline.run_with_trace(
            PreprocessingState(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                plan=plan,
            )
        )
        row_counts, operations = _build_preprocessing_provenance_tables(
            plan=plan,
            input_row_count=input_row_count,
            output_row_count=int(len(preprocessed_state.phospho.index)),
            trace=trace,
        )
        return PreprocessedDatasetBuildTables(
            phospho=preprocessed_state.phospho,
            site_metadata=preprocessed_state.site_metadata,
            sample_metadata=preprocessed_state.sample_metadata,
            total=preprocessed_state.total,
            comparisons=preprocessed_state.comparisons,
            comparison_group_stats=preprocessed_state.comparison_group_stats,
            comparison_pair_stats=preprocessed_state.comparison_pair_stats,
            preprocessing_row_counts=row_counts,
            preprocessing_operations=operations,
            duplicate_site_resolution=preprocessed_state.duplicate_site_resolution,
            metadata_conflicts=preprocessed_state.metadata_conflicts,
        )


def _build_preprocessing_provenance_tables(
    *,
    plan: PreprocessingPlan,
    input_row_count: int,
    output_row_count: int,
    trace: tuple[PreprocessingStageExecution, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    trace_by_stage = {record.stage: record for record in trace}
    row_counts_records: list[dict[str, int | str]] = [
        {
            "stage": _PREPROCESSING_INPUT_STAGE,
            "input_rows": input_row_count,
            "output_rows": input_row_count,
            "dropped_rows": 0,
        }
    ]
    operations_records: list[dict[str, object]] = []

    row_cursor = input_row_count
    step_order = 1
    for stage in _PROVENANCE_CANONICAL_STAGES:
        record = trace_by_stage.get(stage)
        if record is None:
            stage_input_rows = row_cursor
            stage_output_rows = row_cursor
            notes = "stage not scheduled in preprocessing plan"
        else:
            stage_input_rows = int(record.input_rows)
            stage_output_rows = int(record.output_rows)
            notes = "stage executed"
        row_cursor = stage_output_rows
        row_counts_records.append(
            {
                "stage": stage,
                "input_rows": stage_input_rows,
                "output_rows": stage_output_rows,
                "dropped_rows": max(stage_input_rows - stage_output_rows, 0),
            }
        )
        operations_records.append(
            {
                "step_order": step_order,
                "stage": stage,
                "operation": _STAGE_LABEL_TO_OPERATION[stage],
                "parameters": _resolve_stage_parameters(plan=plan, stage=stage),
                "input_rows": stage_input_rows,
                "output_rows": stage_output_rows,
                "notes": notes,
            }
        )
        step_order += 1

    row_counts_records.append(
        {
            "stage": _PREPROCESSING_COMPLETE_STAGE,
            "input_rows": output_row_count,
            "output_rows": output_row_count,
            "dropped_rows": 0,
        }
    )

    row_counts = pd.DataFrame.from_records(
        row_counts_records, columns=_ROW_COUNT_COLUMNS
    )
    operations = pd.DataFrame.from_records(
        operations_records,
        columns=_OPERATION_COLUMNS,
    )
    return row_counts, operations


def _resolve_stage_parameters(
    *, plan: PreprocessingPlan, stage: str
) -> dict[str, object]:
    parameter_names = _STAGE_LABEL_TO_PARAMETERS.get(stage, ())
    return {name: getattr(plan, name) for name in parameter_names}
