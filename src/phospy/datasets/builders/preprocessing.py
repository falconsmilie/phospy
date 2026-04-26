"""Builder-facing adapter for the internal preprocessing subsystem."""

from __future__ import annotations

import pandas as pd

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
)
from phospy.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingState,
    empty_preprocessing_row_audit,
)
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.datasets.processing_state import (
    ComparisonState,
    DatasetProcessingState,
    MissingDataState,
    NormalisationState,
    SiteMatrixState,
    TotalProteinCorrectionState,
)
from phospy.transformations.models import IntensityScaleState

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
_STAGE_LABEL_TO_PARAMETERS: dict[str, tuple[str, ...]] = {
    DATASET_PREPROCESSING_STAGE_NORMALISATION: (),
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
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
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
        row_audit = (
            empty_preprocessing_row_audit()
            if preprocessed_state.row_audit is None
            else preprocessed_state.row_audit
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
            row_audit=row_audit,
            preprocessing_trace=trace,
            duplicate_site_resolution=preprocessed_state.duplicate_site_resolution,
            metadata_conflicts=preprocessed_state.metadata_conflicts,
        )


def build_dataset_processing_state(
    *,
    plan: PreprocessingPlan,
    intensity_scale_state: IntensityScaleState,
) -> DatasetProcessingState:
    """Build compact dataset processing state from the resolved preprocessing plan."""

    comparison_pairs = (
        None
        if plan.comparison_pairs is None
        else tuple((str(left), str(right)) for left, right in plan.comparison_pairs)
    )
    return DatasetProcessingState(
        intensity_scale=intensity_scale_state,
        missing_data=MissingDataState(
            policy=plan.missing_data_policy,
            min_observed_values=plan.missing_data_min_observed_values,
            complete_matrix=True,
            imputed=(
                plan.missing_data_policy
                == DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN
            ),
        ),
        normalisation=NormalisationState(policy=plan.normalisation_policy),
        total_protein_correction=TotalProteinCorrectionState(
            policy=plan.total_protein_correction_policy,
            applied=(
                plan.total_protein_correction_policy
                != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
            ),
        ),
        site_matrix=SiteMatrixState(
            policy=plan.site_matrix_policy,
            constructed=(
                plan.site_matrix_policy
                == DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA
            ),
            missing_data_policy=plan.site_matrix_missing_data_policy,
            minimum_observed_values=plan.site_matrix_minimum_observed_values,
            duplicate_site_strategy=plan.site_matrix_duplicate_site_strategy,
        ),
        comparisons=ComparisonState(
            policy=plan.comparison_building_policy,
            sample_group_column=plan.comparison_sample_group_column,
            pairs=(
                None
                if plan.comparison_building_policy
                == DATASET_COMPARISON_BUILDING_POLICY_NONE
                else comparison_pairs
            ),
        ),
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
            operation = _resolve_stage_operation(plan=plan, stage=stage)
            parameters = _resolve_stage_parameters(plan=plan, stage=stage)
        else:
            stage_input_rows = int(record.input_rows)
            stage_output_rows = int(record.output_rows)
            notes = (
                "stage executed" if record.notes is None else str(record.notes).strip()
            )
            operation = record.operation
            parameters = dict(record.parameters)
        row_cursor = stage_output_rows
        row_counts_records.append(
            {
                "stage": stage,
                "input_rows": stage_input_rows,
                "output_rows": stage_output_rows,
                "dropped_rows": (
                    max(stage_input_rows - stage_output_rows, 0)
                    if record is None
                    else int(max(record.dropped_row_count, 0))
                ),
            }
        )
        operations_records.append(
            {
                "step_order": step_order,
                "stage": stage,
                "operation": operation,
                "parameters": parameters,
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
    if stage == DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM:
        return {"pseudocount": float(plan.intensity_transform_pseudocount)}
    parameter_names = _STAGE_LABEL_TO_PARAMETERS.get(stage, ())
    return {name: getattr(plan, name) for name in parameter_names}


def _resolve_stage_operation(*, plan: PreprocessingPlan, stage: str) -> str:
    if stage == DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM:
        return plan.intensity_transform_policy
    if stage == DATASET_PREPROCESSING_STAGE_NORMALISATION:
        return plan.normalisation_policy
    if stage == DATASET_PREPROCESSING_STAGE_MISSING_DATA:
        return plan.missing_data_policy
    if stage == DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION:
        return plan.total_protein_correction_policy
    if stage == DATASET_PREPROCESSING_STAGE_SITE_MATRIX:
        return plan.site_matrix_policy
    if stage == DATASET_PREPROCESSING_STAGE_COMPARISONS:
        return plan.comparison_building_policy
    return "unsupported_stage"
