"""Builder-facing adapter for the internal preprocessing subsystem."""

from __future__ import annotations

import pandas as pd

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
    resolve_dataset_total_protein_correction_policy,
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
from phospy.datasets.preprocessing.report_schema import (
    PreprocessingOperationRow,
    PreprocessingRowCountRow,
    dataframe_from_operation_rows,
    dataframe_from_row_count_rows,
)
from phospy.datasets.processing_state import (
    ComparisonState,
    DatasetProcessingState,
    MissingDataState,
    NormalisationState,
    SiteMatrixState,
    TotalProteinCorrectionState,
)
from phospy.transformations.models import IntensityScaleState

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
        "site_matrix_duplicate_site_policy",
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
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
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
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None = None,
) -> DatasetProcessingState:
    """Build compact dataset processing state from the resolved preprocessing plan."""

    comparison_pairs = (
        None
        if plan.comparison_pairs is None
        else tuple((str(left), str(right)) for left, right in plan.comparison_pairs)
    )
    resolved_total_policy = resolve_dataset_total_protein_correction_policy(
        plan.total_protein_correction_policy
    )
    total_correction_applied = (
        resolved_total_policy != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )
    correction_diagnostics = _resolve_total_correction_diagnostics(
        preprocessing_trace=preprocessing_trace
    )
    default_formula = (
        "log2_phospho - log2_total"
        if resolved_total_policy
        == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL
        else None
    )
    default_requires_log_scale: bool | None = bool(total_correction_applied)
    default_input_scale = (
        "log2"
        if resolved_total_policy
        == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL
        else None
    )
    default_output_scale = (
        "log2_ratio"
        if resolved_total_policy
        == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL
        else None
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
            policy=resolved_total_policy,
            applied=total_correction_applied,
            formula=_resolve_optional_string_diagnostic(
                correction_diagnostics,
                key="formula",
                default=default_formula,
            ),
            requires_log_scale=_resolve_optional_bool_diagnostic(
                correction_diagnostics,
                key="requires_log_scale",
                default=default_requires_log_scale,
            ),
            input_scale=_resolve_optional_string_diagnostic(
                correction_diagnostics,
                key="input_scale",
                default=default_input_scale,
            ),
            output_scale=_resolve_optional_string_diagnostic(
                correction_diagnostics,
                key="output_scale",
                default=default_output_scale,
            ),
            diagnostics=correction_diagnostics,
        ),
        site_matrix=SiteMatrixState(
            policy=plan.site_matrix_policy,
            constructed=(
                plan.site_matrix_policy
                == DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA
            ),
            missing_data_policy=plan.site_matrix_missing_data_policy,
            minimum_observed_values=plan.site_matrix_minimum_observed_values,
            duplicate_site_policy=plan.site_matrix_duplicate_site_policy,
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


def _resolve_total_correction_diagnostics(
    *,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> dict[str, object] | None:
    if preprocessing_trace is None:
        return None
    for item in preprocessing_trace:
        if item.stage == DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION:
            return dict(item.diagnostics)
    return None


def _resolve_optional_string_diagnostic(
    diagnostics: dict[str, object] | None,
    *,
    key: str,
    default: str | None,
) -> str | None:
    if diagnostics is None:
        return default
    value = diagnostics.get(key, default)
    if value is None:
        return None
    return str(value)


def _resolve_optional_bool_diagnostic(
    diagnostics: dict[str, object] | None,
    *,
    key: str,
    default: bool | None,
) -> bool | None:
    if diagnostics is None or key not in diagnostics:
        return default
    value = diagnostics.get(key)
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    return default


def _build_preprocessing_provenance_tables(
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
        row_count_rows.append(
            PreprocessingRowCountRow(
                stage=stage,
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
                stage=stage,
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
        return str(
            resolve_dataset_total_protein_correction_policy(
                plan.total_protein_correction_policy
            )
        )
    if stage == DATASET_PREPROCESSING_STAGE_SITE_MATRIX:
        return plan.site_matrix_policy
    if stage == DATASET_PREPROCESSING_STAGE_COMPARISONS:
        return plan.comparison_building_policy
    return "unsupported_stage"
