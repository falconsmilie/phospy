"""Builder-facing adapter for the internal preprocessing subsystem."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
)
from phospy.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingState,
)
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.datasets.preprocessing.report_rows import (
    compose_stage_owned_report_tables,
)
from phospy.datasets.preprocessing.report_schema import (
    PreprocessingOperationRow,
    PreprocessingRowCountRow,
    dataframe_from_operation_rows,
    dataframe_from_row_count_rows,
)
from phospy.datasets.processing_state import (
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
    TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
    ComparisonState,
    DatasetProcessingState,
    MissingDataDiagnostics,
    MissingDataState,
    NormalisationState,
    RuvReadinessState,
    SiteMatrixState,
    SiteSequenceResolutionState,
    TotalProteinCorrectionDiagnostics,
    TotalProteinCorrectionState,
)
from phospy.errors.build import DatasetBuildError
from phospy.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    QuantitativeMeaning,
)

_PREPROCESSING_INPUT_STAGE = "preprocessing_input"
_PREPROCESSING_COMPLETE_STAGE = "preprocessing_complete"
_STAGE_LABEL_TO_PARAMETERS: dict[str, tuple[str, ...]] = {
    DATASET_PREPROCESSING_STAGE_NORMALISATION: (),
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION: (
        "site_sequence_resolution_enabled",
        "site_sequence_resolution_fasta_path",
        "site_sequence_resolution_mode",
        "site_sequence_resolution_flank_size",
        "site_sequence_resolution_accession_column",
        "site_sequence_resolution_site_column",
    ),
    DATASET_PREPROCESSING_STAGE_MISSING_DATA: (
        "missing_data_policy",
        "missing_data_min_observed_values",
        "missing_data_q",
        "missing_data_width",
        "missing_data_seed",
        "missing_data_k",
        "missing_data_distance",
        "missing_data_max_missing_fraction_per_row",
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
        report_tables = compose_stage_owned_report_tables(
            preprocessed_state.report_rows
        )
        return PreprocessedDatasetBuildTables(
            phospho=preprocessed_state.phospho,
            site_metadata=preprocessed_state.site_metadata,
            sample_metadata=preprocessed_state.sample_metadata,
            total=preprocessed_state.total,
            comparisons=preprocessed_state.comparisons,
            comparison_group_stats=report_tables.comparison_group_stats,
            comparison_pair_stats=report_tables.comparison_pair_stats,
            preprocessing_row_counts=row_counts,
            preprocessing_operations=operations,
            row_audit=report_tables.row_audit,
            preprocessing_trace=trace,
            duplicate_site_resolution=report_tables.duplicate_site_resolution,
            metadata_conflicts=report_tables.metadata_conflicts,
        )


def build_dataset_processing_state(
    *,
    plan: PreprocessingPlan,
    intensity_scale_state: IntensityScaleState,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None = None,
    final_phospho: pd.DataFrame | None = None,
    final_site_metadata: pd.DataFrame | None = None,
    final_sample_metadata: pd.DataFrame | None = None,
) -> DatasetProcessingState:
    """Build compact dataset processing state from the resolved preprocessing plan."""

    comparison_pairs = (
        None
        if plan.comparison_pairs is None
        else tuple((str(left), str(right)) for left, right in plan.comparison_pairs)
    )
    resolved_total_policy = plan.total_protein_correction_policy
    total_correction_applied = (
        resolved_total_policy != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )
    correction_diagnostics = _resolve_total_correction_diagnostics(
        preprocessing_trace=preprocessing_trace
    )
    missing_data_diagnostics = _resolve_missing_data_diagnostics(
        preprocessing_trace=preprocessing_trace
    )
    site_sequence_resolution_diagnostics = (
        _resolve_site_sequence_resolution_diagnostics(
            preprocessing_trace=preprocessing_trace
        )
    )
    intensity_scale_state = _resolve_quantitative_meaning_state(
        intensity_scale_state=intensity_scale_state,
        total_correction_policy=resolved_total_policy,
        correction_diagnostics=correction_diagnostics,
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
    quantitative_meaning = intensity_scale_state.quantity
    if quantitative_meaning is None:
        raise DatasetBuildError("intensity-scale state is missing quantitative meaning")
    default_quantitative_meaning = quantitative_meaning.value
    correction_diagnostics = _with_default_string_diagnostic(
        correction_diagnostics,
        key="quantitative_meaning",
        default=default_quantitative_meaning,
    )
    correction_diagnostics = _with_default_int_diagnostic(
        correction_diagnostics,
        key="diagnostics_schema_version",
        default=TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
    )
    if missing_data_diagnostics is not None:
        missing_data_diagnostics = _with_default_int_diagnostic(
            missing_data_diagnostics,
            key="diagnostics_schema_version",
            default=MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
        )
        missing_data_diagnostics = _with_default_string_diagnostic(
            missing_data_diagnostics,
            key="missing_data_policy",
            default=str(plan.missing_data_policy),
        )
    typed_correction_diagnostics = (
        None
        if correction_diagnostics is None
        else TotalProteinCorrectionDiagnostics.from_payload(
            correction_diagnostics,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics"
            ),
        )
    )
    typed_missing_data_diagnostics = (
        None
        if missing_data_diagnostics is None
        else MissingDataDiagnostics.from_payload(
            missing_data_diagnostics,
            field_name="dataset processing state missing_data.diagnostics",
        )
    )
    output_missing_cell_count = _resolve_optional_int_diagnostic(
        missing_data_diagnostics,
        key="output_missing_cell_count",
        default=0,
    )
    imputed_cell_count = _resolve_optional_int_diagnostic(
        missing_data_diagnostics,
        key="imputed_cell_count",
        default=0,
    )
    ruv_readiness = _resolve_ruv_readiness_state(
        plan=plan,
        final_phospho=final_phospho,
        final_site_metadata=final_site_metadata,
        final_sample_metadata=final_sample_metadata,
        missing_data_diagnostics=missing_data_diagnostics,
        default_matrix_complete=(output_missing_cell_count == 0),
    )
    return DatasetProcessingState(
        intensity_scale=intensity_scale_state,
        site_sequence_resolution=SiteSequenceResolutionState(
            configured=bool(plan.site_sequence_resolution_enabled),
            mode=(
                str(plan.site_sequence_resolution_mode).strip()
                if plan.site_sequence_resolution_enabled
                else None
            ),
            flank_size=(
                int(plan.site_sequence_resolution_flank_size)
                if plan.site_sequence_resolution_enabled
                else None
            ),
            fasta_sha256=_resolve_optional_string_diagnostic(
                site_sequence_resolution_diagnostics,
                key="fasta_sha256",
                default=None,
            ),
            resolved_site_count=_resolve_optional_int_diagnostic(
                site_sequence_resolution_diagnostics,
                key="resolved_site_count",
                default=0,
            ),
            unresolved_site_count=_resolve_optional_int_diagnostic(
                site_sequence_resolution_diagnostics,
                key="unresolved_site_count",
                default=0,
            ),
            unresolved_counts_by_reason=_resolve_optional_mapping_int_diagnostic(
                site_sequence_resolution_diagnostics,
                key="unresolved_counts_by_reason",
            ),
        ),
        missing_data=MissingDataState(
            policy=plan.missing_data_policy,
            min_observed_values=plan.missing_data_min_observed_values,
            complete_matrix=(output_missing_cell_count == 0),
            imputed=(imputed_cell_count > 0),
            diagnostics=typed_missing_data_diagnostics,
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
            quantitative_meaning=_resolve_optional_string_diagnostic(
                correction_diagnostics,
                key="quantitative_meaning",
                default=default_quantitative_meaning,
            ),
            diagnostics=typed_correction_diagnostics,
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
        ruv_readiness=ruv_readiness,
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


def _resolve_missing_data_diagnostics(
    *,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> dict[str, object] | None:
    if preprocessing_trace is None:
        return None
    for item in preprocessing_trace:
        if item.stage == DATASET_PREPROCESSING_STAGE_MISSING_DATA:
            return dict(item.diagnostics)
    return None


def _resolve_site_sequence_resolution_diagnostics(
    *,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> dict[str, object] | None:
    if preprocessing_trace is None:
        return None
    for item in preprocessing_trace:
        if item.stage == DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION:
            return dict(item.diagnostics)
    return None


def _resolve_optional_string_diagnostic(
    diagnostics: Mapping[str, object] | None,
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
    diagnostics: Mapping[str, object] | None,
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


def _resolve_optional_int_diagnostic(
    diagnostics: Mapping[str, object] | None,
    *,
    key: str,
    default: int,
) -> int:
    if diagnostics is None:
        return int(default)
    value = diagnostics.get(key)
    if value is None:
        return int(default)
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        stripped = value.strip()
        return int(default) if stripped == "" else int(stripped)
    return int(default)


def _resolve_optional_mapping_int_diagnostic(
    diagnostics: Mapping[str, object] | None,
    *,
    key: str,
) -> dict[str, int]:
    if diagnostics is None:
        return {}
    value = diagnostics.get(key)
    if not isinstance(value, Mapping):
        return {}
    resolved: dict[str, int] = {}
    for raw_key, raw_value in value.items():
        normalized_key = str(raw_key)
        if isinstance(raw_value, bool):
            resolved[normalized_key] = int(raw_value)
            continue
        if isinstance(raw_value, int):
            resolved[normalized_key] = int(raw_value)
            continue
        if isinstance(raw_value, float):
            resolved[normalized_key] = int(raw_value)
            continue
        if isinstance(raw_value, str):
            stripped = raw_value.strip()
            if stripped == "":
                continue
            resolved[normalized_key] = int(stripped)
            continue
    return resolved


def _with_default_string_diagnostic(
    diagnostics: dict[str, object] | None,
    *,
    key: str,
    default: str | None,
) -> dict[str, object] | None:
    if diagnostics is None:
        if default is None:
            return None
        return {key: default}
    resolved = dict(diagnostics)
    value = resolved.get(key)
    if value is None and default is not None:
        resolved[key] = default
    return resolved


def _with_default_int_diagnostic(
    diagnostics: dict[str, object] | None,
    *,
    key: str,
    default: int,
) -> dict[str, object] | None:
    if diagnostics is None:
        return None
    resolved = dict(diagnostics)
    value = resolved.get(key)
    if value is None:
        resolved[key] = default
    return resolved


def _resolve_quantitative_meaning_state(
    *,
    intensity_scale_state: IntensityScaleState,
    total_correction_policy: str,
    correction_diagnostics: Mapping[str, object] | None,
) -> IntensityScaleState:
    quantitative_meaning = _resolve_optional_string_diagnostic(
        correction_diagnostics,
        key="quantitative_meaning",
        default=None,
    )
    if quantitative_meaning is not None:
        try:
            return intensity_scale_state.with_quantitative_meaning(
                QuantitativeMeaning(quantitative_meaning)
            )
        except ValueError:
            pass
    if (
        total_correction_policy
        == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL
    ):
        return intensity_scale_state.with_quantitative_meaning(
            QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO
        )
    if intensity_scale_state.kind is IntensityScaleKind.LINEAR:
        return intensity_scale_state.with_quantitative_meaning(
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE
        )
    if intensity_scale_state.kind is IntensityScaleKind.LOG2:
        return intensity_scale_state.with_quantitative_meaning(
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE
        )
    return intensity_scale_state.with_quantitative_meaning(QuantitativeMeaning.UNKNOWN)


def _resolve_ruv_readiness_state(
    *,
    plan: PreprocessingPlan,
    final_phospho: pd.DataFrame | None,
    final_site_metadata: pd.DataFrame | None,
    final_sample_metadata: pd.DataFrame | None,
    missing_data_diagnostics: Mapping[str, object] | None,
    default_matrix_complete: bool,
) -> RuvReadinessState:
    enabled = bool(plan.ruv_readiness_enabled)
    matrix_complete = _resolve_matrix_completeness(
        final_phospho=final_phospho,
        default_matrix_complete=default_matrix_complete,
    )
    control_feature_column = str(plan.ruv_readiness_control_feature_column).strip()
    replicate_group_column = str(plan.ruv_readiness_replicate_group_column).strip()
    batch_column = plan.ruv_readiness_batch_column

    control_feature_count = _count_control_features(
        site_metadata=final_site_metadata,
        control_feature_column=control_feature_column,
    )
    replicate_group_count = _count_distinct_non_missing(
        sample_metadata=final_sample_metadata,
        column=replicate_group_column,
    )
    batch_count = (
        None
        if batch_column is None
        else _count_distinct_non_missing(
            sample_metadata=final_sample_metadata,
            column=batch_column,
        )
    )

    imputation_method_id = _resolve_optional_string_diagnostic(
        missing_data_diagnostics,
        key="imputation_method_id",
        default=None,
    )
    missingness_mask_hash = _resolve_optional_string_diagnostic(
        missing_data_diagnostics,
        key="missingness_mask_hash",
        default=None,
    )
    missingness_mask_preserved = missingness_mask_hash is not None
    reasons: list[str] = []
    if not enabled:
        reasons.append("not configured")
    else:
        if not matrix_complete:
            reasons.append("matrix contains missing values")
        if missing_data_diagnostics is None:
            reasons.append("missing-data diagnostics unavailable")
        if missingness_mask_hash is None:
            reasons.append("missingness_mask_hash unavailable")
        if final_site_metadata is None:
            reasons.append("site metadata unavailable")
        else:
            if control_feature_column not in final_site_metadata.columns:
                reasons.append("control feature column missing")
            if control_feature_count < 1:
                reasons.append("no control features present")
        if final_sample_metadata is None:
            reasons.append("sample metadata unavailable")
        else:
            if replicate_group_column not in final_sample_metadata.columns:
                reasons.append("replicate group column missing")
            if replicate_group_count < 2:
                reasons.append("insufficient replicate groups")
            if batch_column is not None:
                if batch_column not in final_sample_metadata.columns:
                    reasons.append("batch column missing")
                elif (batch_count or 0) < 1:
                    reasons.append("no batch values present")

    return RuvReadinessState(
        enabled=enabled,
        ready=enabled and len(reasons) == 0,
        reasons=tuple(reasons),
        control_feature_column=control_feature_column,
        replicate_group_column=replicate_group_column,
        batch_column=batch_column,
        control_feature_count=control_feature_count,
        replicate_group_count=replicate_group_count,
        batch_count=batch_count,
        requires_complete_matrix=True,
        matrix_complete=matrix_complete,
        imputation_method_id=imputation_method_id,
        missingness_mask_preserved=missingness_mask_preserved,
    )


def _resolve_matrix_completeness(
    *,
    final_phospho: pd.DataFrame | None,
    default_matrix_complete: bool,
) -> bool:
    if final_phospho is None:
        return bool(default_matrix_complete)
    return int(final_phospho.isna().to_numpy().sum()) == 0


def _count_control_features(
    *,
    site_metadata: pd.DataFrame | None,
    control_feature_column: str,
) -> int:
    if site_metadata is None or control_feature_column not in site_metadata.columns:
        return 0
    series = site_metadata.loc[:, control_feature_column]
    return int(series.map(_is_truthy_control_feature).sum())


def _is_truthy_control_feature(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        return normalized in {"1", "true", "t", "yes", "y"}
    return False


def _count_distinct_non_missing(
    *,
    sample_metadata: pd.DataFrame | None,
    column: str,
) -> int:
    if sample_metadata is None or column not in sample_metadata.columns:
        return 0
    values = (
        sample_metadata.loc[:, column]
        .astype("string")
        .str.strip()
        .replace("", pd.NA)
        .dropna()
    )
    return int(values.nunique())


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
    canonical_stages = _resolve_provenance_canonical_stages(plan)
    for stage in canonical_stages:
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


def _resolve_provenance_canonical_stages(plan: PreprocessingPlan) -> tuple[str, ...]:
    if not plan.site_sequence_resolution_enabled:
        return _PROVENANCE_CANONICAL_STAGES
    return (
        DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        *_PROVENANCE_CANONICAL_STAGES,
    )


def _resolve_stage_parameters(
    *, plan: PreprocessingPlan, stage: str
) -> dict[str, object]:
    if stage == DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION:
        return {
            "enabled": bool(plan.site_sequence_resolution_enabled),
            "fasta_path": plan.site_sequence_resolution_fasta_path,
            "mode": plan.site_sequence_resolution_mode,
            "flank_size": int(plan.site_sequence_resolution_flank_size),
            "accession_column": plan.site_sequence_resolution_accession_column,
            "site_column": plan.site_sequence_resolution_site_column,
        }
    if stage == DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION:
        identity = plan.total_protein_correction_identity_policy
        return {
            "total_protein_correction_policy": plan.total_protein_correction_policy,
            "identity_mode": identity.mode,
            "phosphosite_key": identity.phosphosite_key,
            "total_protein_key": identity.total_protein_key,
            "mapping_phosphosite_key": identity.mapping_phosphosite_key,
            "mapping_total_protein_key": identity.mapping_total_protein_key,
            "mapping_table_fingerprint": identity.mapping_table_fingerprint,
            "mapping_table_row_count": (
                None if identity.mapping_table is None else len(identity.mapping_table)
            ),
            "duplicate_policy": identity.duplicate_policy,
            "unmatched_policy": identity.unmatched_policy,
        }
    if stage == DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM:
        return {"pseudocount": float(plan.intensity_transform_pseudocount)}
    parameter_names = _STAGE_LABEL_TO_PARAMETERS.get(stage, ())
    return {name: getattr(plan, name) for name in parameter_names}


def _resolve_stage_operation(*, plan: PreprocessingPlan, stage: str) -> str:
    if stage == DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION:
        return plan.site_sequence_resolution_mode
    if stage == DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM:
        return plan.intensity_transform_policy
    if stage == DATASET_PREPROCESSING_STAGE_NORMALISATION:
        return plan.normalisation_policy
    if stage == DATASET_PREPROCESSING_STAGE_MISSING_DATA:
        return plan.missing_data_policy
    if stage == DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION:
        return str(plan.total_protein_correction_policy)
    if stage == DATASET_PREPROCESSING_STAGE_SITE_MATRIX:
        return plan.site_matrix_policy
    if stage == DATASET_PREPROCESSING_STAGE_COMPARISONS:
        return plan.comparison_building_policy
    return "unsupported_stage"
