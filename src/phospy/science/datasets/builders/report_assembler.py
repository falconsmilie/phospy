"""Preprocessing report assembly for dataset builder executor."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.contracts.configs import DATASET_BATCH_CORRECTION_METHOD_NONE
from phospy.science.datasets.models import (
    DatasetPreprocessingReport,
    SiteSequenceResolutionReport,
)
from phospy.science.datasets.preprocessing.batch_correction import (
    BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
    BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED,
    BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS,
    BATCH_CORRECTION_STATUS_DISABLED,
    BATCH_CORRECTION_STATUS_REJECTED,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.science.datasets.preprocessing.report_schema import (
    PreprocessingOperationRow,
    PreprocessingRowCountRow,
    comparison_group_stats_rows_from_dataframe,
    comparison_pair_stats_rows_from_dataframe,
    duplicate_site_resolution_rows_from_dataframe,
    metadata_conflict_rows_from_dataframe,
    operation_rows_from_dataframe,
    row_audit_rows_from_dataframe,
    row_count_rows_from_dataframe,
)

_FINAL_DATASET_STAGE = "final_dataset_construction"
_PEPTIDE_EVIDENCE_RESOLUTION_STAGE = "peptide_evidence_resolution"
_SITE_SEQUENCE_RESOLUTION_STAGE = DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION
_SITE_SEQUENCE_CONFLICT_POLICY_NOT_APPLIED = "not_applied"


class DatasetPreprocessingReportAssembler:
    """Assemble the final dataset preprocessing report payload."""

    def run(
        self,
        *,
        row_counts: pd.DataFrame | None,
        operations: pd.DataFrame | None,
        row_audit: pd.DataFrame | None,
        duplicate_site_resolution: pd.DataFrame | None,
        metadata_conflicts: pd.DataFrame | None,
        comparison_group_stats: pd.DataFrame | None,
        comparison_pair_stats: pd.DataFrame | None,
        preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
        site_sequence_derivation: dict[str, object] | None,
        input_site_count: int,
        final_dataset_rows: int,
        intensity_scale_label: str,
        intensity_scale_establishment: Mapping[str, object],
        quantitative_meaning: str,
        peptide_evidence_resolution: dict[str, object] | None,
        preprocessing_plan: PreprocessingPlan | None = None,
        sample_metadata: pd.DataFrame | None = None,
        matrix_shape_before: tuple[int, int] | None = None,
        matrix_shape_after: tuple[int, int] | None = None,
        declared_input_intensity_scale_kind: str | None = None,
    ) -> DatasetPreprocessingReport:
        row_count_rows = list(row_count_rows_from_dataframe(row_counts))
        operation_rows = list(operation_rows_from_dataframe(operations))
        row_audit_rows = row_audit_rows_from_dataframe(row_audit)
        duplicate_site_resolution_rows = duplicate_site_resolution_rows_from_dataframe(
            duplicate_site_resolution
        )
        metadata_conflict_rows = metadata_conflict_rows_from_dataframe(
            metadata_conflicts
        )
        comparison_group_stats_rows = comparison_group_stats_rows_from_dataframe(
            comparison_group_stats
        )
        comparison_pair_stats_rows = comparison_pair_stats_rows_from_dataframe(
            comparison_pair_stats
        )
        if peptide_evidence_resolution is not None:
            peptide_observations = _coerce_non_negative_int(
                peptide_evidence_resolution.get("peptide_observations_received"),
                default=0,
            )
            unique_site_ids = _coerce_non_negative_int(
                peptide_evidence_resolution.get("unique_site_ids_produced"),
                default=0,
            )
            excluded_observations = _coerce_non_negative_int(
                peptide_evidence_resolution.get("excluded_observations"),
                default=max(peptide_observations - unique_site_ids, 0),
            )
            row_count_rows.append(
                PreprocessingRowCountRow(
                    stage=_PEPTIDE_EVIDENCE_RESOLUTION_STAGE,
                    input_rows=peptide_observations,
                    output_rows=unique_site_ids,
                    dropped_rows=excluded_observations,
                )
            )
            if not operation_rows:
                step_order = 1
            else:
                step_order = int(max(row.step_order for row in operation_rows)) + 1
            operation_rows.append(
                PreprocessingOperationRow(
                    step_order=step_order,
                    stage=_PEPTIDE_EVIDENCE_RESOLUTION_STAGE,
                    operation="resolve_peptide_evidence_to_site_level",
                    parameters=dict(peptide_evidence_resolution),
                    input_rows=peptide_observations,
                    output_rows=unique_site_ids,
                    notes=(
                        "peptide evidence ambiguity policy and site-resolution summary"
                    ),
                )
            )

        row_count_rows.append(
            PreprocessingRowCountRow(
                stage=_FINAL_DATASET_STAGE,
                input_rows=final_dataset_rows,
                output_rows=final_dataset_rows,
                dropped_rows=0,
            )
        )
        if not operation_rows:
            final_step_order = 1
        else:
            final_step_order = int(max(row.step_order for row in operation_rows)) + 1
        operation_rows.append(
            PreprocessingOperationRow(
                step_order=final_step_order,
                stage=_FINAL_DATASET_STAGE,
                operation="construct_analysis_ready_dataset",
                parameters={
                    "intensity_scale_label": intensity_scale_label,
                    "intensity_scale_establishment": dict(
                        intensity_scale_establishment
                    ),
                    "intensity_transformation_state": {
                        "before_preprocessing": _resolve_input_intensity_scale_label(
                            declared_input_intensity_scale_kind=(
                                declared_input_intensity_scale_kind
                            ),
                            preprocessing_trace=preprocessing_trace,
                            intensity_scale_establishment=intensity_scale_establishment,
                        ),
                        "after_preprocessing": str(intensity_scale_label).strip(),
                    },
                    "preprocessing_diagnostics": _summarize_preprocessing_diagnostics(
                        preprocessing_trace
                    ),
                    "quantitative_meaning": quantitative_meaning,
                },
                input_rows=final_dataset_rows,
                output_rows=final_dataset_rows,
                notes="analysis-ready dataset boundary construction",
            )
        )
        site_sequence_resolution = _build_site_sequence_resolution_report(
            preprocessing_trace=preprocessing_trace,
            site_sequence_derivation=site_sequence_derivation,
            total_sites=int(input_site_count),
            final_sequence_complete_sites=int(final_dataset_rows),
        )
        batch_correction = _build_batch_correction_report(
            plan=preprocessing_plan,
            sample_metadata=sample_metadata,
            matrix_shape_before=matrix_shape_before,
            matrix_shape_after=matrix_shape_after,
        )
        return DatasetPreprocessingReport.from_rows(
            row_count_rows=tuple(row_count_rows),
            operation_rows=tuple(operation_rows),
            row_audit_rows=row_audit_rows,
            duplicate_site_resolution_rows=duplicate_site_resolution_rows,
            metadata_conflict_rows=metadata_conflict_rows,
            comparison_group_stats_rows=comparison_group_stats_rows,
            comparison_pair_stats_rows=comparison_pair_stats_rows,
            site_sequence_resolution=site_sequence_resolution,
            batch_correction=batch_correction,
        )


def _build_site_sequence_resolution_report(
    *,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    site_sequence_derivation: dict[str, object] | None,
    total_sites: int,
    final_sequence_complete_sites: int,
) -> SiteSequenceResolutionReport:
    stage_diagnostics = _resolve_site_sequence_stage_diagnostics(preprocessing_trace)
    if stage_diagnostics is not None:
        (
            provided_by_input,
            resolved_from_fasta,
            unresolved,
        ) = _summarize_stage_sequence_origins(stage_diagnostics)
        conflicts = _coerce_non_negative_int(
            stage_diagnostics.get("existing_sequence_conflict_count"),
            default=0,
        )
        conflict_policy = _resolve_conflict_policy(
            stage_diagnostics.get("conflict_policy")
        )
        return SiteSequenceResolutionReport(
            total_sites=int(max(total_sites, 0)),
            provided_by_input=int(max(provided_by_input, 0)),
            resolved_from_fasta=int(max(resolved_from_fasta, 0)),
            resolved_from_reference=0,
            unresolved=int(max(unresolved, 0)),
            conflicts=int(max(conflicts, 0)),
            conflict_policy=conflict_policy,
            final_sequence_complete_sites=int(max(final_sequence_complete_sites, 0)),
        )

    derivation = (
        {}
        if not isinstance(site_sequence_derivation, Mapping)
        else site_sequence_derivation
    )
    provided_by_input = _coerce_non_negative_int(
        derivation.get("provided_sequence_count"),
        default=0,
    )
    resolved_from_reference = _coerce_non_negative_int(
        derivation.get("derived_sequence_count"),
        default=0,
    )
    unresolved = _coerce_non_negative_int(
        derivation.get("unresolved_sequence_count"),
        default=max(total_sites - provided_by_input - resolved_from_reference, 0),
    )
    conflicts = _coerce_non_negative_int(
        derivation.get("existing_sequence_conflict_count"),
        default=0,
    )
    return SiteSequenceResolutionReport(
        total_sites=int(max(total_sites, 0)),
        provided_by_input=int(max(provided_by_input, 0)),
        resolved_from_fasta=0,
        resolved_from_reference=int(max(resolved_from_reference, 0)),
        unresolved=int(max(unresolved, 0)),
        conflicts=int(max(conflicts, 0)),
        conflict_policy=_SITE_SEQUENCE_CONFLICT_POLICY_NOT_APPLIED,
        final_sequence_complete_sites=int(max(final_sequence_complete_sites, 0)),
    )


def _build_batch_correction_report(
    *,
    plan: PreprocessingPlan | None,
    sample_metadata: pd.DataFrame | None,
    matrix_shape_before: tuple[int, int] | None,
    matrix_shape_after: tuple[int, int] | None,
) -> BatchCorrectionReport | None:
    if plan is None:
        return None
    method = str(plan.batch_correction_method).strip()
    if not method:
        method = DATASET_BATCH_CORRECTION_METHOD_NONE
    batch_column = str(plan.batch_correction_batch_column).strip()
    condition_column = str(plan.batch_correction_condition_column).strip()
    has_batch_column = _sample_metadata_has_column(sample_metadata, batch_column)
    batch_levels = _sample_metadata_levels(sample_metadata, batch_column)
    condition_levels = _sample_metadata_levels(sample_metadata, condition_column)
    number_of_batches = len(batch_levels) if has_batch_column else None
    if method == DATASET_BATCH_CORRECTION_METHOD_NONE:
        status = BATCH_CORRECTION_STATUS_DISABLED
        confounding_status = BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE
        warnings: tuple[str, ...] = ()
        limitations = ("batch correction disabled by preprocessing configuration",)
    else:
        status = BATCH_CORRECTION_STATUS_REJECTED
        confounding_status = BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED
        warnings = ("batch correction was declared but no correction was executed",)
        limitations = (
            "batch correction execution is not implemented yet; matrix values "
            "are unchanged",
        )
    no_op_matrix_shape = (
        matrix_shape_after if matrix_shape_after is not None else matrix_shape_before
    )
    return BatchCorrectionReport(
        status=status,
        policy=BatchCorrectionPolicy(
            method=method,
            batch_column=batch_column,
            condition_column=condition_column,
            design_preservation_policy=(
                BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS
            ),
            preserve_condition_effects=bool(
                plan.batch_correction_preserve_condition_effects
            ),
        ),
        diagnostics=BatchCorrectionDiagnostics(
            number_of_batches=number_of_batches,
            batch_levels=batch_levels,
            condition_levels=condition_levels,
            confounding_check_status=confounding_status,
            matrix_shape_before=no_op_matrix_shape,
            matrix_shape_after=no_op_matrix_shape,
            warnings=warnings,
            limitations=limitations,
        ),
    )


def _sample_metadata_has_column(
    sample_metadata: pd.DataFrame | None,
    column: str | None,
) -> bool:
    if sample_metadata is None:
        return False
    if column is None or str(column).strip() == "":
        return False
    normalized_column = str(column).strip()
    return (
        sum(
            1
            for sample_column in sample_metadata.columns
            if sample_column == normalized_column
        )
        == 1
    )


def _sample_metadata_levels(
    sample_metadata: pd.DataFrame | None,
    column: str | None,
) -> tuple[str, ...]:
    if not _sample_metadata_has_column(sample_metadata, column):
        return ()
    assert sample_metadata is not None
    assert column is not None
    levels: list[str] = []
    seen: set[str] = set()
    for value in sample_metadata.loc[:, str(column).strip()].tolist():
        if _is_missing_metadata_value(value):
            continue
        level = str(value).strip()
        if level == "" or level in seen:
            continue
        seen.add(level)
        levels.append(level)
    return tuple(levels)


def _is_missing_metadata_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _resolve_site_sequence_stage_diagnostics(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> dict[str, object] | None:
    if preprocessing_trace is None:
        return None
    for stage in preprocessing_trace:
        if stage.stage == _SITE_SEQUENCE_RESOLUTION_STAGE:
            return (
                {}
                if not isinstance(stage.diagnostics, Mapping)
                else dict(stage.diagnostics)
            )
    return None


def _summarize_stage_sequence_origins(
    stage_diagnostics: Mapping[str, object],
) -> tuple[int, int, int]:
    row_diagnostics = stage_diagnostics.get("row_diagnostics")
    if not isinstance(row_diagnostics, list):
        provided_by_input = _coerce_non_negative_int(
            stage_diagnostics.get("preserved_existing_count"),
            default=0,
        )
        resolved_from_fasta = _coerce_non_negative_int(
            stage_diagnostics.get("filled_missing_count"),
            default=0,
        ) + _coerce_non_negative_int(
            stage_diagnostics.get("replaced_existing_count"),
            default=0,
        )
        unresolved = _coerce_non_negative_int(
            stage_diagnostics.get("unresolved_site_count"),
            default=0,
        )
        return (provided_by_input, resolved_from_fasta, unresolved)

    provided_by_input = 0
    resolved_from_fasta = 0
    unresolved = 0
    for row in row_diagnostics:
        if not isinstance(row, Mapping):
            continue
        action = str(row.get("action", "")).strip().lower()
        existing_site_sequence = row.get("existing_site_sequence")
        resolved_site_sequence = row.get("resolved_site_sequence")
        has_existing = _has_resolved_site_sequence(existing_site_sequence)
        has_resolved = _has_resolved_site_sequence(resolved_site_sequence)

        if action in {"fill_missing", "replace_existing"} and has_resolved:
            resolved_from_fasta += 1
            continue
        if not has_resolved:
            unresolved += 1
            continue
        if has_existing:
            provided_by_input += 1
            continue
        if action in {"validate_existing", "preserve_existing"}:
            provided_by_input += 1
            continue
        resolved_from_fasta += 1
    return (provided_by_input, resolved_from_fasta, unresolved)


def _resolve_conflict_policy(value: object) -> str:
    if isinstance(value, str):
        normalized = value.strip()
        if normalized:
            return normalized
    return _SITE_SEQUENCE_CONFLICT_POLICY_NOT_APPLIED


def _has_resolved_site_sequence(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    if not text:
        return False
    return text.lower() != "none"


def _coerce_non_negative_int(value: object, *, default: int) -> int:
    if isinstance(value, bool):
        return int(default)
    if isinstance(value, int):
        return max(int(value), 0)
    return int(max(default, 0))


def _resolve_input_intensity_scale_label(
    *,
    declared_input_intensity_scale_kind: str | None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    intensity_scale_establishment: Mapping[str, object],
) -> str | None:
    declared = (
        None
        if declared_input_intensity_scale_kind is None
        else str(declared_input_intensity_scale_kind).strip()
    )
    if declared:
        return declared
    establishment_parameters = intensity_scale_establishment.get("parameters")
    if isinstance(establishment_parameters, Mapping):
        declared_parameter = establishment_parameters.get("declared_scale_kind")
        if isinstance(declared_parameter, str) and declared_parameter.strip():
            return declared_parameter.strip()
    intensity_transform_stage = _resolve_stage(
        preprocessing_trace=preprocessing_trace,
        stage_key=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    )
    if intensity_transform_stage is None:
        return None
    if str(intensity_transform_stage.operation).strip() == "log2":
        return "linear_or_unknown"
    return "unknown"


def _summarize_preprocessing_diagnostics(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> dict[str, object]:
    if preprocessing_trace is None:
        return {"warnings": [], "stages": []}
    warnings: list[str] = []
    stage_summaries: list[dict[str, object]] = []
    for stage in preprocessing_trace:
        diagnostics = (
            {}
            if not isinstance(stage.diagnostics, Mapping)
            else dict(stage.diagnostics)
        )
        stage_warning_messages = _extract_warning_messages(diagnostics)
        warnings.extend(stage_warning_messages)
        stage_summaries.append(
            {
                "stage": stage.stage,
                "operation": stage.operation,
                "dropped_row_count": int(max(stage.dropped_row_count, 0)),
                "imputed_cell_count": int(max(stage.imputed_cell_count, 0)),
                "diagnostic_keys": sorted(str(key) for key in diagnostics.keys()),
                "warning_count": len(stage_warning_messages),
            }
        )
    return {"warnings": list(dict.fromkeys(warnings)), "stages": stage_summaries}


def _extract_warning_messages(diagnostics: Mapping[str, object]) -> list[str]:
    messages: list[str] = []
    for key in ("diagnostic_warnings", "warnings"):
        value = diagnostics.get(key)
        if isinstance(value, list | tuple):
            for item in value:
                if isinstance(item, str) and item.strip():
                    messages.append(item.strip())
    note = diagnostics.get("note")
    if isinstance(note, str) and note.strip():
        messages.append(note.strip())
    return messages


def _resolve_stage(
    *,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    stage_key: str,
) -> PreprocessingStageExecution | None:
    if preprocessing_trace is None:
        return None
    for stage in preprocessing_trace:
        if stage.stage == stage_key:
            return stage
    return None
