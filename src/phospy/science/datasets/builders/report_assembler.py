"""Preprocessing report assembly for dataset builder executor."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.science.datasets.models import (
    DatasetPreprocessingReport,
    SiteSequenceResolutionReport,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
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
        return DatasetPreprocessingReport.from_rows(
            row_count_rows=tuple(row_count_rows),
            operation_rows=tuple(operation_rows),
            row_audit_rows=row_audit_rows,
            duplicate_site_resolution_rows=duplicate_site_resolution_rows,
            metadata_conflict_rows=metadata_conflict_rows,
            comparison_group_stats_rows=comparison_group_stats_rows,
            comparison_pair_stats_rows=comparison_pair_stats_rows,
            site_sequence_resolution=site_sequence_resolution,
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
