"""Internal executor for the dataset builder path.

The public builder lane stays intentionally narrow: establish supported
intensity scale state after applying explicit builder preprocessing policy.
"""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.datasets.builders.contracts import (
    DatasetPreprocessorContract,
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.datasets.builders.preprocessing import (
    DatasetPreprocessor,
    build_dataset_processing_state,
)
from phospy.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    TotalProteinCorrectionIdentityPolicy,
)
from phospy.datasets.preprocessing.report_schema import (
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
from phospy.errors.build import DatasetBuildError
from phospy.errors.transformations import (
    TransformationStateEstablishmentError,
)
from phospy.policy_models import IntensityTransformPolicy
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    JsonValue,
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.scientific_policies import (
    PreprocessingStageOrderPolicy,
    ScientificPolicyRecord,
    build_duplicate_site_resolution_policy,
)
from phospy.transformations.contracts import Transformer
from phospy.transformations.models import IntensityScaleKind
from phospy.transformations.transformers import IdentityTransformer

_FINAL_DATASET_STAGE = "final_dataset_construction"
_SUPPORTED_PREPROCESSING_STAGE_ORDER = (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
)


class DatasetBuildExecutor:
    """Construct `AnalysisReadyPhosphoDataset` from interpreted builder input.

    Default policy uses the identity transformer, which is a pass-through
    establisher for already-prepared quantitative matrices after internal
    preprocessing stages (including optional site-matrix construction).
    """

    def __init__(
        self,
        *,
        transformer: Transformer | None = None,
        intensity_scale_resolver: DatasetIntensityScaleResolver | None = None,
        preprocessor: DatasetPreprocessorContract | None = None,
    ) -> None:
        self._intensity_scale_resolver = (
            intensity_scale_resolver
            or DatasetIntensityScaleResolver(
                transformer=transformer or IdentityTransformer()
            )
        )
        self._preprocessor = preprocessor or DatasetPreprocessor()

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        preprocessed = self._preprocessor.run(
            phospho=request.phospho,
            site_metadata=request.site_metadata,
            sample_metadata=request.sample_metadata,
            total=request.total,
            plan=request.preprocessing_plan,
        )
        resolved = self._intensity_scale_resolver.run(
            phospho=preprocessed.phospho,
            total=preprocessed.total,
            expected_scale_kind=_resolve_expected_intensity_scale_kind(
                request.preprocessing_plan
            ),
        )
        if not resolved.intensity_scale_state.is_established:
            raise TransformationStateEstablishmentError(
                "intensity-scale resolver returned a non-established "
                "intensity scale state; this violates the dataset boundary "
                "contract"
            )
        processing_state = build_dataset_processing_state(
            plan=request.preprocessing_plan,
            intensity_scale_state=resolved.intensity_scale_state,
            preprocessing_trace=preprocessed.preprocessing_trace,
            final_phospho=resolved.phospho,
            final_site_metadata=preprocessed.site_metadata,
            final_sample_metadata=preprocessed.sample_metadata,
        )
        intensity_scale_state = processing_state.intensity_scale
        quantitative_meaning = intensity_scale_state.quantity
        if quantitative_meaning is None:
            raise DatasetBuildError(
                "intensity-scale state is missing quantitative meaning"
            )
        report = _build_dataset_preprocessing_report(
            row_counts=preprocessed.preprocessing_row_counts,
            operations=preprocessed.preprocessing_operations,
            row_audit=preprocessed.row_audit,
            duplicate_site_resolution=preprocessed.duplicate_site_resolution,
            metadata_conflicts=preprocessed.metadata_conflicts,
            comparison_group_stats=preprocessed.comparison_group_stats,
            comparison_pair_stats=preprocessed.comparison_pair_stats,
            final_dataset_rows=int(len(resolved.phospho.index)),
            intensity_scale_label=intensity_scale_state.label,
            quantitative_meaning=quantitative_meaning.value,
        )
        provenance = _build_dataset_run_provenance(
            request=request,
            preprocessed=preprocessed,
            resolved_phospho=resolved.phospho,
            resolved_total=resolved.total,
            preprocessing_trace=preprocessed.preprocessing_trace,
            intensity_scale_label=intensity_scale_state.label,
            quantitative_meaning=quantitative_meaning.value,
        )
        return AnalysisReadyPhosphoDataset._from_owned(
            phospho=resolved.phospho,
            site_metadata=preprocessed.site_metadata,
            sample_metadata=preprocessed.sample_metadata,
            total=resolved.total,
            comparisons=preprocessed.comparisons,
            organism=request.organism,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            preprocessing_report=report,
            provenance=provenance,
        )


def _build_dataset_preprocessing_report(
    *,
    row_counts: pd.DataFrame | None,
    operations: pd.DataFrame | None,
    row_audit: pd.DataFrame | None,
    duplicate_site_resolution: pd.DataFrame | None,
    metadata_conflicts: pd.DataFrame | None,
    comparison_group_stats: pd.DataFrame | None,
    comparison_pair_stats: pd.DataFrame | None,
    final_dataset_rows: int,
    intensity_scale_label: str,
    quantitative_meaning: str,
) -> DatasetPreprocessingReport:
    row_count_rows = list(row_count_rows_from_dataframe(row_counts))
    operation_rows = list(operation_rows_from_dataframe(operations))
    row_audit_rows = row_audit_rows_from_dataframe(row_audit)
    duplicate_site_resolution_rows = duplicate_site_resolution_rows_from_dataframe(
        duplicate_site_resolution
    )
    metadata_conflict_rows = metadata_conflict_rows_from_dataframe(metadata_conflicts)
    comparison_group_stats_rows = comparison_group_stats_rows_from_dataframe(
        comparison_group_stats
    )
    comparison_pair_stats_rows = comparison_pair_stats_rows_from_dataframe(
        comparison_pair_stats
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
                "quantitative_meaning": quantitative_meaning,
            },
            input_rows=final_dataset_rows,
            output_rows=final_dataset_rows,
            notes="analysis-ready dataset boundary construction",
        )
    )
    return DatasetPreprocessingReport.from_rows(
        row_count_rows=tuple(row_count_rows),
        operation_rows=tuple(operation_rows),
        row_audit_rows=row_audit_rows,
        duplicate_site_resolution_rows=duplicate_site_resolution_rows,
        metadata_conflict_rows=metadata_conflict_rows,
        comparison_group_stats_rows=comparison_group_stats_rows,
        comparison_pair_stats_rows=comparison_pair_stats_rows,
    )


def _resolve_expected_intensity_scale_kind(
    preprocessing_plan: PreprocessingPlan,
) -> IntensityScaleKind:
    if preprocessing_plan.intensity_transform_policy is IntensityTransformPolicy.LOG2:
        return IntensityScaleKind.LOG2
    return IntensityScaleKind.LINEAR


def _build_dataset_run_provenance(
    *,
    request: InterpretedDatasetBuildRequest,
    preprocessed: PreprocessedDatasetBuildTables,
    resolved_phospho: pd.DataFrame,
    resolved_total: pd.DataFrame | None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    intensity_scale_label: str,
    quantitative_meaning: str,
) -> RunProvenance:
    input_tables = _collect_fingerprints(
        (
            ("dataset.phospho", request.phospho),
            ("dataset.site_metadata", request.site_metadata),
            ("dataset.sample_metadata", request.sample_metadata),
            ("dataset.total", request.total),
        )
    )
    output_tables = _collect_fingerprints(
        (
            ("dataset.phospho", resolved_phospho),
            ("dataset.site_metadata", preprocessed.site_metadata),
            ("dataset.sample_metadata", preprocessed.sample_metadata),
            ("dataset.total", resolved_total),
            ("dataset.comparisons", preprocessed.comparisons),
        )
    )
    return RunProvenance(
        environment=collect_environment_provenance(),
        input_tables=input_tables,
        preprocessing_stages=_stage_trace_to_provenance(preprocessing_trace),
        reference=None,
        workflow_name="dataset_builder",
        workflow_parameters={
            "preprocessing_plan": _preprocessing_plan_to_payload(
                request.preprocessing_plan
            ),
            "intensity_scale_label": intensity_scale_label,
            "quantitative_meaning": quantitative_meaning,
        },
        random_state=None,
        random_seed_policy=None,
        output_tables=output_tables,
        scientific_policies=_dataset_scientific_policies(request.preprocessing_plan),
    )


def _collect_fingerprints(
    entries: tuple[tuple[str, pd.DataFrame | None], ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _stage_trace_to_provenance(
    trace: tuple[PreprocessingStageExecution, ...] | None,
) -> tuple[PreprocessingStageProvenance, ...]:
    if trace is None:
        return ()
    return tuple(
        PreprocessingStageProvenance(
            stage=item.stage,
            operation=item.operation,
            parameters=dict(item.parameters),
            input_shape=item.input_shape,
            output_shape=item.output_shape,
            input_hash=item.input_hash,
            output_hash=item.output_hash,
            dropped_row_ids=item.dropped_row_ids,
            dropped_row_count=int(item.dropped_row_count),
            schema_version=int(item.schema_version),
            consumed_input_tables=tuple(item.consumed_input_tables),
            produced_output_tables=tuple(item.produced_output_tables),
            backend=item.backend,
            random_seed=item.random_seed,
            is_deterministic=bool(item.is_deterministic),
            imputed_cell_count=int(item.imputed_cell_count),
            imputed_row_ids=item.imputed_row_ids,
            notes=item.notes,
            diagnostics=_to_json_mapping(item.diagnostics),
        )
        for item in trace
    )


def _to_json_mapping(values: Mapping[str, object]) -> dict[str, JsonValue]:
    return {str(key): _to_json_value(value) for key, value in values.items()}


def _to_json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _to_json_value(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_to_json_value(item) for item in value]
    if isinstance(value, list):
        return [_to_json_value(item) for item in value]
    return str(value)


def _preprocessing_plan_to_payload(plan: PreprocessingPlan) -> dict[str, object]:
    payload: dict[str, object] = {
        "intensity_transform_policy": plan.intensity_transform_policy.value,
        "intensity_transform_pseudocount": float(plan.intensity_transform_pseudocount),
        "normalisation_policy": plan.normalisation_policy.value,
        "missing_data_policy": plan.missing_data_policy.value,
        "missing_data_min_observed_values": plan.missing_data_min_observed_values,
        "missing_data_q": plan.missing_data_q,
        "missing_data_width": plan.missing_data_width,
        "missing_data_seed": plan.missing_data_seed,
        "missing_data_k": plan.missing_data_k,
        "missing_data_distance": plan.missing_data_distance,
        "missing_data_max_missing_fraction_per_row": (
            plan.missing_data_max_missing_fraction_per_row
        ),
        "site_sequence_resolution_enabled": plan.site_sequence_resolution_enabled,
        "site_sequence_resolution_fasta_path": (
            plan.site_sequence_resolution_fasta_path
        ),
        "site_sequence_resolution_mode": plan.site_sequence_resolution_mode.value,
        "site_sequence_resolution_flank_size": int(
            plan.site_sequence_resolution_flank_size
        ),
        "site_sequence_resolution_accession_column": (
            plan.site_sequence_resolution_accession_column
        ),
        "site_sequence_resolution_site_column": (
            plan.site_sequence_resolution_site_column
        ),
        "total_protein_correction_policy": plan.total_protein_correction_policy.value,
        "total_protein_correction_identity_policy": (
            _total_correction_identity_policy_to_payload(
                plan.total_protein_correction_identity_policy
            )
        ),
        "site_matrix_policy": plan.site_matrix_policy.value,
        "comparison_building_policy": plan.comparison_building_policy.value,
        "site_matrix_duplicate_site_policy": plan.site_matrix_duplicate_site_policy.value,
        "site_matrix_missing_data_policy": plan.site_matrix_missing_data_policy.value,
        "site_matrix_minimum_observed_values": plan.site_matrix_minimum_observed_values,
        "comparison_sample_group_column": plan.comparison_sample_group_column,
        "comparison_pairs": (
            None if plan.comparison_pairs is None else list(plan.comparison_pairs)
        ),
        "ruv_readiness_enabled": bool(plan.ruv_readiness_enabled),
        "ruv_readiness_control_feature_column": (
            plan.ruv_readiness_control_feature_column
        ),
        "ruv_readiness_replicate_group_column": (
            plan.ruv_readiness_replicate_group_column
        ),
        "ruv_readiness_batch_column": plan.ruv_readiness_batch_column,
        "stage_order": list(plan.stage_order),
    }
    return payload


def _total_correction_identity_policy_to_payload(
    policy: TotalProteinCorrectionIdentityPolicy,
) -> dict[str, object]:
    return {
        "mode": str(policy.mode),
        "matching_policy": str(policy.matching_policy),
        "phosphosite_key": policy.phosphosite_key,
        "total_protein_key": policy.total_protein_key,
        "mapping_phosphosite_key": policy.mapping_phosphosite_key,
        "mapping_total_protein_key": policy.mapping_total_protein_key,
        "mapping_table_fingerprint": policy.mapping_table_fingerprint,
        "mapping_table_row_count": (
            None if policy.mapping_table is None else int(len(policy.mapping_table))
        ),
        "duplicate_policy": str(policy.duplicate_policy),
        "unmatched_policy": str(policy.unmatched_policy),
    }


def _dataset_scientific_policies(
    preprocessing_plan: PreprocessingPlan,
) -> tuple[ScientificPolicyRecord, ...]:
    policies = [
        PreprocessingStageOrderPolicy(
            configured_stage_order=tuple(
                str(stage) for stage in preprocessing_plan.stage_order
            ),
            default_stage_order=tuple(
                str(stage) for stage in DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT
            ),
            supported_stage_order=_SUPPORTED_PREPROCESSING_STAGE_ORDER,
        ).record,
    ]
    if DATASET_PREPROCESSING_STAGE_SITE_MATRIX in preprocessing_plan.stage_order:
        policies.append(
            build_duplicate_site_resolution_policy(
                duplicate_site_policy=(
                    preprocessing_plan.site_matrix_duplicate_site_policy.value
                )
            )
        )
    return tuple(policies)
