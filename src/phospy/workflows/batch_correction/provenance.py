"""Provenance assembly for the batch-correction workflow."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Protocol, cast, runtime_checkable

import pandas as pd

from phospy.provenance import fingerprint_matrix
from phospy.provenance.models import (
    BatchCorrectionProvenance,
    BatchCorrectionRejectedEntity,
    JsonValue,
    TableFingerprint,
)
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteEligibility,
    ControlSiteMapping,
)
from phospy.validation.datasets.batch_correction import ResolvedBatchDesignMetadata
from phospy.workflows.batch_correction.contracts import (
    BatchCorrectionExecutorResultContract,
    BatchCorrectionWorkflowRequest,
)
from phospy.workflows.batch_correction.interpreter import ResolvedBatchCorrectionPlan


@runtime_checkable
class _PayloadProvider(Protocol):
    def to_payload(self) -> Mapping[str, object]: ...


class BatchCorrectionProvenanceRecorder:
    """Assemble typed workflow provenance from validated collaborators."""

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
        dataset_metadata: ResolvedBatchDesignMetadata,
        control_site_mapping: ControlSiteMapping,
        missingness_policy: object,
        plan: ResolvedBatchCorrectionPlan,
        executor_result: BatchCorrectionExecutorResultContract,
        **_: object,
    ) -> BatchCorrectionProvenance:
        output_mask = executor_result.output_observation_mask.astype("int8")
        observation_masks = (
            fingerprint_matrix(
                output_mask,
                name="batch_correction.workflow.output_observation_mask",
            ),
        )
        if request.upstream_observation_mask is not None:
            observation_masks = (
                fingerprint_matrix(
                    request.upstream_observation_mask.astype("int8"),
                    name="batch_correction.workflow.upstream_observation_mask",
                ),
                *observation_masks,
            )
        return BatchCorrectionProvenance(
            requested_method=request.config.method.value,
            resolved_parameters=_json_mapping(
                {
                    "config": _config_payload(request),
                    "interpreter_plan": plan.to_payload(),
                    "interpreter_seed_data": dict(plan.provenance_seed_data),
                    "executor": dict(executor_result.provenance_payload),
                }
            ),
            preprocessing_stage_order=tuple(plan.stage_order),
            control_site_source=_json_mapping(
                {
                    "source": request.config.control_site_source.value,
                    "mode": request.config.control_site_mode.value,
                    "selected_control_count": len(plan.eligible_control_site_rows),
                }
            ),
            selected_site_key_rows=tuple(
                row.site_key for row in plan.eligible_control_site_rows
            ),
            batch_metadata=_json_mapping(
                {
                    "sample_order": list(dataset_metadata.sample_order),
                    "batch_by_sample": dict(dataset_metadata.batch_by_sample),
                    "condition_by_sample": dict(dataset_metadata.condition_by_sample),
                }
            ),
            replicate_metadata=_replicate_metadata(dataset_metadata),
            design_metadata=_json_mapping(
                {
                    "resolved_design_matrix": _frame_payload(
                        plan.resolved_design_matrix
                    ),
                    "batch_terms": list(plan.batch_terms),
                    "condition_terms_to_preserve": (
                        list(plan.condition_terms_to_preserve)
                    ),
                }
            ),
            missing_value_policy=_json_mapping(_payload(missingness_policy)),
            observation_masks=observation_masks,
            input_matrix_fingerprint=fingerprint_matrix(
                request.phospho,
                name="batch_correction.workflow.input",
            ),
            output_matrix_fingerprint=fingerprint_matrix(
                executor_result.corrected_matrix,
                name="batch_correction.workflow.corrected",
            ),
            diagnostics=_json_mapping(
                {
                    "interpreter_diagnostic_requirements": (
                        plan.diagnostic_requirements.to_payload()
                    ),
                    "executor": executor_result.diagnostics.to_payload(),
                }
            ),
            warnings=tuple(str(warning) for warning in executor_result.warnings),
            rejected_entities=(
                *_control_site_rejected_entities(control_site_mapping),
                *_rejected_entities(executor_result),
            ),
        )


def _config_payload(request: BatchCorrectionWorkflowRequest) -> dict[str, object]:
    config = request.config
    return {
        "method": config.method.value,
        "batch_column": config.batch_column,
        "condition_columns": list(config.condition_columns),
        "replicate_column": config.replicate_column,
        "control_site_source": config.control_site_source.value,
        "control_site_mode": config.control_site_mode.value,
        "missing_value_policy": config.missing_value_policy.value,
        "imputation_policy": config.imputation_policy.value,
        "n_unwanted_factors": config.n_unwanted_factors,
        "requested_stage_order": config.stage_order.value,
        "stage_order": config.stage_order.value,
        "diagnostics_enabled": config.diagnostics_enabled,
    }


def _replicate_metadata(
    dataset_metadata: ResolvedBatchDesignMetadata,
) -> Mapping[str, JsonValue] | None:
    if dataset_metadata.replicate_by_sample is None:
        return None
    return _json_mapping(
        {
            "replicate_by_sample": dict(dataset_metadata.replicate_by_sample),
            "replicate_labels": list(dataset_metadata.replicate_labels or ()),
        }
    )


def _rejected_entities(
    executor_result: BatchCorrectionExecutorResultContract,
) -> tuple[BatchCorrectionRejectedEntity, ...]:
    entities: list[BatchCorrectionRejectedEntity] = []
    for row in executor_result.rejected_rows:
        entities.append(
            BatchCorrectionRejectedEntity(
                entity_type="row",
                identifier=str(row),
                reason="executor_rejected_row",
            )
        )
    for feature_id, sample_id in executor_result.rejected_cells:
        entities.append(
            BatchCorrectionRejectedEntity(
                entity_type="cell",
                identifier=f"{feature_id}|{sample_id}",
                reason="executor_rejected_cell",
                details={"feature_id": feature_id, "sample_id": sample_id},
            )
        )
    for row in executor_result.withheld_rows:
        entities.append(
            BatchCorrectionRejectedEntity(
                entity_type="row",
                identifier=str(row),
                reason="executor_withheld_row",
            )
        )
    for feature_id, sample_id in executor_result.withheld_cells:
        entities.append(
            BatchCorrectionRejectedEntity(
                entity_type="cell",
                identifier=f"{feature_id}|{sample_id}",
                reason="executor_withheld_cell",
                details={"feature_id": feature_id, "sample_id": sample_id},
            )
        )
    return tuple(entities)


def _control_site_rejected_entities(
    control_site_mapping: ControlSiteMapping,
) -> tuple[BatchCorrectionRejectedEntity, ...]:
    entities: list[BatchCorrectionRejectedEntity] = []
    for scope, rows in (
        ("row_eligibility", control_site_mapping.row_eligibility),
        ("unmapped_annotations", control_site_mapping.unmapped_annotations),
    ):
        for row in rows:
            if row.is_control or not _is_reportable_control_rejection(row):
                continue
            identifier = row.site_key or "<missing site_key>"
            entities.append(
                BatchCorrectionRejectedEntity(
                    entity_type="site",
                    identifier=str(identifier),
                    reason=_control_rejection_reason(row),
                    details={
                        "scope": scope,
                        "control_status": row.control_status.value,
                        "valid": row.valid,
                        "reasons": list(row.reasons),
                        "row_position": row.row_position,
                        "annotation_indices": list(row.annotation_indices),
                        "exclusion_reason": row.exclusion_reason,
                    },
                )
            )
    return tuple(entities)


def _is_reportable_control_rejection(row: ControlSiteEligibility) -> bool:
    annotation_count = row.annotation_count
    return bool(
        getattr(row, "reasons", ())
        or row.control_status.value in {"excluded", "invalid", "unknown"}
        or annotation_count > 0
    )


def _control_rejection_reason(row: ControlSiteEligibility) -> str:
    reasons = tuple(str(reason) for reason in row.reasons)
    if reasons:
        return reasons[0]
    status = row.control_status.value
    if status == "excluded":
        exclusion_reason = row.exclusion_reason
        if exclusion_reason:
            return str(exclusion_reason)
        return "excluded_control_site"
    if status == "non_control":
        return "not_marked_as_control"
    if status == "unknown":
        return "unknown_control_status"
    if status == "invalid":
        return "invalid_control_site_annotation"
    return "not_eligible_control_site"


def _frame_payload(frame: pd.DataFrame) -> dict[str, object]:
    return {
        "index": [str(value) for value in frame.index.tolist()],
        "columns": [str(value) for value in frame.columns.tolist()],
        "data": [
            [float(value) for value in row]
            for row in frame.to_numpy(dtype="float64").tolist()
        ],
    }


def _payload(value: object) -> Mapping[str, object]:
    if isinstance(value, _PayloadProvider):
        payload = value.to_payload()
        if isinstance(payload, Mapping):
            return payload
    if isinstance(value, Mapping):
        return value
    return {"value": repr(value)}


def _json_mapping(value: Mapping[str, object]) -> Mapping[str, JsonValue]:
    return cast(Mapping[str, JsonValue], _json_value(value))


def _json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, str | int | float | bool):
        return value
    if isinstance(value, TableFingerprint):
        return _json_value(
            {
                "name": value.name,
                "rows": value.rows,
                "columns": value.columns,
                "index_name": value.index_name,
                "column_names": list(value.column_names),
                "dtypes": list(value.dtypes),
                "exact_hash_algorithm": value.exact_hash_algorithm,
                "exact_hash_value": value.exact_hash_value,
                "tolerance_hash_algorithm": value.tolerance_hash_algorithm,
                "tolerance_hash_value": value.tolerance_hash_value,
                "index_structure": value.index_structure,
                "column_index_structure": value.column_index_structure,
            }
        )
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
        return [_json_value(item) for item in value]
    if isinstance(value, _PayloadProvider):
        return _json_value(value.to_payload())
    return str(value)


__all__ = ["BatchCorrectionProvenanceRecorder"]
