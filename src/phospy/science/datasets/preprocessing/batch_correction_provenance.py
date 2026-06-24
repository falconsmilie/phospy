"""Native dataset batch-correction provenance assembly."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import cast

import pandas as pd

from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_matrix
from phospy.provenance.models import (
    BatchCorrectionProvenance,
    BatchCorrectionRejectedEntity,
    JsonValue,
)
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.batch_correction_metadata import (
    ResolvedBatchCorrectionMetadata,
    levels_in_sample_order,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan


def build_native_batch_correction_provenance(
    *,
    input_matrix: pd.DataFrame,
    output_matrix: pd.DataFrame | None,
    plan: PreprocessingPlan,
    report: BatchCorrectionReport,
    metadata: ResolvedBatchCorrectionMetadata | None,
    diagnostics: Mapping[str, object],
    warnings: Sequence[str] = (),
    observation_mask: pd.DataFrame | None = None,
    corrected_cell_status: pd.DataFrame | None = None,
    control_site_source: Mapping[str, object] | None = None,
    selected_site_key_rows: Sequence[str] = (),
    rejected_entities: Sequence[BatchCorrectionRejectedEntity] = (),
    source: str,
) -> BatchCorrectionProvenance:
    """Build provenance for native dataset correction without altering results."""

    environment = collect_environment_provenance()
    resolved_observation_mask = _resolve_observation_mask(
        input_matrix=input_matrix,
        output_matrix=output_matrix,
        observation_mask=observation_mask,
    )
    observation_masks = (
        fingerprint_matrix(
            resolved_observation_mask.astype("int8"),
            name="batch_correction.native.observation_mask",
        ),
    )
    if corrected_cell_status is not None:
        observation_masks = (
            *observation_masks,
            fingerprint_matrix(
                corrected_cell_status,
                name="batch_correction.native.corrected_cell_status",
            ),
        )

    missing_value_policy = _missing_value_policy_payload(plan)
    imputation_policy = _imputation_policy_payload(plan)
    return BatchCorrectionProvenance(
        requested_method=str(report.method),
        resolved_parameters=_json_mapping(
            {
                "source": source,
                "report": report.to_payload(),
                "plan": _plan_payload(plan),
                "diagnostics": dict(diagnostics),
            }
        ),
        preprocessing_stage_order=tuple(str(stage) for stage in plan.stage_order),
        control_site_source=_json_mapping(
            control_site_source
            if control_site_source is not None
            else {
                "source_type": "not_applicable",
                "reason": "linear_residualize_batch does not use control sites",
            }
        ),
        selected_site_key_rows=tuple(str(row) for row in selected_site_key_rows),
        batch_metadata=_json_mapping(_batch_metadata_payload(plan, metadata, report)),
        replicate_metadata=None,
        design_metadata=_json_mapping(_design_metadata_payload(plan, metadata, report)),
        missing_value_policy=_json_mapping(missing_value_policy),
        observation_masks=observation_masks,
        input_matrix_fingerprint=fingerprint_matrix(
            input_matrix,
            name="batch_correction.native.input",
        ),
        output_matrix_fingerprint=(
            None
            if output_matrix is None
            else fingerprint_matrix(
                output_matrix,
                name="batch_correction.native.corrected",
            )
        ),
        diagnostics=_json_mapping(
            {
                "status": report.status,
                "report": report.to_payload(),
                "stage_diagnostics": dict(diagnostics),
            }
        ),
        warnings=tuple(str(warning) for warning in warnings),
        rejected_entities=tuple(rejected_entities),
        phospy_version=environment.package_version,
        dependency_versions=environment.dependency_versions,
        imputation_policy=_json_mapping(imputation_policy),
    )


def _resolve_observation_mask(
    *,
    input_matrix: pd.DataFrame,
    output_matrix: pd.DataFrame | None,
    observation_mask: pd.DataFrame | None,
) -> pd.DataFrame:
    if observation_mask is not None:
        return observation_mask.copy(deep=True)
    reference = input_matrix if output_matrix is None else output_matrix
    return pd.DataFrame(
        input_matrix.notna().reindex(
            index=reference.index,
            columns=reference.columns,
            fill_value=False,
        ),
        index=reference.index.copy(),
        columns=reference.columns.copy(),
    )


def _plan_payload(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "batch_correction_method": plan.batch_correction_method,
        "batch_correction_batch_column": plan.batch_correction_batch_column,
        "batch_correction_condition_column": plan.batch_correction_condition_column,
        "batch_correction_preserve_condition_effects": (
            plan.batch_correction_preserve_condition_effects
        ),
        "missing_data_policy": plan.missing_data_policy.value,
        "missing_data_seed": plan.missing_data_seed,
        "stage_order": list(plan.stage_order),
        "resolved_stage_order": [
            {
                "stage": item.stage,
                "order_index": int(item.order_index),
                "rationale": item.rationale,
            }
            for item in plan.stage_order_resolution
        ],
    }


def _batch_metadata_payload(
    plan: PreprocessingPlan,
    metadata: ResolvedBatchCorrectionMetadata | None,
    report: BatchCorrectionReport,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "column": plan.batch_correction_batch_column,
        "levels": list(report.batch_levels),
    }
    if metadata is None:
        payload["available"] = False
        return payload
    payload.update(
        {
            "available": True,
            "sample_order": list(metadata.sample_order),
            "batch_by_sample": dict(metadata.batch_by_sample),
            "batch_labels": list(metadata.batch_labels),
        }
    )
    return payload


def _design_metadata_payload(
    plan: PreprocessingPlan,
    metadata: ResolvedBatchCorrectionMetadata | None,
    report: BatchCorrectionReport,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "condition_columns": [plan.batch_correction_condition_column],
        "condition_levels": list(report.condition_levels),
        "preserve_condition_effects": bool(
            plan.batch_correction_preserve_condition_effects
        ),
        "design_preservation_policy": report.design_preservation_policy,
        "confounding_check_status": report.confounding_check_status,
    }
    if metadata is None:
        payload["available"] = False
        return payload
    payload.update(
        {
            "available": True,
            "sample_order": list(metadata.sample_order),
            "condition_by_sample": dict(metadata.condition_by_sample),
            "condition_labels": list(metadata.condition_labels),
            "batch_levels_in_sample_order": list(
                levels_in_sample_order(
                    metadata.batch_by_sample,
                    sample_order=metadata.sample_order,
                )
            ),
            "condition_levels_in_sample_order": list(
                levels_in_sample_order(
                    metadata.condition_by_sample,
                    sample_order=metadata.sample_order,
                )
            ),
        }
    )
    return payload


def _missing_value_policy_payload(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "policy": "reject_missing_at_batch_correction",
        "upstream_missing_data_policy": plan.missing_data_policy.value,
        "requires_complete_matrix": True,
        "imputation_policy": _imputation_policy_payload(plan),
    }


def _imputation_policy_payload(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "policy": plan.missing_data_policy.value,
        "seed": plan.missing_data_seed,
        "min_observed_values": plan.missing_data_min_observed_values,
        "q": plan.missing_data_q,
        "width": plan.missing_data_width,
        "k": plan.missing_data_k,
        "distance": plan.missing_data_distance,
        "max_missing_fraction_per_row": (
            plan.missing_data_max_missing_fraction_per_row
        ),
    }


def _json_mapping(value: Mapping[str, object]) -> Mapping[str, JsonValue]:
    return cast(Mapping[str, JsonValue], _json_value(value))


def _json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, str | int | float | bool):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
        return [_json_value(item) for item in value]
    return str(value)


__all__ = ["build_native_batch_correction_provenance"]
