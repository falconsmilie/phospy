"""Batch-correction provenance payload serialization."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import (
    BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1,
    BatchCorrectionProvenance,
    BatchCorrectionRejectedEntity,
)
from phospy.provenance.serialization._payload import (
    optional_mapping,
    optional_str,
    require_int,
    require_mapping,
    require_sequence,
    require_str,
    to_json_safe,
    to_json_value,
)
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
    table_fingerprints_from_payload,
)


def batch_correction_provenance_to_payload(
    provenance: BatchCorrectionProvenance,
) -> dict[str, object]:
    """Serialize batch-correction provenance to a JSON-safe payload."""

    payload = {
        "schema_version": int(provenance.schema_version),
        "requested_method": provenance.requested_method,
        "resolved_parameters": to_json_safe(provenance.resolved_parameters),
        "preprocessing_stage_order": list(provenance.preprocessing_stage_order),
        "control_site_source": to_json_safe(provenance.control_site_source),
        "selected_site_key_rows": list(provenance.selected_site_key_rows),
        "batch_metadata": to_json_safe(provenance.batch_metadata),
        "replicate_metadata": (
            None
            if provenance.replicate_metadata is None
            else to_json_safe(provenance.replicate_metadata)
        ),
        "design_metadata": to_json_safe(provenance.design_metadata),
        "missing_value_policy": to_json_safe(provenance.missing_value_policy),
        "observation_masks": [
            table_fingerprint_to_payload(item) for item in provenance.observation_masks
        ],
        "input_matrix_fingerprint": table_fingerprint_to_payload(
            provenance.input_matrix_fingerprint
        ),
        "output_matrix_fingerprint": (
            None
            if provenance.output_matrix_fingerprint is None
            else table_fingerprint_to_payload(provenance.output_matrix_fingerprint)
        ),
        "diagnostics": to_json_safe(provenance.diagnostics),
        "warnings": list(provenance.warnings),
        "rejected_entities": [
            _batch_correction_rejected_entity_to_payload(item)
            for item in provenance.rejected_entities
        ],
        "phospy_version": provenance.phospy_version,
        "python_version": provenance.python_version,
        "dependency_versions": to_json_safe(provenance.dependency_versions),
    }
    if provenance.imputation_policy:
        payload["imputation_policy"] = to_json_safe(provenance.imputation_policy)
    return payload


def batch_correction_provenance_from_payload(
    payload: Mapping[str, object],
) -> BatchCorrectionProvenance:
    """Deserialize batch-correction provenance from a decoded payload."""

    payload = require_mapping(
        payload,
        field_name="batch_correction_provenance",
    )
    schema_version = require_int(
        payload.get("schema_version"),
        field_name="batch_correction_provenance.schema_version",
    )
    if schema_version != BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1:
        raise PhosPyInputError(
            f"Unsupported batch-correction provenance schema version: {schema_version}"
        )
    resolved_parameters = require_mapping(
        payload.get("resolved_parameters"),
        field_name="batch_correction_provenance.resolved_parameters",
    )
    control_site_source = require_mapping(
        payload.get("control_site_source"),
        field_name="batch_correction_provenance.control_site_source",
    )
    batch_metadata = require_mapping(
        payload.get("batch_metadata"),
        field_name="batch_correction_provenance.batch_metadata",
    )
    replicate_metadata = optional_mapping(
        payload.get("replicate_metadata"),
        field_name="batch_correction_provenance.replicate_metadata",
    )
    design_metadata = require_mapping(
        payload.get("design_metadata"),
        field_name="batch_correction_provenance.design_metadata",
    )
    missing_value_policy = require_mapping(
        payload.get("missing_value_policy"),
        field_name="batch_correction_provenance.missing_value_policy",
    )
    diagnostics = require_mapping(
        payload.get("diagnostics"),
        field_name="batch_correction_provenance.diagnostics",
    )
    dependency_versions = require_mapping(
        payload.get("dependency_versions"),
        field_name="batch_correction_provenance.dependency_versions",
    )
    imputation_policy = optional_mapping(
        payload.get("imputation_policy"),
        field_name="batch_correction_provenance.imputation_policy",
    )
    output_matrix_payload = payload.get("output_matrix_fingerprint")
    return BatchCorrectionProvenance(
        schema_version=schema_version,
        requested_method=require_str(
            payload.get("requested_method"),
            field_name="batch_correction_provenance.requested_method",
        ),
        resolved_parameters={
            key: to_json_value(value) for key, value in resolved_parameters.items()
        },
        preprocessing_stage_order=tuple(
            require_str(
                item,
                field_name="batch_correction_provenance.preprocessing_stage_order[]",
            )
            for item in require_sequence(
                payload.get("preprocessing_stage_order"),
                field_name="batch_correction_provenance.preprocessing_stage_order",
            )
        ),
        control_site_source={
            key: to_json_value(value) for key, value in control_site_source.items()
        },
        selected_site_key_rows=tuple(
            require_str(
                item,
                field_name="batch_correction_provenance.selected_site_key_rows[]",
            )
            for item in require_sequence(
                payload.get("selected_site_key_rows"),
                field_name="batch_correction_provenance.selected_site_key_rows",
            )
        ),
        batch_metadata={
            key: to_json_value(value) for key, value in batch_metadata.items()
        },
        replicate_metadata=None
        if replicate_metadata is None
        else {key: to_json_value(value) for key, value in replicate_metadata.items()},
        design_metadata={
            key: to_json_value(value) for key, value in design_metadata.items()
        },
        missing_value_policy={
            key: to_json_value(value) for key, value in missing_value_policy.items()
        },
        observation_masks=table_fingerprints_from_payload(
            payload.get("observation_masks"),
            field_name="batch_correction_provenance.observation_masks",
        ),
        input_matrix_fingerprint=table_fingerprint_from_payload(
            require_mapping(
                payload.get("input_matrix_fingerprint"),
                field_name="batch_correction_provenance.input_matrix_fingerprint",
            )
        ),
        output_matrix_fingerprint=None
        if output_matrix_payload is None
        else table_fingerprint_from_payload(
            require_mapping(
                output_matrix_payload,
                field_name="batch_correction_provenance.output_matrix_fingerprint",
            )
        ),
        diagnostics={key: to_json_value(value) for key, value in diagnostics.items()},
        warnings=tuple(
            require_str(
                item,
                field_name="batch_correction_provenance.warnings[]",
            )
            for item in require_sequence(
                payload.get("warnings"),
                field_name="batch_correction_provenance.warnings",
            )
        ),
        rejected_entities=tuple(
            _batch_correction_rejected_entity_from_payload(
                require_mapping(
                    item,
                    field_name=(
                        f"batch_correction_provenance.rejected_entities[{position}]"
                    ),
                )
            )
            for position, item in enumerate(
                require_sequence(
                    payload.get("rejected_entities"),
                    field_name="batch_correction_provenance.rejected_entities",
                )
            )
        ),
        phospy_version=require_str(
            payload.get("phospy_version"),
            field_name="batch_correction_provenance.phospy_version",
        ),
        python_version=optional_str(
            payload.get("python_version"),
            field_name="batch_correction_provenance.python_version",
        )
        or "unknown",
        dependency_versions={
            key: (
                None
                if value is None
                else require_str(
                    value,
                    field_name=(
                        f"batch_correction_provenance.dependency_versions['{key}']"
                    ),
                )
            )
            for key, value in dependency_versions.items()
        },
        imputation_policy={}
        if imputation_policy is None
        else {key: to_json_value(value) for key, value in imputation_policy.items()},
    )


def _batch_correction_rejected_entity_to_payload(
    entity: BatchCorrectionRejectedEntity,
) -> dict[str, object]:
    return {
        "entity_type": entity.entity_type,
        "identifier": entity.identifier,
        "reason": entity.reason,
        "details": None if entity.details is None else to_json_safe(entity.details),
    }


def _batch_correction_rejected_entity_from_payload(
    payload: Mapping[str, object],
) -> BatchCorrectionRejectedEntity:
    details = optional_mapping(
        payload.get("details"),
        field_name="batch_correction_rejected_entity.details",
    )
    return BatchCorrectionRejectedEntity(
        entity_type=require_str(
            payload.get("entity_type"),
            field_name="batch_correction_rejected_entity.entity_type",
        ),
        identifier=require_str(
            payload.get("identifier"),
            field_name="batch_correction_rejected_entity.identifier",
        ),
        reason=require_str(
            payload.get("reason"),
            field_name="batch_correction_rejected_entity.reason",
        ),
        details=None
        if details is None
        else {key: to_json_value(value) for key, value in details.items()},
    )
