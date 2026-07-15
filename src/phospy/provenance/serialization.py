"""Serialization helpers for provenance payloads."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import (
    BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1,
    ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
    PREPROCESSING_STAGE_DETERMINISM_EXTERNALLY_NONDETERMINISTIC,
    PREPROCESSING_STAGE_DETERMINISM_PURE,
    PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    BatchCorrectionProvenance,
    BatchCorrectionRejectedEntity,
    DeterminismKind,
    EnvironmentProvenance,
    JsonValue,
    PreprocessingStageProvenance,
    ReferenceProvenance,
    ReproducibilityCaveat,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord
from phospy.science.references.models import ReferenceContext

if TYPE_CHECKING:
    from phospy.science.references.identifiers import (
        ReferenceIdentifierNormalisationRecord,
        ReferenceIdentifierNormalisationReport,
    )


_LEGACY_PROVENANCE_SCHEMA_ERROR = (
    "Legacy provenance schemas are no longer supported. Regenerate the result "
    "with the current PhosPy version."
)
_LEGACY_TABLE_FINGERPRINT_FIELDS = frozenset({"hash_algorithm", "hash_value"})
_LEGACY_PREPROCESSING_STAGE_FIELDS = frozenset({"is_deterministic"})
_DETERMINISM_ALIASES = {
    "pure": DeterminismKind.DETERMINISTIC,
    "external_dependency": DeterminismKind.EXTERNALLY_NONDETERMINISTIC,
}


def to_payload(provenance: RunProvenance) -> dict[str, object]:
    """Serialize run provenance to a JSON-safe payload."""

    return {
        "environment": _environment_to_payload(provenance.environment),
        "input_tables": [
            _table_fingerprint_to_payload(item) for item in provenance.input_tables
        ],
        "preprocessing_stages": [
            _stage_to_payload(item) for item in provenance.preprocessing_stages
        ],
        "reference": (
            None
            if provenance.reference is None
            else _reference_to_payload(provenance.reference)
        ),
        "reference_context": (
            None
            if provenance.reference_context is None
            else provenance.reference_context.to_payload()
        ),
        "workflow_name": provenance.workflow_name,
        "workflow_parameters": _to_json_safe(provenance.workflow_parameters),
        "random_state": provenance.random_state,
        "random_seed_policy": provenance.random_seed_policy,
        "output_tables": [
            _table_fingerprint_to_payload(item) for item in provenance.output_tables
        ],
        "scientific_policies": [
            item.to_payload() for item in provenance.scientific_policies
        ],
    }


def from_payload(payload: Mapping[str, object]) -> RunProvenance:
    """Deserialize run provenance from a decoded payload."""

    environment_payload = _require_mapping(
        payload.get("environment"),
        field_name="provenance.environment",
    )
    input_tables_payload = _require_sequence(
        payload.get("input_tables"),
        field_name="provenance.input_tables",
    )
    stages_payload = _require_sequence(
        payload.get("preprocessing_stages"),
        field_name="provenance.preprocessing_stages",
    )
    output_tables_payload = _require_sequence(
        payload.get("output_tables"),
        field_name="provenance.output_tables",
    )
    scientific_policies_payload = _require_sequence(
        payload.get("scientific_policies", []),
        field_name="provenance.scientific_policies",
    )
    reference_raw = payload.get("reference")
    if reference_raw is None:
        reference = None
    else:
        reference = _reference_from_payload(
            _require_mapping(reference_raw, field_name="provenance.reference")
        )
    reference_context_payload = _optional_mapping(
        payload.get("reference_context"),
        field_name="provenance.reference_context",
    )
    workflow_parameters = _require_mapping(
        payload.get("workflow_parameters"),
        field_name="provenance.workflow_parameters",
    )
    return RunProvenance(
        environment=_environment_from_payload(environment_payload),
        input_tables=tuple(
            _table_fingerprint_from_payload(
                _require_mapping(
                    item,
                    field_name=f"provenance.input_tables[{position}]",
                )
            )
            for position, item in enumerate(input_tables_payload)
        ),
        preprocessing_stages=tuple(
            _stage_from_payload(
                _require_mapping(
                    item,
                    field_name=f"provenance.preprocessing_stages[{position}]",
                )
            )
            for position, item in enumerate(stages_payload)
        ),
        reference=reference,
        workflow_name=_optional_str(
            payload.get("workflow_name"),
            field_name="provenance.workflow_name",
        ),
        workflow_parameters={
            str(key): _to_json_value(value)
            for key, value in workflow_parameters.items()
        },
        random_state=_optional_int(
            payload.get("random_state"),
            field_name="provenance.random_state",
        ),
        random_seed_policy=_optional_str(
            payload.get("random_seed_policy"),
            field_name="provenance.random_seed_policy",
        ),
        output_tables=tuple(
            _table_fingerprint_from_payload(
                _require_mapping(
                    item,
                    field_name=f"provenance.output_tables[{position}]",
                )
            )
            for position, item in enumerate(output_tables_payload)
        ),
        scientific_policies=tuple(
            ScientificPolicyRecord.from_payload(
                {
                    str(key): value
                    for key, value in _require_mapping(
                        item,
                        field_name=f"provenance.scientific_policies[{position}]",
                    ).items()
                }
            )
            for position, item in enumerate(scientific_policies_payload)
        ),
        reference_context=None
        if reference_context_payload is None
        else ReferenceContext.from_payload(reference_context_payload),
    )


def batch_correction_provenance_to_payload(
    provenance: BatchCorrectionProvenance,
) -> dict[str, object]:
    """Serialize batch-correction provenance to a JSON-safe payload."""

    payload = {
        "schema_version": int(provenance.schema_version),
        "requested_method": provenance.requested_method,
        "resolved_parameters": _to_json_safe(provenance.resolved_parameters),
        "preprocessing_stage_order": list(provenance.preprocessing_stage_order),
        "control_site_source": _to_json_safe(provenance.control_site_source),
        "selected_site_key_rows": list(provenance.selected_site_key_rows),
        "batch_metadata": _to_json_safe(provenance.batch_metadata),
        "replicate_metadata": (
            None
            if provenance.replicate_metadata is None
            else _to_json_safe(provenance.replicate_metadata)
        ),
        "design_metadata": _to_json_safe(provenance.design_metadata),
        "missing_value_policy": _to_json_safe(provenance.missing_value_policy),
        "observation_masks": [
            _table_fingerprint_to_payload(item) for item in provenance.observation_masks
        ],
        "input_matrix_fingerprint": _table_fingerprint_to_payload(
            provenance.input_matrix_fingerprint
        ),
        "output_matrix_fingerprint": (
            None
            if provenance.output_matrix_fingerprint is None
            else _table_fingerprint_to_payload(provenance.output_matrix_fingerprint)
        ),
        "diagnostics": _to_json_safe(provenance.diagnostics),
        "warnings": list(provenance.warnings),
        "rejected_entities": [
            _batch_correction_rejected_entity_to_payload(item)
            for item in provenance.rejected_entities
        ],
        "phospy_version": provenance.phospy_version,
        "python_version": provenance.python_version,
        "dependency_versions": _to_json_safe(provenance.dependency_versions),
    }
    if provenance.imputation_policy:
        payload["imputation_policy"] = _to_json_safe(provenance.imputation_policy)
    return payload


def batch_correction_provenance_from_payload(
    payload: Mapping[str, object],
) -> BatchCorrectionProvenance:
    """Deserialize batch-correction provenance from a decoded payload."""

    schema_version = _require_int(
        payload.get("schema_version"),
        field_name="batch_correction_provenance.schema_version",
    )
    if schema_version != BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1:
        raise PhosPyInputError(
            f"Unsupported batch-correction provenance schema version: {schema_version}"
        )
    resolved_parameters = _require_mapping(
        payload.get("resolved_parameters"),
        field_name="batch_correction_provenance.resolved_parameters",
    )
    control_site_source = _require_mapping(
        payload.get("control_site_source"),
        field_name="batch_correction_provenance.control_site_source",
    )
    batch_metadata = _require_mapping(
        payload.get("batch_metadata"),
        field_name="batch_correction_provenance.batch_metadata",
    )
    replicate_metadata = _optional_mapping(
        payload.get("replicate_metadata"),
        field_name="batch_correction_provenance.replicate_metadata",
    )
    design_metadata = _require_mapping(
        payload.get("design_metadata"),
        field_name="batch_correction_provenance.design_metadata",
    )
    missing_value_policy = _require_mapping(
        payload.get("missing_value_policy"),
        field_name="batch_correction_provenance.missing_value_policy",
    )
    diagnostics = _require_mapping(
        payload.get("diagnostics"),
        field_name="batch_correction_provenance.diagnostics",
    )
    dependency_versions = _require_mapping(
        payload.get("dependency_versions"),
        field_name="batch_correction_provenance.dependency_versions",
    )
    imputation_policy = _optional_mapping(
        payload.get("imputation_policy"),
        field_name="batch_correction_provenance.imputation_policy",
    )
    output_matrix_payload = payload.get("output_matrix_fingerprint")
    return BatchCorrectionProvenance(
        schema_version=schema_version,
        requested_method=_require_str(
            payload.get("requested_method"),
            field_name="batch_correction_provenance.requested_method",
        ),
        resolved_parameters={
            str(key): _to_json_value(value)
            for key, value in resolved_parameters.items()
        },
        preprocessing_stage_order=tuple(
            _require_str(
                item,
                field_name="batch_correction_provenance.preprocessing_stage_order[]",
            )
            for item in _require_sequence(
                payload.get("preprocessing_stage_order"),
                field_name="batch_correction_provenance.preprocessing_stage_order",
            )
        ),
        control_site_source={
            str(key): _to_json_value(value)
            for key, value in control_site_source.items()
        },
        selected_site_key_rows=tuple(
            _require_str(
                item,
                field_name="batch_correction_provenance.selected_site_key_rows[]",
            )
            for item in _require_sequence(
                payload.get("selected_site_key_rows"),
                field_name="batch_correction_provenance.selected_site_key_rows",
            )
        ),
        batch_metadata={
            str(key): _to_json_value(value) for key, value in batch_metadata.items()
        },
        replicate_metadata=None
        if replicate_metadata is None
        else {
            str(key): _to_json_value(value) for key, value in replicate_metadata.items()
        },
        design_metadata={
            str(key): _to_json_value(value) for key, value in design_metadata.items()
        },
        missing_value_policy={
            str(key): _to_json_value(value)
            for key, value in missing_value_policy.items()
        },
        observation_masks=_table_fingerprints_from_payload(
            payload.get("observation_masks"),
            field_name="batch_correction_provenance.observation_masks",
        ),
        input_matrix_fingerprint=_table_fingerprint_from_payload(
            _require_mapping(
                payload.get("input_matrix_fingerprint"),
                field_name="batch_correction_provenance.input_matrix_fingerprint",
            )
        ),
        output_matrix_fingerprint=None
        if output_matrix_payload is None
        else _table_fingerprint_from_payload(
            _require_mapping(
                output_matrix_payload,
                field_name="batch_correction_provenance.output_matrix_fingerprint",
            )
        ),
        diagnostics={
            str(key): _to_json_value(value) for key, value in diagnostics.items()
        },
        warnings=tuple(
            _require_str(
                item,
                field_name="batch_correction_provenance.warnings[]",
            )
            for item in _require_sequence(
                payload.get("warnings"),
                field_name="batch_correction_provenance.warnings",
            )
        ),
        rejected_entities=tuple(
            _batch_correction_rejected_entity_from_payload(
                _require_mapping(
                    item,
                    field_name=(
                        f"batch_correction_provenance.rejected_entities[{position}]"
                    ),
                )
            )
            for position, item in enumerate(
                _require_sequence(
                    payload.get("rejected_entities"),
                    field_name="batch_correction_provenance.rejected_entities",
                )
            )
        ),
        phospy_version=_require_str(
            payload.get("phospy_version"),
            field_name="batch_correction_provenance.phospy_version",
        ),
        python_version=_optional_str(
            payload.get("python_version"),
            field_name="batch_correction_provenance.python_version",
        )
        or "unknown",
        dependency_versions={
            str(key): (
                None
                if value is None
                else _require_str(
                    value,
                    field_name=(
                        f"batch_correction_provenance.dependency_versions['{str(key)}']"
                    ),
                )
            )
            for key, value in dependency_versions.items()
        },
        imputation_policy={}
        if imputation_policy is None
        else {
            str(key): _to_json_value(value) for key, value in imputation_policy.items()
        },
    )


def _batch_correction_rejected_entity_to_payload(
    entity: BatchCorrectionRejectedEntity,
) -> dict[str, object]:
    return {
        "entity_type": entity.entity_type,
        "identifier": entity.identifier,
        "reason": entity.reason,
        "details": None if entity.details is None else _to_json_safe(entity.details),
    }


def _batch_correction_rejected_entity_from_payload(
    payload: Mapping[str, object],
) -> BatchCorrectionRejectedEntity:
    details = _optional_mapping(
        payload.get("details"),
        field_name="batch_correction_rejected_entity.details",
    )
    return BatchCorrectionRejectedEntity(
        entity_type=_require_str(
            payload.get("entity_type"),
            field_name="batch_correction_rejected_entity.entity_type",
        ),
        identifier=_require_str(
            payload.get("identifier"),
            field_name="batch_correction_rejected_entity.identifier",
        ),
        reason=_require_str(
            payload.get("reason"),
            field_name="batch_correction_rejected_entity.reason",
        ),
        details=None
        if details is None
        else {str(key): _to_json_value(value) for key, value in details.items()},
    )


def _table_fingerprint_to_payload(fingerprint: TableFingerprint) -> dict[str, object]:
    return {
        "name": fingerprint.name,
        "rows": int(fingerprint.rows),
        "columns": int(fingerprint.columns),
        "index_name": fingerprint.index_name,
        "column_names": list(fingerprint.column_names),
        "dtypes": list(fingerprint.dtypes),
        "exact_hash_algorithm": fingerprint.exact_hash_algorithm,
        "exact_hash_value": fingerprint.exact_hash_value,
        "tolerance_hash_algorithm": fingerprint.tolerance_hash_algorithm,
        "tolerance_hash_value": fingerprint.tolerance_hash_value,
        "index_structure": (
            None
            if fingerprint.index_structure is None
            else _to_json_safe(fingerprint.index_structure)
        ),
        "column_index_structure": (
            None
            if fingerprint.column_index_structure is None
            else _to_json_safe(fingerprint.column_index_structure)
        ),
    }


def table_fingerprint_to_payload(
    fingerprint: TableFingerprint,
) -> dict[str, object]:
    """Serialize a table fingerprint to a JSON-safe payload."""

    return _table_fingerprint_to_payload(fingerprint)


def _table_fingerprint_from_payload(payload: Mapping[str, object]) -> TableFingerprint:
    _reject_legacy_provenance_fields(
        payload,
        field_name="table_fingerprint",
        legacy_fields=_LEGACY_TABLE_FINGERPRINT_FIELDS,
    )
    index_structure = _optional_mapping(
        payload.get("index_structure"),
        field_name="table_fingerprint.index_structure",
    )
    column_index_structure = _optional_mapping(
        payload.get("column_index_structure"),
        field_name="table_fingerprint.column_index_structure",
    )
    exact_hash_algorithm = _require_str(
        payload.get("exact_hash_algorithm"),
        field_name="table_fingerprint.exact_hash_algorithm",
    )
    exact_hash_value = _require_str(
        payload.get("exact_hash_value"),
        field_name="table_fingerprint.exact_hash_value",
    )
    tolerance_hash_algorithm = _require_str(
        payload.get("tolerance_hash_algorithm"),
        field_name="table_fingerprint.tolerance_hash_algorithm",
    )
    tolerance_hash_value = _require_str(
        payload.get("tolerance_hash_value"),
        field_name="table_fingerprint.tolerance_hash_value",
    )
    return TableFingerprint(
        name=_require_str(payload.get("name"), field_name="table_fingerprint.name"),
        rows=_require_int(payload.get("rows"), field_name="table_fingerprint.rows"),
        columns=_require_int(
            payload.get("columns"),
            field_name="table_fingerprint.columns",
        ),
        index_name=_optional_str(
            payload.get("index_name"),
            field_name="table_fingerprint.index_name",
        ),
        column_names=tuple(
            _require_str(item, field_name="table_fingerprint.column_names[]")
            for item in _require_sequence(
                payload.get("column_names"),
                field_name="table_fingerprint.column_names",
            )
        ),
        dtypes=tuple(
            _require_str(item, field_name="table_fingerprint.dtypes[]")
            for item in _require_sequence(
                payload.get("dtypes"),
                field_name="table_fingerprint.dtypes",
            )
        ),
        exact_hash_algorithm=exact_hash_algorithm,
        exact_hash_value=exact_hash_value,
        tolerance_hash_algorithm=tolerance_hash_algorithm,
        tolerance_hash_value=tolerance_hash_value,
        index_structure=None
        if index_structure is None
        else {
            str(key): _to_json_value(value) for key, value in index_structure.items()
        },
        column_index_structure=None
        if column_index_structure is None
        else {
            str(key): _to_json_value(value)
            for key, value in column_index_structure.items()
        },
    )


def table_fingerprint_from_payload(
    payload: Mapping[str, object],
) -> TableFingerprint:
    """Deserialize a table fingerprint from a decoded payload."""

    return _table_fingerprint_from_payload(payload)


def _environment_to_payload(environment: EnvironmentProvenance) -> dict[str, object]:
    return {
        "schema_version": int(environment.schema_version),
        "package_name": environment.package_name,
        "package_version": environment.package_version,
        "python_version": environment.python_version,
        "dependency_versions": _to_json_safe(environment.dependency_versions),
        "platform": _to_json_safe(environment.platform),
        "blas_lapack": _to_json_safe(environment.blas_lapack),
        "thread_environment": _to_json_safe(environment.thread_environment),
        "timezone": environment.timezone,
        "locale": _to_json_safe(environment.locale),
        "constraints_fingerprint": _to_json_safe(environment.constraints_fingerprint),
    }


def _environment_from_payload(payload: Mapping[str, object]) -> EnvironmentProvenance:
    dependency_versions = _require_mapping(
        payload.get("dependency_versions"),
        field_name="provenance.environment.dependency_versions",
    )
    platform_payload = _require_mapping(
        payload.get("platform", {}),
        field_name="provenance.environment.platform",
    )
    blas_lapack_payload = _require_mapping(
        payload.get("blas_lapack", {}),
        field_name="provenance.environment.blas_lapack",
    )
    thread_environment_payload = _require_mapping(
        payload.get("thread_environment", {}),
        field_name="provenance.environment.thread_environment",
    )
    locale_payload = _require_mapping(
        payload.get("locale", {}),
        field_name="provenance.environment.locale",
    )
    constraints_fingerprint_payload = _require_mapping(
        payload.get("constraints_fingerprint", {}),
        field_name="provenance.environment.constraints_fingerprint",
    )
    return EnvironmentProvenance(
        schema_version=_require_current_environment_schema_version(payload),
        package_name=_require_str(
            payload.get("package_name"),
            field_name="provenance.environment.package_name",
        ),
        package_version=_require_str(
            payload.get("package_version"),
            field_name="provenance.environment.package_version",
        ),
        python_version=_require_str(
            payload.get("python_version"),
            field_name="provenance.environment.python_version",
        ),
        dependency_versions={
            str(key): (
                None
                if value is None
                else _require_str(
                    value,
                    field_name=(
                        f"provenance.environment.dependency_versions['{str(key)}']"
                    ),
                )
            )
            for key, value in dependency_versions.items()
        },
        platform={
            str(key): _require_str(
                value,
                field_name=f"provenance.environment.platform['{str(key)}']",
            )
            for key, value in platform_payload.items()
        },
        blas_lapack={
            str(key): _to_json_value(value)
            for key, value in blas_lapack_payload.items()
        },
        thread_environment={
            str(key): (
                None
                if value is None
                else _require_str(
                    value,
                    field_name=f"provenance.environment.thread_environment['{str(key)}']",
                )
            )
            for key, value in thread_environment_payload.items()
        },
        timezone=_optional_str(
            payload.get("timezone"),
            field_name="provenance.environment.timezone",
        ),
        locale={
            str(key): (
                None
                if value is None
                else _require_str(
                    value,
                    field_name=f"provenance.environment.locale['{str(key)}']",
                )
            )
            for key, value in locale_payload.items()
        },
        constraints_fingerprint={
            str(key): (
                None
                if value is None
                else _require_str(
                    value,
                    field_name=(
                        f"provenance.environment.constraints_fingerprint['{str(key)}']"
                    ),
                )
            )
            for key, value in constraints_fingerprint_payload.items()
        },
    )


def _stage_to_payload(stage: PreprocessingStageProvenance) -> dict[str, object]:
    determinism = _resolve_determinism(
        stage.determinism,
        field_name="preprocessing_stage.determinism",
    )
    phospho_input_hash = _require_str(
        stage.phospho_input_hash,
        field_name="preprocessing_stage.phospho_input_hash",
    )
    phospho_output_hash = _require_str(
        stage.phospho_output_hash,
        field_name="preprocessing_stage.phospho_output_hash",
    )
    payload = {
        "stage": stage.stage,
        "operation": stage.operation,
        "parameters": _to_json_safe(stage.parameters),
        "input_shape": [int(stage.input_shape[0]), int(stage.input_shape[1])],
        "output_shape": [int(stage.output_shape[0]), int(stage.output_shape[1])],
        "input_hash": stage.input_hash,
        "output_hash": stage.output_hash,
        "phospho_input_hash": phospho_input_hash,
        "phospho_output_hash": phospho_output_hash,
        "schema_version": int(stage.schema_version),
        "consumed_input_tables": [
            _table_fingerprint_to_payload(item) for item in stage.consumed_input_tables
        ],
        "produced_output_tables": [
            _table_fingerprint_to_payload(item) for item in stage.produced_output_tables
        ],
        "backend": stage.backend,
        "random_seed": stage.random_seed,
        "determinism": determinism.value,
        "reproducibility_caveats": [
            _reproducibility_caveat_to_payload(item)
            for item in stage.reproducibility_caveats
        ],
        "dropped_row_ids": list(stage.dropped_row_ids),
        "dropped_row_count": int(stage.dropped_row_count),
        "imputed_cell_count": int(stage.imputed_cell_count),
        "imputed_row_ids": list(stage.imputed_row_ids),
        "notes": stage.notes,
        "diagnostics": (
            None if stage.diagnostics is None else _to_json_safe(stage.diagnostics)
        ),
    }
    if stage.batch_correction_provenance is not None:
        payload["batch_correction_provenance"] = batch_correction_provenance_to_payload(
            stage.batch_correction_provenance
        )
    return payload


def _stage_from_payload(payload: Mapping[str, object]) -> PreprocessingStageProvenance:
    _reject_legacy_provenance_fields(
        payload,
        field_name="preprocessing_stage",
        legacy_fields=_LEGACY_PREPROCESSING_STAGE_FIELDS,
    )
    schema_version = _require_current_stage_schema_version(payload)
    parameters = _require_mapping(
        payload.get("parameters"),
        field_name="preprocessing_stage.parameters",
    )
    input_shape = _require_shape(
        payload.get("input_shape"),
        field_name="preprocessing_stage.input_shape",
    )
    output_shape = _require_shape(
        payload.get("output_shape"),
        field_name="preprocessing_stage.output_shape",
    )
    dropped_row_ids = tuple(
        _require_str(item, field_name="preprocessing_stage.dropped_row_ids[]")
        for item in _require_sequence(
            payload.get("dropped_row_ids"),
            field_name="preprocessing_stage.dropped_row_ids",
        )
    )
    imputed_row_ids = tuple(
        _require_str(item, field_name="preprocessing_stage.imputed_row_ids[]")
        for item in _require_sequence(
            payload.get("imputed_row_ids", []),
            field_name="preprocessing_stage.imputed_row_ids",
        )
    )
    diagnostics_raw = payload.get("diagnostics")
    diagnostics = (
        None
        if diagnostics_raw is None
        else {
            str(key): _to_json_value(value)
            for key, value in _require_mapping(
                diagnostics_raw,
                field_name="preprocessing_stage.diagnostics",
            ).items()
        }
    )
    consumed_input_tables = _table_fingerprints_from_payload(
        payload.get("consumed_input_tables"),
        field_name="preprocessing_stage.consumed_input_tables",
    )
    produced_output_tables = _table_fingerprints_from_payload(
        payload.get("produced_output_tables"),
        field_name="preprocessing_stage.produced_output_tables",
    )
    random_seed = _optional_int(
        payload.get("random_seed"),
        field_name="preprocessing_stage.random_seed",
    )
    determinism = _resolve_determinism(
        payload.get("determinism"),
        field_name="preprocessing_stage.determinism",
    )
    reproducibility_caveats = tuple(
        _reproducibility_caveat_from_payload(
            _require_mapping(
                item,
                field_name=(f"preprocessing_stage.reproducibility_caveats[{position}]"),
            )
        )
        for position, item in enumerate(
            _require_sequence(
                payload.get("reproducibility_caveats", []),
                field_name="preprocessing_stage.reproducibility_caveats",
            )
        )
    )
    input_hash = _require_str(
        payload.get("input_hash"),
        field_name="preprocessing_stage.input_hash",
    )
    output_hash = _require_str(
        payload.get("output_hash"),
        field_name="preprocessing_stage.output_hash",
    )
    phospho_input_hash = _require_str(
        payload.get("phospho_input_hash"),
        field_name="preprocessing_stage.phospho_input_hash",
    )
    phospho_output_hash = _require_str(
        payload.get("phospho_output_hash"),
        field_name="preprocessing_stage.phospho_output_hash",
    )
    batch_correction_provenance_raw = payload.get("batch_correction_provenance")
    batch_correction_provenance = (
        None
        if batch_correction_provenance_raw is None
        else batch_correction_provenance_from_payload(
            _require_mapping(
                batch_correction_provenance_raw,
                field_name="preprocessing_stage.batch_correction_provenance",
            )
        )
    )
    return PreprocessingStageProvenance(
        stage=_require_str(
            payload.get("stage"), field_name="preprocessing_stage.stage"
        ),
        operation=_require_str(
            payload.get("operation"),
            field_name="preprocessing_stage.operation",
        ),
        parameters={
            str(key): _to_json_value(value) for key, value in parameters.items()
        },
        input_shape=input_shape,
        output_shape=output_shape,
        input_hash=input_hash,
        output_hash=output_hash,
        phospho_input_hash=phospho_input_hash,
        phospho_output_hash=phospho_output_hash,
        dropped_row_ids=dropped_row_ids,
        dropped_row_count=_require_int(
            payload.get("dropped_row_count"),
            field_name="preprocessing_stage.dropped_row_count",
        ),
        schema_version=schema_version,
        consumed_input_tables=consumed_input_tables,
        produced_output_tables=produced_output_tables,
        backend=_optional_str(
            payload.get("backend"),
            field_name="preprocessing_stage.backend",
        ),
        random_seed=random_seed,
        determinism=determinism,
        reproducibility_caveats=reproducibility_caveats,
        imputed_cell_count=_require_int(
            payload.get("imputed_cell_count", 0),
            field_name="preprocessing_stage.imputed_cell_count",
        ),
        imputed_row_ids=imputed_row_ids,
        notes=_optional_str(
            payload.get("notes"),
            field_name="preprocessing_stage.notes",
        ),
        diagnostics=diagnostics,
        batch_correction_provenance=batch_correction_provenance,
    )


def _reference_to_payload(reference: ReferenceProvenance) -> dict[str, object]:
    return {
        "source_type": reference.source_type,
        "organism": reference.organism,
        "bundle_id": reference.bundle_id,
        "source_name": reference.source_name,
        "source_version": reference.source_version,
        "retrieved_at": reference.retrieved_at,
        "identifier_namespace": reference.identifier_namespace,
        "sequence_window": (
            None
            if reference.sequence_window is None
            else _to_json_safe(reference.sequence_window)
        ),
        "manifest": (
            None if reference.manifest is None else _to_json_safe(reference.manifest)
        ),
        "reference_context": (
            None
            if reference.reference_context is None
            else reference.reference_context.to_payload()
        ),
        "table_fingerprints": [
            _table_fingerprint_to_payload(item) for item in reference.table_fingerprints
        ],
        "identifier_normalisation": (
            None
            if reference.identifier_normalisation is None
            else _reference_identifier_normalisation_to_payload(
                reference.identifier_normalisation
            )
        ),
    }


def _reference_from_payload(payload: Mapping[str, object]) -> ReferenceProvenance:
    table_payload = _require_sequence(
        payload.get("table_fingerprints"),
        field_name="reference_provenance.table_fingerprints",
    )
    sequence_window_payload = _optional_mapping(
        payload.get("sequence_window"),
        field_name="reference_provenance.sequence_window",
    )
    manifest_payload = _optional_mapping(
        payload.get("manifest"),
        field_name="reference_provenance.manifest",
    )
    reference_context_payload = _optional_mapping(
        payload.get("reference_context"),
        field_name="reference_provenance.reference_context",
    )
    return ReferenceProvenance(
        source_type=_require_str(
            payload.get("source_type"),
            field_name="reference_provenance.source_type",
        ),
        organism=_require_str(
            payload.get("organism"),
            field_name="reference_provenance.organism",
        ),
        bundle_id=_optional_str(
            payload.get("bundle_id"),
            field_name="reference_provenance.bundle_id",
        ),
        source_name=_optional_str(
            payload.get("source_name"),
            field_name="reference_provenance.source_name",
        ),
        source_version=_optional_str(
            payload.get("source_version"),
            field_name="reference_provenance.source_version",
        ),
        retrieved_at=_optional_str(
            payload.get("retrieved_at"),
            field_name="reference_provenance.retrieved_at",
        ),
        identifier_namespace=_optional_str(
            payload.get("identifier_namespace"),
            field_name="reference_provenance.identifier_namespace",
        ),
        sequence_window=None
        if sequence_window_payload is None
        else {
            str(key): _to_json_value(value)
            for key, value in sequence_window_payload.items()
        },
        manifest=None
        if manifest_payload is None
        else {
            str(key): _to_json_value(value) for key, value in manifest_payload.items()
        },
        reference_context=None
        if reference_context_payload is None
        else ReferenceContext.from_payload(reference_context_payload),
        table_fingerprints=tuple(
            _table_fingerprint_from_payload(
                _require_mapping(
                    item,
                    field_name=f"reference_provenance.table_fingerprints[{position}]",
                )
            )
            for position, item in enumerate(table_payload)
        ),
        identifier_normalisation=_optional_reference_identifier_normalisation_from_payload(
            payload.get("identifier_normalisation"),
            field_name="reference_provenance.identifier_normalisation",
        ),
    )


def _reference_identifier_normalisation_to_payload(
    report: ReferenceIdentifierNormalisationReport,
) -> dict[str, object]:
    return {
        "schema_version": int(report.schema_version),
        "original_row_count": int(report.original_row_count),
        "normalised_row_count": int(report.normalised_row_count),
        "invalid_identifier_count": int(report.invalid_identifier_count),
        "changed_identifier_count": int(report.changed_identifier_count),
        "duplicate_identifier_count": int(report.duplicate_identifier_count),
        "conflict_count": int(report.conflict_count),
        "records": [
            _reference_identifier_normalisation_record_to_payload(record)
            for record in report.records
        ],
    }


def _reference_identifier_normalisation_record_to_payload(
    record: ReferenceIdentifierNormalisationRecord,
) -> dict[str, object]:
    return {
        "table_name": record.table_name,
        "column_name": record.column_name,
        "row_position": int(record.row_position),
        "identifier_kind": record.identifier_kind,
        "original_value": record.original_value,
        "normalised_value": record.normalised_value,
        "status": record.status,
        "reason": record.reason,
    }


def _optional_reference_identifier_normalisation_from_payload(
    value: object,
    *,
    field_name: str,
) -> ReferenceIdentifierNormalisationReport | None:
    if value is None:
        return None
    return _reference_identifier_normalisation_from_payload(
        _require_mapping(value, field_name=field_name)
    )


def _reference_identifier_normalisation_from_payload(
    payload: Mapping[str, object],
) -> ReferenceIdentifierNormalisationReport:
    from phospy.science.references.identifiers import (
        ReferenceIdentifierNormalisationReport,
    )

    records_payload = _require_sequence(
        payload.get("records"),
        field_name="reference_identifier_normalisation.records",
    )
    return ReferenceIdentifierNormalisationReport(
        schema_version=_require_int(
            payload.get("schema_version"),
            field_name="reference_identifier_normalisation.schema_version",
        ),
        original_row_count=_require_int(
            payload.get("original_row_count"),
            field_name="reference_identifier_normalisation.original_row_count",
        ),
        normalised_row_count=_require_int(
            payload.get("normalised_row_count"),
            field_name="reference_identifier_normalisation.normalised_row_count",
        ),
        invalid_identifier_count=_require_int(
            payload.get("invalid_identifier_count"),
            field_name="reference_identifier_normalisation.invalid_identifier_count",
        ),
        changed_identifier_count=_require_int(
            payload.get("changed_identifier_count"),
            field_name="reference_identifier_normalisation.changed_identifier_count",
        ),
        duplicate_identifier_count=_require_int(
            payload.get("duplicate_identifier_count"),
            field_name=(
                "reference_identifier_normalisation.duplicate_identifier_count"
            ),
        ),
        conflict_count=_require_int(
            payload.get("conflict_count"),
            field_name="reference_identifier_normalisation.conflict_count",
        ),
        records=tuple(
            _reference_identifier_normalisation_record_from_payload(
                _require_mapping(
                    item,
                    field_name=(
                        f"reference_identifier_normalisation.records[{position}]"
                    ),
                )
            )
            for position, item in enumerate(records_payload)
        ),
    )


def _reference_identifier_normalisation_record_from_payload(
    payload: Mapping[str, object],
) -> ReferenceIdentifierNormalisationRecord:
    from phospy.science.references.identifiers import (
        ReferenceIdentifierNormalisationRecord,
    )

    return ReferenceIdentifierNormalisationRecord(
        table_name=_require_str(
            payload.get("table_name"),
            field_name="reference_identifier_normalisation_record.table_name",
        ),
        column_name=_require_str(
            payload.get("column_name"),
            field_name="reference_identifier_normalisation_record.column_name",
        ),
        row_position=_require_int(
            payload.get("row_position"),
            field_name="reference_identifier_normalisation_record.row_position",
        ),
        identifier_kind=_require_str(
            payload.get("identifier_kind"),
            field_name="reference_identifier_normalisation_record.identifier_kind",
        ),
        original_value=_require_raw_str(
            payload.get("original_value"),
            field_name="reference_identifier_normalisation_record.original_value",
        ),
        normalised_value=_optional_raw_str(
            payload.get("normalised_value"),
            field_name="reference_identifier_normalisation_record.normalised_value",
        ),
        status=_require_str(
            payload.get("status"),
            field_name="reference_identifier_normalisation_record.status",
        ),
        reason=_optional_str(
            payload.get("reason"),
            field_name="reference_identifier_normalisation_record.reason",
        ),
    )


def _to_json_safe(value: object) -> object:
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, object], value)
        return {str(key): _to_json_safe(item) for key, item in mapping.items()}
    if isinstance(value, tuple):
        return [_to_json_safe(item) for item in cast(tuple[object, ...], value)]
    if isinstance(value, list):
        return [_to_json_safe(item) for item in cast(list[object], value)]
    return value


def _to_json_value(value: object) -> JsonValue:
    return cast(JsonValue, _to_json_safe(value))


def _require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, object], value)
        return {str(key): item for key, item in mapping.items()}
    raise PhosPyInputError(f"{field_name} must be an object")


def _require_sequence(value: object, *, field_name: str) -> list[object]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(cast(Sequence[object], value))
    raise PhosPyInputError(f"{field_name} must be an array")


def _require_str(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def _optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_str(value, field_name=field_name)


def _require_raw_str(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    return value


def _optional_raw_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_raw_str(value, field_name=field_name)


def _optional_mapping(value: object, *, field_name: str) -> Mapping[str, object] | None:
    if value is None:
        return None
    return _require_mapping(value, field_name=field_name)


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an int")
    return int(value)


def _optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return _require_int(value, field_name=field_name)


def _require_current_environment_schema_version(
    payload: Mapping[str, object],
) -> int:
    if "schema_version" not in payload:
        _raise_legacy_provenance_schema()
    schema_version = _require_int(
        payload.get("schema_version"),
        field_name="provenance.environment.schema_version",
    )
    if schema_version != ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2:
        _raise_legacy_provenance_schema()
    return schema_version


def _require_current_stage_schema_version(payload: Mapping[str, object]) -> int:
    if "schema_version" not in payload:
        _raise_legacy_provenance_schema()
    schema_version = _require_int(
        payload.get("schema_version"),
        field_name="preprocessing_stage.schema_version",
    )
    if schema_version != PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3:
        _raise_legacy_provenance_schema()
    return schema_version


def _table_fingerprints_from_payload(
    value: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if value is None:
        _raise_legacy_provenance_schema()
    payload = _require_sequence(value, field_name=field_name)
    return tuple(
        _table_fingerprint_from_payload(
            _require_mapping(item, field_name=f"{field_name}[{position}]")
        )
        for position, item in enumerate(payload)
    )


def _reproducibility_caveat_to_payload(
    caveat: ReproducibilityCaveat,
) -> dict[str, object]:
    return caveat.to_payload()


def _reproducibility_caveat_from_payload(
    payload: Mapping[str, object],
) -> ReproducibilityCaveat:
    details = _require_mapping(
        payload.get("details", {}),
        field_name="reproducibility_caveat.details",
    )
    return ReproducibilityCaveat(
        code=_require_str(
            payload.get("code"),
            field_name="reproducibility_caveat.code",
        ),
        severity=_require_str(
            payload.get("severity"),
            field_name="reproducibility_caveat.severity",
        ),
        message=_require_str(
            payload.get("message"),
            field_name="reproducibility_caveat.message",
        ),
        details={str(key): _to_json_value(value) for key, value in details.items()},
    )


def _resolve_determinism(value: object, *, field_name: str) -> DeterminismKind:
    if isinstance(value, DeterminismKind):
        return value
    normalized = _require_str(value, field_name=field_name)
    alias = _DETERMINISM_ALIASES.get(normalized)
    if alias is not None:
        return alias
    try:
        return DeterminismKind(normalized)
    except ValueError as exc:
        raise PhosPyInputError(
            f"{field_name} must be one of: "
            f"{PREPROCESSING_STAGE_DETERMINISM_PURE!r}, "
            f"{PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC!r}, "
            f"{PREPROCESSING_STAGE_DETERMINISM_EXTERNALLY_NONDETERMINISTIC!r}"
        ) from exc


def _reject_legacy_provenance_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    legacy_fields: frozenset[str],
) -> None:
    present = sorted(str(key) for key in legacy_fields if str(key) in payload)
    if present:
        _raise_legacy_provenance_schema()


def _raise_legacy_provenance_schema() -> None:
    raise PhosPyInputError(_LEGACY_PROVENANCE_SCHEMA_ERROR)


def _require_shape(value: object, *, field_name: str) -> tuple[int, int]:
    sequence = _require_sequence(value, field_name=field_name)
    if len(sequence) != 2:
        raise PhosPyInputError(f"{field_name} must contain exactly two integers")
    return (
        _require_int(sequence[0], field_name=f"{field_name}[0]"),
        _require_int(sequence[1], field_name=f"{field_name}[1]"),
    )


__all__ = [
    "batch_correction_provenance_from_payload",
    "batch_correction_provenance_to_payload",
    "from_payload",
    "table_fingerprint_from_payload",
    "table_fingerprint_to_payload",
    "to_payload",
]
