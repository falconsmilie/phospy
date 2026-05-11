"""Serialization helpers for provenance payloads."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import (
    ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V1,
    PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY,
    PREPROCESSING_STAGE_DETERMINISM_PURE,
    PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC,
    EnvironmentProvenance,
    JsonValue,
    PreprocessingStageProvenance,
    ReferenceProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.scientific_policies import ScientificPolicyRecord

if TYPE_CHECKING:
    from phospy.references.identifiers import (
        ReferenceIdentifierNormalisationRecord,
        ReferenceIdentifierNormalisationReport,
    )


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
    )


def _table_fingerprint_to_payload(fingerprint: TableFingerprint) -> dict[str, object]:
    legacy_hash_algorithm = (
        fingerprint.tolerance_hash_algorithm
        if fingerprint.hash_algorithm is None
        else fingerprint.hash_algorithm
    )
    legacy_hash_value = (
        fingerprint.tolerance_hash_value
        if fingerprint.hash_value is None
        else fingerprint.hash_value
    )
    exact_hash_algorithm = (
        fingerprint.hash_algorithm
        if fingerprint.exact_hash_algorithm is None
        else fingerprint.exact_hash_algorithm
    )
    exact_hash_value = (
        fingerprint.hash_value
        if fingerprint.exact_hash_value is None
        else fingerprint.exact_hash_value
    )
    tolerance_hash_algorithm = (
        fingerprint.hash_algorithm
        if fingerprint.tolerance_hash_algorithm is None
        else fingerprint.tolerance_hash_algorithm
    )
    tolerance_hash_value = (
        fingerprint.hash_value
        if fingerprint.tolerance_hash_value is None
        else fingerprint.tolerance_hash_value
    )
    return {
        "name": fingerprint.name,
        "rows": int(fingerprint.rows),
        "columns": int(fingerprint.columns),
        "index_name": fingerprint.index_name,
        "column_names": list(fingerprint.column_names),
        "dtypes": list(fingerprint.dtypes),
        # Compatibility aliases retained for existing consumers.
        "hash_algorithm": legacy_hash_algorithm,
        "hash_value": legacy_hash_value,
        "exact_hash_algorithm": exact_hash_algorithm,
        "exact_hash_value": exact_hash_value,
        "tolerance_hash_algorithm": tolerance_hash_algorithm,
        "tolerance_hash_value": tolerance_hash_value,
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


def _table_fingerprint_from_payload(payload: Mapping[str, object]) -> TableFingerprint:
    index_structure = _optional_mapping(
        payload.get("index_structure"),
        field_name="table_fingerprint.index_structure",
    )
    column_index_structure = _optional_mapping(
        payload.get("column_index_structure"),
        field_name="table_fingerprint.column_index_structure",
    )
    hash_algorithm = _optional_str(
        payload.get("hash_algorithm"),
        field_name="table_fingerprint.hash_algorithm",
    )
    hash_value = _optional_str(
        payload.get("hash_value"),
        field_name="table_fingerprint.hash_value",
    )
    exact_hash_algorithm = _optional_str(
        payload.get("exact_hash_algorithm"),
        field_name="table_fingerprint.exact_hash_algorithm",
    )
    exact_hash_value = _optional_str(
        payload.get("exact_hash_value"),
        field_name="table_fingerprint.exact_hash_value",
    )
    tolerance_hash_algorithm = _optional_str(
        payload.get("tolerance_hash_algorithm"),
        field_name="table_fingerprint.tolerance_hash_algorithm",
    )
    tolerance_hash_value = _optional_str(
        payload.get("tolerance_hash_value"),
        field_name="table_fingerprint.tolerance_hash_value",
    )
    resolved_tolerance_hash_algorithm = _resolve_fingerprint_hash_field(
        primary=tolerance_hash_algorithm,
        fallback=(
            hash_algorithm if hash_algorithm is not None else exact_hash_algorithm
        ),
        field_name="table_fingerprint.tolerance_hash_algorithm",
    )
    resolved_tolerance_hash_value = _resolve_fingerprint_hash_field(
        primary=tolerance_hash_value,
        fallback=(hash_value if hash_value is not None else exact_hash_value),
        field_name="table_fingerprint.tolerance_hash_value",
    )
    resolved_exact_hash_algorithm = _resolve_fingerprint_hash_field(
        primary=exact_hash_algorithm,
        fallback=(
            hash_algorithm if hash_algorithm is not None else tolerance_hash_algorithm
        ),
        field_name="table_fingerprint.exact_hash_algorithm",
    )
    resolved_exact_hash_value = _resolve_fingerprint_hash_field(
        primary=exact_hash_value,
        fallback=(hash_value if hash_value is not None else tolerance_hash_value),
        field_name="table_fingerprint.exact_hash_value",
    )
    resolved_hash_algorithm = _resolve_fingerprint_hash_field(
        primary=hash_algorithm,
        fallback=(
            tolerance_hash_algorithm
            if tolerance_hash_algorithm is not None
            else exact_hash_algorithm
        ),
        field_name="table_fingerprint.hash_algorithm",
    )
    resolved_hash_value = _resolve_fingerprint_hash_field(
        primary=hash_value,
        fallback=(
            tolerance_hash_value
            if tolerance_hash_value is not None
            else exact_hash_value
        ),
        field_name="table_fingerprint.hash_value",
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
        hash_algorithm=resolved_hash_algorithm,
        hash_value=resolved_hash_value,
        exact_hash_algorithm=resolved_exact_hash_algorithm,
        exact_hash_value=resolved_exact_hash_value,
        tolerance_hash_algorithm=resolved_tolerance_hash_algorithm,
        tolerance_hash_value=resolved_tolerance_hash_value,
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


def _resolve_fingerprint_hash_field(
    *,
    primary: str | None,
    fallback: str | None,
    field_name: str,
) -> str:
    if primary is not None:
        return primary
    if fallback is not None:
        return fallback
    raise PhosPyInputError(f"{field_name} must be provided")


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
        schema_version=_require_int(
            payload.get(
                "schema_version",
                ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V1,
            ),
            field_name="provenance.environment.schema_version",
        ),
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
    return {
        "stage": stage.stage,
        "operation": stage.operation,
        "parameters": _to_json_safe(stage.parameters),
        "input_shape": [int(stage.input_shape[0]), int(stage.input_shape[1])],
        "output_shape": [int(stage.output_shape[0]), int(stage.output_shape[1])],
        "input_hash": stage.input_hash,
        "output_hash": stage.output_hash,
        "phospho_input_hash": (
            stage.input_hash
            if stage.phospho_input_hash is None
            else stage.phospho_input_hash
        ),
        "phospho_output_hash": (
            stage.output_hash
            if stage.phospho_output_hash is None
            else stage.phospho_output_hash
        ),
        "schema_version": int(stage.schema_version),
        "consumed_input_tables": [
            _table_fingerprint_to_payload(item) for item in stage.consumed_input_tables
        ],
        "produced_output_tables": [
            _table_fingerprint_to_payload(item) for item in stage.produced_output_tables
        ],
        "backend": stage.backend,
        "random_seed": stage.random_seed,
        "determinism": determinism,
        "is_deterministic": _determinism_is_deterministic(determinism),
        "dropped_row_ids": list(stage.dropped_row_ids),
        "dropped_row_count": int(stage.dropped_row_count),
        "imputed_cell_count": int(stage.imputed_cell_count),
        "imputed_row_ids": list(stage.imputed_row_ids),
        "notes": stage.notes,
        "diagnostics": (
            None if stage.diagnostics is None else _to_json_safe(stage.diagnostics)
        ),
    }


def _stage_from_payload(payload: Mapping[str, object]) -> PreprocessingStageProvenance:
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
    consumed_input_tables = _optional_table_fingerprints(
        payload.get("consumed_input_tables"),
        field_name="preprocessing_stage.consumed_input_tables",
    )
    produced_output_tables = _optional_table_fingerprints(
        payload.get("produced_output_tables"),
        field_name="preprocessing_stage.produced_output_tables",
    )
    random_seed = _optional_int(
        payload.get("random_seed"),
        field_name="preprocessing_stage.random_seed",
    )
    deterministic_alias_raw = payload.get("is_deterministic")
    deterministic_alias = (
        None
        if deterministic_alias_raw is None
        else _require_bool(
            deterministic_alias_raw,
            field_name="preprocessing_stage.is_deterministic",
        )
    )
    determinism = _resolve_stage_determinism_from_payload(
        payload=payload,
        random_seed=random_seed,
        is_deterministic=deterministic_alias,
    )
    is_deterministic = _determinism_is_deterministic(determinism)
    input_hash = _resolve_primary_hash(
        payload=payload,
        primary_field_name="input_hash",
        alias_field_name="phospho_input_hash",
    )
    output_hash = _resolve_primary_hash(
        payload=payload,
        primary_field_name="output_hash",
        alias_field_name="phospho_output_hash",
    )
    phospho_input_hash = _optional_str(
        payload.get("phospho_input_hash"),
        field_name="preprocessing_stage.phospho_input_hash",
    )
    phospho_output_hash = _optional_str(
        payload.get("phospho_output_hash"),
        field_name="preprocessing_stage.phospho_output_hash",
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
        phospho_input_hash=(
            input_hash if phospho_input_hash is None else phospho_input_hash
        ),
        phospho_output_hash=(
            output_hash if phospho_output_hash is None else phospho_output_hash
        ),
        dropped_row_ids=dropped_row_ids,
        dropped_row_count=_require_int(
            payload.get("dropped_row_count"),
            field_name="preprocessing_stage.dropped_row_count",
        ),
        schema_version=_require_int(
            payload.get("schema_version", 1),
            field_name="preprocessing_stage.schema_version",
        ),
        consumed_input_tables=consumed_input_tables,
        produced_output_tables=produced_output_tables,
        backend=_optional_str(
            payload.get("backend"),
            field_name="preprocessing_stage.backend",
        ),
        random_seed=random_seed,
        determinism=determinism,
        is_deterministic=is_deterministic,
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
    from phospy.references.identifiers import ReferenceIdentifierNormalisationReport

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
    from phospy.references.identifiers import ReferenceIdentifierNormalisationRecord

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


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return bool(value)


def _optional_table_fingerprints(
    value: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if value is None:
        return ()
    payload = _require_sequence(value, field_name=field_name)
    return tuple(
        _table_fingerprint_from_payload(
            _require_mapping(item, field_name=f"{field_name}[{position}]")
        )
        for position, item in enumerate(payload)
    )


def _resolve_primary_hash(
    *,
    payload: Mapping[str, object],
    primary_field_name: str,
    alias_field_name: str,
) -> str:
    raw_value = payload.get(primary_field_name)
    if raw_value is None:
        raw_value = payload.get(alias_field_name)
        if raw_value is None:
            raise PhosPyInputError(
                f"preprocessing_stage.{primary_field_name} is required"
            )
    return _require_str(
        raw_value, field_name=f"preprocessing_stage.{primary_field_name}"
    )


def _resolve_stage_determinism_from_payload(
    *,
    payload: Mapping[str, object],
    random_seed: int | None,
    is_deterministic: bool | None,
) -> str:
    determinism_raw = payload.get("determinism")
    if determinism_raw is None:
        if random_seed is not None:
            return PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC
        if is_deterministic is False:
            return PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY
        return PREPROCESSING_STAGE_DETERMINISM_PURE
    return _resolve_determinism(
        determinism_raw,
        field_name="preprocessing_stage.determinism",
    )


def _resolve_determinism(value: object, *, field_name: str) -> str:
    normalized = _require_str(value, field_name=field_name)
    if normalized not in {
        PREPROCESSING_STAGE_DETERMINISM_PURE,
        PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC,
        PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY,
    }:
        raise PhosPyInputError(
            f"{field_name} must be one of: "
            f"{PREPROCESSING_STAGE_DETERMINISM_PURE!r}, "
            f"{PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC!r}, "
            f"{PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY!r}"
        )
    return normalized


def _determinism_is_deterministic(determinism: str) -> bool:
    return determinism == PREPROCESSING_STAGE_DETERMINISM_PURE


def _require_shape(value: object, *, field_name: str) -> tuple[int, int]:
    sequence = _require_sequence(value, field_name=field_name)
    if len(sequence) != 2:
        raise PhosPyInputError(f"{field_name} must contain exactly two integers")
    return (
        _require_int(sequence[0], field_name=f"{field_name}[0]"),
        _require_int(sequence[1], field_name=f"{field_name}[1]"),
    )


__all__ = ["from_payload", "to_payload"]
