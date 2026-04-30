"""Serialization helpers for provenance payloads."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import (
    EnvironmentProvenance,
    PreprocessingStageProvenance,
    ReferenceProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.scientific_policies import ScientificPolicyRecord


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
            str(key): value for key, value in workflow_parameters.items()
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
    return {
        "name": fingerprint.name,
        "rows": int(fingerprint.rows),
        "columns": int(fingerprint.columns),
        "index_name": fingerprint.index_name,
        "column_names": list(fingerprint.column_names),
        "dtypes": list(fingerprint.dtypes),
        "hash_algorithm": fingerprint.hash_algorithm,
        "hash_value": fingerprint.hash_value,
    }


def _table_fingerprint_from_payload(payload: Mapping[str, object]) -> TableFingerprint:
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
        hash_algorithm=_require_str(
            payload.get("hash_algorithm"),
            field_name="table_fingerprint.hash_algorithm",
        ),
        hash_value=_require_str(
            payload.get("hash_value"),
            field_name="table_fingerprint.hash_value",
        ),
    )


def _environment_to_payload(environment: EnvironmentProvenance) -> dict[str, object]:
    return {
        "package_name": environment.package_name,
        "package_version": environment.package_version,
        "python_version": environment.python_version,
        "dependency_versions": _to_json_safe(environment.dependency_versions),
        "platform": _to_json_safe(environment.platform),
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
    return EnvironmentProvenance(
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
    )


def _stage_to_payload(stage: PreprocessingStageProvenance) -> dict[str, object]:
    return {
        "stage": stage.stage,
        "operation": stage.operation,
        "parameters": _to_json_safe(stage.parameters),
        "input_shape": [int(stage.input_shape[0]), int(stage.input_shape[1])],
        "output_shape": [int(stage.output_shape[0]), int(stage.output_shape[1])],
        "input_hash": stage.input_hash,
        "output_hash": stage.output_hash,
        "schema_version": int(stage.schema_version),
        "consumed_input_tables": [
            _table_fingerprint_to_payload(item) for item in stage.consumed_input_tables
        ],
        "produced_output_tables": [
            _table_fingerprint_to_payload(item) for item in stage.produced_output_tables
        ],
        "backend": stage.backend,
        "random_seed": stage.random_seed,
        "is_deterministic": bool(stage.is_deterministic),
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
            str(key): value
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
    return PreprocessingStageProvenance(
        stage=_require_str(
            payload.get("stage"), field_name="preprocessing_stage.stage"
        ),
        operation=_require_str(
            payload.get("operation"),
            field_name="preprocessing_stage.operation",
        ),
        parameters={str(key): value for key, value in parameters.items()},
        input_shape=input_shape,
        output_shape=output_shape,
        input_hash=_require_str(
            payload.get("input_hash"),
            field_name="preprocessing_stage.input_hash",
        ),
        output_hash=_require_str(
            payload.get("output_hash"),
            field_name="preprocessing_stage.output_hash",
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
        random_seed=_optional_int(
            payload.get("random_seed"),
            field_name="preprocessing_stage.random_seed",
        ),
        is_deterministic=_require_bool(
            payload.get("is_deterministic", True),
            field_name="preprocessing_stage.is_deterministic",
        ),
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
        "table_fingerprints": [
            _table_fingerprint_to_payload(item) for item in reference.table_fingerprints
        ],
    }


def _reference_from_payload(payload: Mapping[str, object]) -> ReferenceProvenance:
    table_payload = _require_sequence(
        payload.get("table_fingerprints"),
        field_name="reference_provenance.table_fingerprints",
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
        table_fingerprints=tuple(
            _table_fingerprint_from_payload(
                _require_mapping(
                    item,
                    field_name=f"reference_provenance.table_fingerprints[{position}]",
                )
            )
            for position, item in enumerate(table_payload)
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


def _require_shape(value: object, *, field_name: str) -> tuple[int, int]:
    sequence = _require_sequence(value, field_name=field_name)
    if len(sequence) != 2:
        raise PhosPyInputError(f"{field_name} must contain exactly two integers")
    return (
        _require_int(sequence[0], field_name=f"{field_name}[0]"),
        _require_int(sequence[1], field_name=f"{field_name}[1]"),
    )


__all__ = ["from_payload", "to_payload"]
