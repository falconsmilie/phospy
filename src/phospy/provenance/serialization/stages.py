"""Preprocessing-stage provenance payload serialization."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import (
    PREPROCESSING_STAGE_DETERMINISM_EXTERNALLY_NONDETERMINISTIC,
    PREPROCESSING_STAGE_DETERMINISM_PURE,
    PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    DeterminismKind,
    PreprocessingStageProvenance,
    ReproducibilityCaveat,
)
from phospy.provenance.serialization._payload import (
    optional_int,
    optional_str,
    raise_legacy_provenance_schema,
    reject_legacy_provenance_fields,
    require_int,
    require_mapping,
    require_sequence,
    require_shape,
    require_str,
    to_json_safe,
    to_json_value,
)
from phospy.provenance.serialization.batch_correction import (
    batch_correction_provenance_from_payload,
    batch_correction_provenance_to_payload,
)
from phospy.provenance.serialization.tables import (
    table_fingerprint_to_payload,
    table_fingerprints_from_payload,
)

_LEGACY_PREPROCESSING_STAGE_FIELDS = frozenset({"is_deterministic"})

_DETERMINISM_ALIASES = {
    "pure": DeterminismKind.DETERMINISTIC,
    "external_dependency": DeterminismKind.EXTERNALLY_NONDETERMINISTIC,
}


def stage_to_payload(stage: PreprocessingStageProvenance) -> dict[str, object]:
    determinism = _resolve_determinism(
        stage.determinism,
        field_name="preprocessing_stage.determinism",
    )
    phospho_input_hash = require_str(
        stage.phospho_input_hash,
        field_name="preprocessing_stage.phospho_input_hash",
    )
    phospho_output_hash = require_str(
        stage.phospho_output_hash,
        field_name="preprocessing_stage.phospho_output_hash",
    )
    payload = {
        "stage": stage.stage,
        "operation": stage.operation,
        "parameters": to_json_safe(stage.parameters),
        "input_shape": [int(stage.input_shape[0]), int(stage.input_shape[1])],
        "output_shape": [int(stage.output_shape[0]), int(stage.output_shape[1])],
        "input_hash": stage.input_hash,
        "output_hash": stage.output_hash,
        "phospho_input_hash": phospho_input_hash,
        "phospho_output_hash": phospho_output_hash,
        "schema_version": int(stage.schema_version),
        "consumed_input_tables": [
            table_fingerprint_to_payload(item) for item in stage.consumed_input_tables
        ],
        "produced_output_tables": [
            table_fingerprint_to_payload(item) for item in stage.produced_output_tables
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
            None if stage.diagnostics is None else to_json_safe(stage.diagnostics)
        ),
    }
    if stage.batch_correction_provenance is not None:
        payload["batch_correction_provenance"] = batch_correction_provenance_to_payload(
            stage.batch_correction_provenance
        )
    return payload


def stage_from_payload(payload: Mapping[str, object]) -> PreprocessingStageProvenance:
    reject_legacy_provenance_fields(
        payload,
        field_name="preprocessing_stage",
        legacy_fields=_LEGACY_PREPROCESSING_STAGE_FIELDS,
    )
    schema_version = _require_current_stage_schema_version(payload)
    parameters = require_mapping(
        payload.get("parameters"),
        field_name="preprocessing_stage.parameters",
    )
    input_shape = require_shape(
        payload.get("input_shape"),
        field_name="preprocessing_stage.input_shape",
    )
    output_shape = require_shape(
        payload.get("output_shape"),
        field_name="preprocessing_stage.output_shape",
    )
    dropped_row_ids = tuple(
        require_str(item, field_name="preprocessing_stage.dropped_row_ids[]")
        for item in require_sequence(
            payload.get("dropped_row_ids"),
            field_name="preprocessing_stage.dropped_row_ids",
        )
    )
    imputed_row_ids = tuple(
        require_str(item, field_name="preprocessing_stage.imputed_row_ids[]")
        for item in require_sequence(
            payload.get("imputed_row_ids", []),
            field_name="preprocessing_stage.imputed_row_ids",
        )
    )
    diagnostics_raw = payload.get("diagnostics")
    diagnostics = (
        None
        if diagnostics_raw is None
        else {
            key: to_json_value(value)
            for key, value in require_mapping(
                diagnostics_raw,
                field_name="preprocessing_stage.diagnostics",
            ).items()
        }
    )
    consumed_input_tables = table_fingerprints_from_payload(
        payload.get("consumed_input_tables"),
        field_name="preprocessing_stage.consumed_input_tables",
    )
    produced_output_tables = table_fingerprints_from_payload(
        payload.get("produced_output_tables"),
        field_name="preprocessing_stage.produced_output_tables",
    )
    random_seed = optional_int(
        payload.get("random_seed"),
        field_name="preprocessing_stage.random_seed",
    )
    determinism = _resolve_determinism(
        payload.get("determinism"),
        field_name="preprocessing_stage.determinism",
    )
    reproducibility_caveats = tuple(
        _reproducibility_caveat_from_payload(
            require_mapping(
                item,
                field_name=(f"preprocessing_stage.reproducibility_caveats[{position}]"),
            )
        )
        for position, item in enumerate(
            require_sequence(
                payload.get("reproducibility_caveats", []),
                field_name="preprocessing_stage.reproducibility_caveats",
            )
        )
    )
    input_hash = require_str(
        payload.get("input_hash"),
        field_name="preprocessing_stage.input_hash",
    )
    output_hash = require_str(
        payload.get("output_hash"),
        field_name="preprocessing_stage.output_hash",
    )
    phospho_input_hash = require_str(
        payload.get("phospho_input_hash"),
        field_name="preprocessing_stage.phospho_input_hash",
    )
    phospho_output_hash = require_str(
        payload.get("phospho_output_hash"),
        field_name="preprocessing_stage.phospho_output_hash",
    )
    batch_correction_provenance_raw = payload.get("batch_correction_provenance")
    batch_correction_provenance = (
        None
        if batch_correction_provenance_raw is None
        else batch_correction_provenance_from_payload(
            require_mapping(
                batch_correction_provenance_raw,
                field_name="preprocessing_stage.batch_correction_provenance",
            )
        )
    )
    return PreprocessingStageProvenance(
        stage=require_str(payload.get("stage"), field_name="preprocessing_stage.stage"),
        operation=require_str(
            payload.get("operation"),
            field_name="preprocessing_stage.operation",
        ),
        parameters={key: to_json_value(value) for key, value in parameters.items()},
        input_shape=input_shape,
        output_shape=output_shape,
        input_hash=input_hash,
        output_hash=output_hash,
        phospho_input_hash=phospho_input_hash,
        phospho_output_hash=phospho_output_hash,
        dropped_row_ids=dropped_row_ids,
        dropped_row_count=require_int(
            payload.get("dropped_row_count"),
            field_name="preprocessing_stage.dropped_row_count",
        ),
        schema_version=schema_version,
        consumed_input_tables=consumed_input_tables,
        produced_output_tables=produced_output_tables,
        backend=optional_str(
            payload.get("backend"),
            field_name="preprocessing_stage.backend",
        ),
        random_seed=random_seed,
        determinism=determinism,
        reproducibility_caveats=reproducibility_caveats,
        imputed_cell_count=require_int(
            payload.get("imputed_cell_count", 0),
            field_name="preprocessing_stage.imputed_cell_count",
        ),
        imputed_row_ids=imputed_row_ids,
        notes=optional_str(
            payload.get("notes"),
            field_name="preprocessing_stage.notes",
        ),
        diagnostics=diagnostics,
        batch_correction_provenance=batch_correction_provenance,
    )


def _require_current_stage_schema_version(payload: Mapping[str, object]) -> int:
    if "schema_version" not in payload:
        raise_legacy_provenance_schema()
    schema_version = require_int(
        payload.get("schema_version"),
        field_name="preprocessing_stage.schema_version",
    )
    if schema_version != PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3:
        raise_legacy_provenance_schema()
    return schema_version


def _reproducibility_caveat_to_payload(
    caveat: ReproducibilityCaveat,
) -> dict[str, object]:
    return caveat.to_payload()


def _reproducibility_caveat_from_payload(
    payload: Mapping[str, object],
) -> ReproducibilityCaveat:
    details = require_mapping(
        payload.get("details", {}),
        field_name="reproducibility_caveat.details",
    )
    return ReproducibilityCaveat(
        code=require_str(
            payload.get("code"),
            field_name="reproducibility_caveat.code",
        ),
        severity=require_str(
            payload.get("severity"),
            field_name="reproducibility_caveat.severity",
        ),
        message=require_str(
            payload.get("message"),
            field_name="reproducibility_caveat.message",
        ),
        details={key: to_json_value(value) for key, value in details.items()},
    )


def _resolve_determinism(value: object, *, field_name: str) -> DeterminismKind:
    if isinstance(value, DeterminismKind):
        return value
    normalized = require_str(value, field_name=field_name)
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
