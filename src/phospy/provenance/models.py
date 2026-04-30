"""Typed machine-readable provenance models."""

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.scientific_policies import ScientificPolicyRecord

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V1 = 1
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2 = 2

JsonPrimitive = str | int | float | bool | None
JsonValue = (
    JsonPrimitive | list["JsonValue"] | tuple["JsonValue", ...] | dict[str, "JsonValue"]
)


def _empty_platform_provenance() -> dict[str, str]:
    return {}


@dataclass(frozen=True, slots=True)
class TableFingerprint:
    """Deterministic table fingerprint metadata."""

    name: str
    rows: int
    columns: int
    index_name: str | None
    column_names: tuple[str, ...]
    dtypes: tuple[str, ...]
    hash_algorithm: str
    hash_value: str


@dataclass(frozen=True, slots=True)
class EnvironmentProvenance:
    """Runtime environment fingerprint for reproducibility."""

    package_name: str
    package_version: str
    python_version: str
    dependency_versions: dict[str, str | None]
    platform: dict[str, str] = field(default_factory=_empty_platform_provenance)


@dataclass(frozen=True, slots=True)
class PreprocessingStageProvenance:
    """Executed preprocessing-stage provenance record."""

    stage: str
    operation: str
    parameters: dict[str, object]
    input_shape: tuple[int, int]
    output_shape: tuple[int, int]
    input_hash: str
    output_hash: str
    dropped_row_ids: tuple[str, ...]
    dropped_row_count: int
    schema_version: int = PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V1
    consumed_input_tables: tuple[TableFingerprint, ...] = ()
    produced_output_tables: tuple[TableFingerprint, ...] = ()
    backend: str | None = None
    random_seed: int | None = None
    is_deterministic: bool = True
    imputed_cell_count: int = 0
    imputed_row_ids: tuple[str, ...] = ()
    notes: str | None = None
    diagnostics: dict[str, object] | None = None


@dataclass(frozen=True, slots=True)
class ReferenceProvenance:
    """Resolved reference identity and table fingerprints."""

    source_type: str
    organism: str
    bundle_id: str | None
    table_fingerprints: tuple[TableFingerprint, ...]


@dataclass(frozen=True, slots=True)
class RunProvenance:
    """Machine-readable run provenance payload."""

    environment: EnvironmentProvenance
    input_tables: tuple[TableFingerprint, ...]
    preprocessing_stages: tuple[PreprocessingStageProvenance, ...]
    reference: ReferenceProvenance | None
    workflow_name: str | None
    workflow_parameters: dict[str, object]
    random_state: int | None
    random_seed_policy: str | None
    output_tables: tuple[TableFingerprint, ...]
    scientific_policies: tuple[ScientificPolicyRecord, ...] = ()


__all__ = [
    "EnvironmentProvenance",
    "JsonValue",
    "PreprocessingStageProvenance",
    "ReferenceProvenance",
    "RunProvenance",
    "TableFingerprint",
]
