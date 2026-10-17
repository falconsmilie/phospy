"""Workflow and preprocessing-stage provenance models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import (
    freeze_json_mapping,
    freeze_optional_json_mapping,
    thaw_json_mapping,
)
from phospy.provenance.models._shared import (
    JsonValue,
    _provenance_string_tuple,
    _required_provenance_text,
    _required_shape,
)
from phospy.provenance.models.batch_correction import (
    BatchCorrectionProvenance as _BatchCorrectionProvenance,
)
from phospy.provenance.models.tables import (
    TableFingerprint,
    _required_table_fingerprint_tuple,
)

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V1 = 1

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2 = 2

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3 = 3


class DeterminismKind(str, Enum):
    """Declared execution determinism for preprocessing stage provenance."""

    DETERMINISTIC = "deterministic"
    SEEDED_STOCHASTIC = "seeded_stochastic"
    EXTERNALLY_NONDETERMINISTIC = "externally_nondeterministic"


PREPROCESSING_STAGE_DETERMINISM_DETERMINISTIC = DeterminismKind.DETERMINISTIC.value

PREPROCESSING_STAGE_DETERMINISM_PURE = PREPROCESSING_STAGE_DETERMINISM_DETERMINISTIC

PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC = (
    DeterminismKind.SEEDED_STOCHASTIC.value
)

PREPROCESSING_STAGE_DETERMINISM_EXTERNALLY_NONDETERMINISTIC = (
    DeterminismKind.EXTERNALLY_NONDETERMINISTIC.value
)

PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY = (
    PREPROCESSING_STAGE_DETERMINISM_EXTERNALLY_NONDETERMINISTIC
)

PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE = (
    "preprocessing_external_nondeterminism"
)

_REPRODUCIBILITY_CAVEAT_SEVERITIES = frozenset({"warning", "error"})


@dataclass(frozen=True, slots=True)
class ReproducibilityCaveat:
    """Machine-readable reproducibility caveat for provenance records."""

    code: str
    severity: str
    message: str
    details: Mapping[str, JsonValue] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "code",
            _required_provenance_text(
                self.code, field_name="reproducibility_caveat.code"
            ),
        )
        severity = _required_provenance_text(
            self.severity,
            field_name="reproducibility_caveat.severity",
        )
        if severity not in _REPRODUCIBILITY_CAVEAT_SEVERITIES:
            supported = ", ".join(sorted(_REPRODUCIBILITY_CAVEAT_SEVERITIES))
            raise PhosPyInputError(
                "reproducibility_caveat.severity must be one of: " + supported
            )
        object.__setattr__(self, "severity", severity)
        object.__setattr__(
            self,
            "message",
            _required_provenance_text(
                self.message,
                field_name="reproducibility_caveat.message",
            ),
        )
        object.__setattr__(
            self,
            "details",
            freeze_json_mapping(
                self.details,
                field_name="reproducibility_caveat.details",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible caveat payload."""

        return {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "details": thaw_json_mapping(
                self.details,
                field_name="reproducibility_caveat.details",
            ),
        }


@dataclass(frozen=True, slots=True)
class PreprocessingStageProvenance:
    """Executed preprocessing-stage provenance record."""

    stage: str
    operation: str
    parameters: Mapping[str, object]
    input_shape: tuple[int, int]
    output_shape: tuple[int, int]
    input_hash: str
    output_hash: str
    dropped_row_ids: tuple[str, ...]
    dropped_row_count: int
    phospho_input_hash: str | None = None
    phospho_output_hash: str | None = None
    schema_version: int = PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3
    consumed_input_tables: tuple[TableFingerprint, ...] = ()
    produced_output_tables: tuple[TableFingerprint, ...] = ()
    backend: str | None = None
    random_seed: int | None = None
    determinism: DeterminismKind | str = DeterminismKind.DETERMINISTIC
    reproducibility_caveats: tuple[ReproducibilityCaveat, ...] = ()
    imputed_cell_count: int = 0
    imputed_row_ids: tuple[str, ...] = ()
    notes: str | None = None
    diagnostics: dict[str, JsonValue] | None = None
    batch_correction_provenance: _BatchCorrectionProvenance | None = None
    _frozen_parameters: Mapping[str, JsonValue] = field(
        init=False,
        repr=False,
        compare=False,
    )
    _frozen_diagnostics: Mapping[str, JsonValue] | None = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __getattribute__(self, name: str) -> object:
        if name == "parameters":
            try:
                frozen = object.__getattribute__(self, "_frozen_parameters")
            except AttributeError:
                return object.__getattribute__(self, name)
            return thaw_json_mapping(
                frozen, field_name="preprocessing_stage.parameters"
            )
        if name == "diagnostics":
            try:
                frozen_optional = object.__getattribute__(
                    self,
                    "_frozen_diagnostics",
                )
            except AttributeError:
                return object.__getattribute__(self, name)
            if frozen_optional is None:
                return None
            return thaw_json_mapping(
                frozen_optional,
                field_name="preprocessing_stage.diagnostics",
            )
        return object.__getattribute__(self, name)

    def __post_init__(self) -> None:
        raw_parameters = object.__getattribute__(self, "parameters")
        frozen_parameters = freeze_json_mapping(
            raw_parameters,
            field_name="preprocessing_stage.parameters",
        )
        object.__setattr__(self, "_frozen_parameters", frozen_parameters)
        object.__setattr__(
            self,
            "parameters",
            frozen_parameters,
        )
        object.__setattr__(
            self,
            "input_shape",
            _required_shape(
                self.input_shape, field_name="preprocessing_stage.input_shape"
            ),
        )
        object.__setattr__(
            self,
            "output_shape",
            _required_shape(
                self.output_shape,
                field_name="preprocessing_stage.output_shape",
            ),
        )
        object.__setattr__(
            self,
            "dropped_row_ids",
            _provenance_string_tuple(
                self.dropped_row_ids,
                field_name="preprocessing_stage.dropped_row_ids",
            ),
        )
        object.__setattr__(
            self,
            "consumed_input_tables",
            _required_table_fingerprint_tuple(
                self.consumed_input_tables,
                field_name="preprocessing_stage.consumed_input_tables",
            ),
        )
        object.__setattr__(
            self,
            "produced_output_tables",
            _required_table_fingerprint_tuple(
                self.produced_output_tables,
                field_name="preprocessing_stage.produced_output_tables",
            ),
        )
        object.__setattr__(
            self,
            "reproducibility_caveats",
            _required_reproducibility_caveat_tuple(
                self.reproducibility_caveats,
                field_name="preprocessing_stage.reproducibility_caveats",
            ),
        )
        object.__setattr__(
            self,
            "imputed_row_ids",
            _provenance_string_tuple(
                self.imputed_row_ids,
                field_name="preprocessing_stage.imputed_row_ids",
            ),
        )
        raw_diagnostics = object.__getattribute__(self, "diagnostics")
        frozen_diagnostics = freeze_optional_json_mapping(
            raw_diagnostics,
            field_name="preprocessing_stage.diagnostics",
        )
        object.__setattr__(self, "_frozen_diagnostics", frozen_diagnostics)
        object.__setattr__(self, "diagnostics", frozen_diagnostics)


def _required_reproducibility_caveat_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[ReproducibilityCaveat, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    caveats = tuple(values)
    for caveat in caveats:
        if not isinstance(caveat, ReproducibilityCaveat):
            raise PhosPyInputError(
                f"{field_name} must contain only ReproducibilityCaveat values"
            )
    return caveats


def _required_preprocessing_stage_tuple(
    values: object,
) -> tuple[PreprocessingStageProvenance, ...]:
    field_name = "run_provenance.preprocessing_stages"
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    stages = tuple(values)
    for stage in stages:
        if not isinstance(stage, PreprocessingStageProvenance):
            raise PhosPyInputError(
                f"{field_name} must contain only PreprocessingStageProvenance values"
            )
    return stages
