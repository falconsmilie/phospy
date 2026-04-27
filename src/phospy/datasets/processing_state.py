"""Dataset preprocessing-state summary models."""

from __future__ import annotations

import math
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from typing import TypeAlias

from phospy.errors.input import PhosPyInputError
from phospy.transformations.models import IntensityScaleState

JsonPrimitive: TypeAlias = None | str | bool | int | float
JsonValue: TypeAlias = JsonPrimitive | list["JsonValue"] | dict[str, "JsonValue"]

TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1 = 1
_V1_KNOWN_DIAGNOSTICS_FIELDS = frozenset(
    (
        "diagnostics_schema_version",
        "policy",
        "requested_policy",
        "resolved_policy",
        "formula",
        "requires_log_scale",
        "input_scale",
        "output_scale",
        "quantitative_meaning",
        "output_quantity",
        "matched_rows",
        "total_table_hash",
        "input_phospho_hash",
        "output_phospho_hash",
    )
)


class TotalProteinCorrectionDiagnostics(Mapping[str, JsonValue]):
    """Typed diagnostics contract for total-protein correction state."""

    diagnostics_schema_version: int

    @classmethod
    def from_payload(
        cls,
        payload: object,
        *,
        field_name: str,
    ) -> TotalProteinCorrectionDiagnostics:
        mapping = _require_mapping(payload, field_name=field_name)
        return TotalProteinCorrectionDiagnosticsV1.from_mapping(
            mapping,
            field_name=field_name,
        )

    def to_payload(self) -> dict[str, JsonValue]:
        """Return normalized diagnostics payload suitable for bundle JSON."""

        raise NotImplementedError


@dataclass(frozen=True, slots=True, eq=False)
class TotalProteinCorrectionDiagnosticsV1(TotalProteinCorrectionDiagnostics):
    """Versioned total-protein correction diagnostics payload (schema v1)."""

    policy: str | None = None
    requested_policy: str | None = None
    resolved_policy: str | None = None
    formula: str | None = None
    requires_log_scale: bool | None = None
    input_scale: str | None = None
    output_scale: str | None = None
    quantitative_meaning: str | None = None
    output_quantity: str | None = None
    matched_rows: int | None = None
    total_table_hash: str | None = None
    input_phospho_hash: str | None = None
    output_phospho_hash: str | None = None
    extra: Mapping[str, JsonValue] = field(default_factory=dict)
    diagnostics_schema_version: int = (
        TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
    )
    _payload: dict[str, JsonValue] = field(init=False, repr=False, compare=False)

    @classmethod
    def from_mapping(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> TotalProteinCorrectionDiagnosticsV1:
        if "diagnostics_schema_version" in payload:
            return cls._from_versioned_payload(payload, field_name=field_name)
        return cls._from_legacy_payload(payload, field_name=field_name)

    @classmethod
    def _from_versioned_payload(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> TotalProteinCorrectionDiagnosticsV1:
        version = _require_int(
            payload.get("diagnostics_schema_version"),
            field_name=f"{field_name}.diagnostics_schema_version",
        )
        if version != TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1:
            raise PhosPyInputError(
                f"{field_name}.diagnostics_schema_version={version!r} is unsupported; "
                f"expected {TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        extra = {
            str(key): value
            for key, value in payload.items()
            if str(key) not in _V1_KNOWN_DIAGNOSTICS_FIELDS
        }
        return cls(
            policy=_require_optional_str(
                payload.get("policy"),
                field_name=f"{field_name}.policy",
            ),
            requested_policy=_require_optional_str(
                payload.get("requested_policy"),
                field_name=f"{field_name}.requested_policy",
            ),
            resolved_policy=_require_optional_str(
                payload.get("resolved_policy"),
                field_name=f"{field_name}.resolved_policy",
            ),
            formula=_require_optional_str(
                payload.get("formula"),
                field_name=f"{field_name}.formula",
            ),
            requires_log_scale=_require_optional_bool(
                payload.get("requires_log_scale"),
                field_name=f"{field_name}.requires_log_scale",
            ),
            input_scale=_require_optional_str(
                payload.get("input_scale"),
                field_name=f"{field_name}.input_scale",
            ),
            output_scale=_require_optional_str(
                payload.get("output_scale"),
                field_name=f"{field_name}.output_scale",
            ),
            quantitative_meaning=_require_optional_str(
                payload.get("quantitative_meaning"),
                field_name=f"{field_name}.quantitative_meaning",
            ),
            output_quantity=_require_optional_str(
                payload.get("output_quantity"),
                field_name=f"{field_name}.output_quantity",
            ),
            matched_rows=_require_optional_non_negative_int(
                payload.get("matched_rows"),
                field_name=f"{field_name}.matched_rows",
            ),
            total_table_hash=_require_optional_str(
                payload.get("total_table_hash"),
                field_name=f"{field_name}.total_table_hash",
            ),
            input_phospho_hash=_require_optional_str(
                payload.get("input_phospho_hash"),
                field_name=f"{field_name}.input_phospho_hash",
            ),
            output_phospho_hash=_require_optional_str(
                payload.get("output_phospho_hash"),
                field_name=f"{field_name}.output_phospho_hash",
            ),
            extra=_validate_json_safe_mapping(
                extra,
                field_name=f"{field_name}.extra",
            ),
            diagnostics_schema_version=version,
        )

    @classmethod
    def _from_legacy_payload(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str,
    ) -> TotalProteinCorrectionDiagnosticsV1:
        normalized = _validate_json_safe_mapping(payload, field_name=field_name)
        known = {
            key: normalized.get(key)
            for key in _V1_KNOWN_DIAGNOSTICS_FIELDS
            if key in normalized
        }
        extra = {
            key: value
            for key, value in normalized.items()
            if key not in _V1_KNOWN_DIAGNOSTICS_FIELDS
        }
        return cls(
            policy=_require_optional_str(
                known.get("policy"),
                field_name=f"{field_name}.policy",
            ),
            requested_policy=_require_optional_str(
                known.get("requested_policy"),
                field_name=f"{field_name}.requested_policy",
            ),
            resolved_policy=_require_optional_str(
                known.get("resolved_policy"),
                field_name=f"{field_name}.resolved_policy",
            ),
            formula=_require_optional_str(
                known.get("formula"),
                field_name=f"{field_name}.formula",
            ),
            requires_log_scale=_require_optional_bool(
                known.get("requires_log_scale"),
                field_name=f"{field_name}.requires_log_scale",
            ),
            input_scale=_require_optional_str(
                known.get("input_scale"),
                field_name=f"{field_name}.input_scale",
            ),
            output_scale=_require_optional_str(
                known.get("output_scale"),
                field_name=f"{field_name}.output_scale",
            ),
            quantitative_meaning=_require_optional_str(
                known.get("quantitative_meaning"),
                field_name=f"{field_name}.quantitative_meaning",
            ),
            output_quantity=_require_optional_str(
                known.get("output_quantity"),
                field_name=f"{field_name}.output_quantity",
            ),
            matched_rows=_require_optional_non_negative_int(
                known.get("matched_rows"),
                field_name=f"{field_name}.matched_rows",
            ),
            total_table_hash=_require_optional_str(
                known.get("total_table_hash"),
                field_name=f"{field_name}.total_table_hash",
            ),
            input_phospho_hash=_require_optional_str(
                known.get("input_phospho_hash"),
                field_name=f"{field_name}.input_phospho_hash",
            ),
            output_phospho_hash=_require_optional_str(
                known.get("output_phospho_hash"),
                field_name=f"{field_name}.output_phospho_hash",
            ),
            extra=extra,
        )

    def __post_init__(self) -> None:
        if (
            self.diagnostics_schema_version
            != TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
        ):
            raise PhosPyInputError(
                "dataset processing state total_protein_correction diagnostics "
                f"schema version must be {TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1}"
            )
        policy = _require_optional_str(
            self.policy,
            field_name="dataset processing state total_protein_correction.diagnostics.policy",
        )
        requested_policy = _require_optional_str(
            self.requested_policy,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "requested_policy"
            ),
        )
        resolved_policy = _require_optional_str(
            self.resolved_policy,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "resolved_policy"
            ),
        )
        formula = _require_optional_str(
            self.formula,
            field_name="dataset processing state total_protein_correction.diagnostics.formula",
        )
        requires_log_scale = _require_optional_bool(
            self.requires_log_scale,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "requires_log_scale"
            ),
        )
        input_scale = _require_optional_str(
            self.input_scale,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "input_scale"
            ),
        )
        output_scale = _require_optional_str(
            self.output_scale,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "output_scale"
            ),
        )
        quantitative_meaning = _require_optional_str(
            self.quantitative_meaning,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "quantitative_meaning"
            ),
        )
        output_quantity = _require_optional_str(
            self.output_quantity,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "output_quantity"
            ),
        )
        matched_rows = _require_optional_non_negative_int(
            self.matched_rows,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "matched_rows"
            ),
        )
        total_table_hash = _require_optional_str(
            self.total_table_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "total_table_hash"
            ),
        )
        input_phospho_hash = _require_optional_str(
            self.input_phospho_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "input_phospho_hash"
            ),
        )
        output_phospho_hash = _require_optional_str(
            self.output_phospho_hash,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics."
                "output_phospho_hash"
            ),
        )
        normalized_extra = _validate_json_safe_mapping(
            self.extra,
            field_name=(
                "dataset processing state total_protein_correction.diagnostics.extra"
            ),
        )
        collisions = sorted(
            set(normalized_extra).intersection(_V1_KNOWN_DIAGNOSTICS_FIELDS)
        )
        if collisions:
            joined = ", ".join(collisions)
            raise PhosPyInputError(
                "dataset processing state total_protein_correction.diagnostics.extra "
                f"contains reserved key(s): {joined}"
            )
        if quantitative_meaning is None and output_quantity is not None:
            quantitative_meaning = output_quantity

        payload: dict[str, JsonValue] = {
            "diagnostics_schema_version": self.diagnostics_schema_version
        }
        _set_optional_payload_value(payload, "policy", policy)
        _set_optional_payload_value(payload, "requested_policy", requested_policy)
        _set_optional_payload_value(payload, "resolved_policy", resolved_policy)
        _set_optional_payload_value(payload, "formula", formula)
        _set_optional_payload_value(payload, "requires_log_scale", requires_log_scale)
        _set_optional_payload_value(payload, "input_scale", input_scale)
        _set_optional_payload_value(payload, "output_scale", output_scale)
        _set_optional_payload_value(
            payload, "quantitative_meaning", quantitative_meaning
        )
        _set_optional_payload_value(payload, "output_quantity", output_quantity)
        _set_optional_payload_value(payload, "matched_rows", matched_rows)
        _set_optional_payload_value(payload, "total_table_hash", total_table_hash)
        _set_optional_payload_value(payload, "input_phospho_hash", input_phospho_hash)
        _set_optional_payload_value(payload, "output_phospho_hash", output_phospho_hash)
        payload.update(normalized_extra)

        object.__setattr__(self, "policy", policy)
        object.__setattr__(self, "requested_policy", requested_policy)
        object.__setattr__(self, "resolved_policy", resolved_policy)
        object.__setattr__(self, "formula", formula)
        object.__setattr__(self, "requires_log_scale", requires_log_scale)
        object.__setattr__(self, "input_scale", input_scale)
        object.__setattr__(self, "output_scale", output_scale)
        object.__setattr__(self, "quantitative_meaning", quantitative_meaning)
        object.__setattr__(self, "output_quantity", output_quantity)
        object.__setattr__(self, "matched_rows", matched_rows)
        object.__setattr__(self, "total_table_hash", total_table_hash)
        object.__setattr__(self, "input_phospho_hash", input_phospho_hash)
        object.__setattr__(self, "output_phospho_hash", output_phospho_hash)
        object.__setattr__(self, "extra", normalized_extra)
        object.__setattr__(self, "_payload", payload)

    def __getitem__(self, key: str) -> JsonValue:
        return self._payload[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._payload)

    def __len__(self) -> int:
        return len(self._payload)

    def to_payload(self) -> dict[str, JsonValue]:
        return dict(self._payload)


@dataclass(frozen=True, slots=True)
class MissingDataState:
    """Missing-data policy state at the analysis-ready dataset boundary."""

    policy: str
    min_observed_values: int | None
    complete_matrix: bool
    imputed: bool


@dataclass(frozen=True, slots=True)
class NormalisationState:
    """Normalisation policy state at the analysis-ready dataset boundary."""

    policy: str


@dataclass(frozen=True, slots=True)
class TotalProteinCorrectionState:
    """Total-protein correction state at the analysis-ready dataset boundary."""

    policy: str
    applied: bool
    formula: str | None = None
    requires_log_scale: bool | None = False
    input_scale: str | None = None
    output_scale: str | None = None
    quantitative_meaning: str | None = None
    diagnostics: TotalProteinCorrectionDiagnostics | None = None

    def __post_init__(self) -> None:
        if self.diagnostics is None:
            return
        if isinstance(self.diagnostics, TotalProteinCorrectionDiagnostics):
            return
        normalized = TotalProteinCorrectionDiagnostics.from_payload(
            self.diagnostics,
            field_name="dataset processing state total_protein_correction.diagnostics",
        )
        object.__setattr__(self, "diagnostics", normalized)


@dataclass(frozen=True, slots=True)
class SiteMatrixState:
    """Site-matrix construction state at the analysis-ready dataset boundary."""

    policy: str
    constructed: bool
    missing_data_policy: str
    minimum_observed_values: int | None
    duplicate_site_policy: str


@dataclass(frozen=True, slots=True)
class ComparisonState:
    """Comparison-building state at the analysis-ready dataset boundary."""

    policy: str
    sample_group_column: str
    pairs: tuple[tuple[str, str], ...] | None


@dataclass(frozen=True, slots=True)
class DatasetProcessingState:
    """Compact summary of preprocessing state at the analysis-ready boundary."""

    intensity_scale: IntensityScaleState
    missing_data: MissingDataState
    normalisation: NormalisationState
    total_protein_correction: TotalProteinCorrectionState
    site_matrix: SiteMatrixState
    comparisons: ComparisonState


def _require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be an object")
    return value


def _require_optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def _require_optional_bool(value: object, *, field_name: str) -> bool | None:
    if value is None:
        return None
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return value


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an int")
    return value


def _require_optional_non_negative_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    parsed = _require_int(value, field_name=field_name)
    if parsed < 0:
        raise PhosPyInputError(f"{field_name} must be >= 0")
    return parsed


def _validate_json_safe_mapping(
    value: Mapping[str, object] | Mapping[str, JsonValue],
    *,
    field_name: str,
) -> dict[str, JsonValue]:
    normalized: dict[str, JsonValue] = {}
    for key, item in value.items():
        if not isinstance(key, str):
            raise PhosPyInputError(
                f"{field_name} must contain only string keys; got key "
                f"{key!r} ({type(key).__name__})"
            )
        normalized[key] = _validate_json_safe_value(item, path=f"{field_name}.{key}")
    return normalized


def _validate_json_safe_value(value: object, *, path: str) -> JsonValue:
    if value is None:
        return None
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise PhosPyInputError(f"{path} must contain only finite float values")
        return value
    if isinstance(value, list):
        return [
            _validate_json_safe_value(item, path=f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, tuple):
        return [
            _validate_json_safe_value(item, path=f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, Mapping):
        nested: dict[str, JsonValue] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise PhosPyInputError(
                    f"{path} must contain only string keys; got key "
                    f"{key!r} ({type(key).__name__})"
                )
            nested[key] = _validate_json_safe_value(item, path=f"{path}.{key}")
        return nested
    raise PhosPyInputError(
        f"{path} contains unsupported value type "
        f"{type(value).__module__}.{type(value).__name__}; expected JSON-safe "
        "scalars, arrays, or objects"
    )


def _set_optional_payload_value(
    payload: dict[str, JsonValue],
    key: str,
    value: JsonValue,
) -> None:
    if value is not None:
        payload[key] = value
