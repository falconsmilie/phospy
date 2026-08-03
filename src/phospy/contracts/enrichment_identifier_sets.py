"""Typed enrichment identifier-set provenance contracts."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from typing import TypeVar, cast

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import ContractValidationError
from phospy.provenance.models import InputIntensityScaleEvidence, TableFingerprint
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)

_EnumT = TypeVar("_EnumT", bound=Enum)


class EnrichmentIdentifierSetSourceType(str, Enum):
    """Source category for an enrichment identifier set."""

    MANUAL = "manual"
    RAW_IDENTIFIER_LIST = "raw_identifier_list"
    PHOSPY_DERIVED_QUANTITATIVE = "phospy_derived_quantitative"


class EnrichmentDerivedSetSourceResultKind(str, Enum):
    """Kind of PhosPy quantitative result axis used to derive a set."""

    CONTRAST = "contrast"
    PROFILE = "profile"


class EnrichmentDerivedSetThresholdDirection(str, Enum):
    """Threshold direction used before ORA consumes the derived identifiers."""

    GREATER_THAN = "greater_than"
    GREATER_THAN_OR_EQUAL = "greater_than_or_equal"
    LESS_THAN = "less_than"
    LESS_THAN_OR_EQUAL = "less_than_or_equal"
    ABSOLUTE_GREATER_THAN = "absolute_greater_than"
    ABSOLUTE_GREATER_THAN_OR_EQUAL = "absolute_greater_than_or_equal"


class EnrichmentDerivedSetMissingValueRule(str, Enum):
    """How missing source values were handled while deriving a set."""

    DROP_MISSING = "drop_missing"
    TREAT_MISSING_AS_NOT_SELECTED = "treat_missing_as_not_selected"
    ERROR_ON_MISSING = "error_on_missing"


class EnrichmentDerivedSetValueScale(str, Enum):
    """Scale of the thresholded source value."""

    LINEAR = "linear"
    LOG2 = "log2"
    NEGATIVE_LOG10 = "negative_log10"
    PROBABILITY = "probability"
    Z_SCORE = "z_score"
    UNITLESS = "unitless"


class EnrichmentDerivedSetValueMeaning(str, Enum):
    """Scientific meaning of the thresholded source value."""

    CONTRAST_LOG2_FOLD_CHANGE = "contrast_log2_fold_change"
    DIFFERENTIAL_EFFECT_SIZE = "differential_effect_size"
    P_VALUE = "p_value"
    ADJUSTED_P_VALUE = "adjusted_p_value"
    ACTIVITY_SCORE = "activity_score"
    PROFILE_SCORE = "profile_score"
    MEMBERSHIP_SCORE = "membership_score"


@dataclass(frozen=True, slots=True)
class EnrichmentDerivedQuantitativeSetProvenance:
    """Typed derivation provenance for PhosPy-derived ORA identifier sets."""

    source_result_fingerprint: TableFingerprint
    source_result_kind: EnrichmentDerivedSetSourceResultKind
    source_profile_or_contrast: str
    identifier_namespace: str
    threshold: float
    direction: EnrichmentDerivedSetThresholdDirection
    missing_value_rule: EnrichmentDerivedSetMissingValueRule
    quantitative_scale: EnrichmentDerivedSetValueScale
    quantitative_meaning: EnrichmentDerivedSetValueMeaning
    software_version: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "source_result_fingerprint",
            _require_table_fingerprint(
                self.source_result_fingerprint,
                field_name=(
                    "enrichment derived quantitative set provenance "
                    "source_result_fingerprint"
                ),
            ),
        )
        object.__setattr__(
            self,
            "source_result_kind",
            _coerce_enum(
                self.source_result_kind,
                EnrichmentDerivedSetSourceResultKind,
                field_name=(
                    "enrichment derived quantitative set provenance source_result_kind"
                ),
            ),
        )
        object.__setattr__(
            self,
            "source_profile_or_contrast",
            _require_non_empty_text(
                self.source_profile_or_contrast,
                field_name=(
                    "enrichment derived quantitative set provenance "
                    "source_profile_or_contrast"
                ),
            ),
        )
        object.__setattr__(
            self,
            "identifier_namespace",
            _require_non_empty_text(
                self.identifier_namespace,
                field_name=(
                    "enrichment derived quantitative set provenance "
                    "identifier_namespace"
                ),
            ),
        )
        object.__setattr__(
            self,
            "threshold",
            _require_finite_float(
                self.threshold,
                field_name=("enrichment derived quantitative set provenance threshold"),
            ),
        )
        object.__setattr__(
            self,
            "direction",
            _coerce_enum(
                self.direction,
                EnrichmentDerivedSetThresholdDirection,
                field_name=("enrichment derived quantitative set provenance direction"),
            ),
        )
        object.__setattr__(
            self,
            "missing_value_rule",
            _coerce_enum(
                self.missing_value_rule,
                EnrichmentDerivedSetMissingValueRule,
                field_name=(
                    "enrichment derived quantitative set provenance missing_value_rule"
                ),
            ),
        )
        object.__setattr__(
            self,
            "quantitative_scale",
            _coerce_enum(
                self.quantitative_scale,
                EnrichmentDerivedSetValueScale,
                field_name=(
                    "enrichment derived quantitative set provenance quantitative_scale"
                ),
            ),
        )
        object.__setattr__(
            self,
            "quantitative_meaning",
            _coerce_enum(
                self.quantitative_meaning,
                EnrichmentDerivedSetValueMeaning,
                field_name=(
                    "enrichment derived quantitative set provenance "
                    "quantitative_meaning"
                ),
            ),
        )
        object.__setattr__(
            self,
            "software_version",
            _require_non_empty_text(
                self.software_version,
                field_name=(
                    "enrichment derived quantitative set provenance software_version"
                ),
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a compact JSON-compatible derivation payload."""

        return {
            "source_result_fingerprint": table_fingerprint_to_payload(
                self.source_result_fingerprint
            ),
            "source_result_kind": self.source_result_kind.value,
            "source_profile_or_contrast": self.source_profile_or_contrast,
            "identifier_namespace": self.identifier_namespace,
            "threshold": float(self.threshold),
            "direction": self.direction.value,
            "missing_value_rule": self.missing_value_rule.value,
            "quantitative_scale": self.quantitative_scale.value,
            "quantitative_meaning": self.quantitative_meaning.value,
            "software_version": self.software_version,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> EnrichmentDerivedQuantitativeSetProvenance:
        """Restore typed derivation provenance from its JSON payload."""

        mapping = _require_mapping_payload(
            payload,
            field_name="enrichment derived quantitative set provenance",
        )
        _reject_unexpected_payload_keys(
            mapping,
            expected_keys=frozenset(
                {
                    "source_result_fingerprint",
                    "source_result_kind",
                    "source_profile_or_contrast",
                    "identifier_namespace",
                    "threshold",
                    "direction",
                    "missing_value_rule",
                    "quantitative_scale",
                    "quantitative_meaning",
                    "software_version",
                }
            ),
            field_name="enrichment derived quantitative set provenance",
        )
        return cls(
            source_result_fingerprint=_table_fingerprint_from_payload(
                _require_mapping_payload(
                    mapping.get("source_result_fingerprint"),
                    field_name=(
                        "enrichment derived quantitative set provenance "
                        "source_result_fingerprint"
                    ),
                )
            ),
            source_result_kind=mapping.get("source_result_kind"),  # type: ignore[arg-type]
            source_profile_or_contrast=mapping.get(  # type: ignore[arg-type]
                "source_profile_or_contrast"
            ),
            identifier_namespace=mapping.get("identifier_namespace"),  # type: ignore[arg-type]
            threshold=mapping.get("threshold"),  # type: ignore[arg-type]
            direction=mapping.get("direction"),  # type: ignore[arg-type]
            missing_value_rule=mapping.get("missing_value_rule"),  # type: ignore[arg-type]
            quantitative_scale=mapping.get("quantitative_scale"),  # type: ignore[arg-type]
            quantitative_meaning=mapping.get("quantitative_meaning"),  # type: ignore[arg-type]
            software_version=mapping.get("software_version"),  # type: ignore[arg-type]
        )


@dataclass(frozen=True, slots=True)
class EnrichmentIdentifierSetProvenance:
    """Typed provenance for selected or background enrichment identifiers."""

    source_type: EnrichmentIdentifierSetSourceType
    source_label: str
    identifier_count: int
    upstream_workflow_id: str | None = None
    upstream_result_id: str | None = None
    input_intensity_scale_evidence: InputIntensityScaleEvidence | None = None
    derived_quantitative_provenance: (
        EnrichmentDerivedQuantitativeSetProvenance | None
    ) = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "source_type",
            _coerce_enrichment_identifier_set_source_type(self.source_type),
        )
        object.__setattr__(
            self,
            "source_label",
            _require_non_empty_text(
                self.source_label,
                field_name="enrichment identifier-set provenance source_label",
            ),
        )
        object.__setattr__(
            self,
            "identifier_count",
            _require_non_negative_int(
                self.identifier_count,
                field_name="enrichment identifier-set provenance identifier_count",
            ),
        )
        object.__setattr__(
            self,
            "upstream_workflow_id",
            _normalise_optional_text(
                self.upstream_workflow_id,
                field_name=(
                    "enrichment identifier-set provenance upstream_workflow_id"
                ),
            ),
        )
        object.__setattr__(
            self,
            "upstream_result_id",
            _normalise_optional_text(
                self.upstream_result_id,
                field_name="enrichment identifier-set provenance upstream_result_id",
            ),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_evidence",
            _coerce_input_intensity_scale_evidence(
                self.input_intensity_scale_evidence,
            ),
        )
        object.__setattr__(
            self,
            "derived_quantitative_provenance",
            _coerce_derived_quantitative_provenance(
                self.derived_quantitative_provenance,
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a compact JSON-compatible provenance payload."""

        payload: dict[str, object] = {
            "source_type": self.source_type.value,
            "source_label": self.source_label,
            "identifier_count": int(self.identifier_count),
        }
        if self.upstream_workflow_id is not None:
            payload["upstream_workflow_id"] = self.upstream_workflow_id
        if self.upstream_result_id is not None:
            payload["upstream_result_id"] = self.upstream_result_id
        if self.input_intensity_scale_evidence is not None:
            payload["input_intensity_scale_evidence"] = (
                self.input_intensity_scale_evidence.to_payload()
            )
        if self.derived_quantitative_provenance is not None:
            payload["derived_quantitative_provenance"] = (
                self.derived_quantitative_provenance.to_payload()
            )
        return payload

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> EnrichmentIdentifierSetProvenance:
        """Restore typed identifier-set provenance from its JSON payload."""

        mapping = _require_mapping_payload(
            payload,
            field_name="enrichment identifier-set provenance",
        )
        _reject_unexpected_payload_keys(
            mapping,
            expected_keys=frozenset(
                {
                    "source_type",
                    "source_label",
                    "identifier_count",
                    "upstream_workflow_id",
                    "upstream_result_id",
                    "input_intensity_scale_evidence",
                    "derived_quantitative_provenance",
                }
            ),
            field_name="enrichment identifier-set provenance",
        )
        return cls(
            source_type=mapping.get("source_type"),  # type: ignore[arg-type]
            source_label=mapping.get("source_label"),  # type: ignore[arg-type]
            identifier_count=mapping.get("identifier_count"),  # type: ignore[arg-type]
            upstream_workflow_id=mapping.get("upstream_workflow_id"),  # type: ignore[arg-type]
            upstream_result_id=mapping.get("upstream_result_id"),  # type: ignore[arg-type]
            input_intensity_scale_evidence=mapping.get(  # type: ignore[arg-type]
                "input_intensity_scale_evidence"
            ),
            derived_quantitative_provenance=mapping.get(  # type: ignore[arg-type]
                "derived_quantitative_provenance"
            ),
        )


def _coerce_enrichment_identifier_set_source_type(
    value: object,
) -> EnrichmentIdentifierSetSourceType:
    if isinstance(value, EnrichmentIdentifierSetSourceType):
        return value
    try:
        return EnrichmentIdentifierSetSourceType(str(value).strip())
    except ValueError as exc:
        supported = ", ".join(item.value for item in EnrichmentIdentifierSetSourceType)
        raise ContractValidationError(
            "enrichment identifier-set provenance source_type must be one of: "
            + supported
        ) from exc


def _coerce_input_intensity_scale_evidence(
    value: object | None,
) -> InputIntensityScaleEvidence | None:
    if value is None or isinstance(value, InputIntensityScaleEvidence):
        return value
    if isinstance(value, Mapping):
        return InputIntensityScaleEvidence(
            input_intensity_scale=cast(str, value.get("input_intensity_scale")),
            input_intensity_scale_evidence_level=cast(
                str,
                value.get("input_intensity_scale_evidence_level"),
            ),
            input_intensity_scale_source=cast(
                str,
                value.get("input_intensity_scale_source"),
            ),
            input_intensity_scale_source_detail=value.get(
                "input_intensity_scale_source_detail"
            ),
        )
    raise ContractValidationError(
        "input_intensity_scale_evidence must be InputIntensityScaleEvidence, "
        "mapping, or None"
    )


def _coerce_derived_quantitative_provenance(
    value: object | None,
) -> EnrichmentDerivedQuantitativeSetProvenance | None:
    if value is None or isinstance(value, EnrichmentDerivedQuantitativeSetProvenance):
        return value
    if isinstance(value, Mapping):
        return EnrichmentDerivedQuantitativeSetProvenance.from_payload(
            _require_mapping_payload(
                value,
                field_name=(
                    "enrichment identifier-set provenance "
                    "derived_quantitative_provenance"
                ),
            )
        )
    raise ContractValidationError(
        "derived_quantitative_provenance must be "
        "EnrichmentDerivedQuantitativeSetProvenance, mapping, or None"
    )


def _coerce_enum(
    value: object,
    enum_type: type[_EnumT],
    *,
    field_name: str,
) -> _EnumT:
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(str(value).strip())
    except ValueError as exc:
        supported = ", ".join(item.value for item in enum_type)
        raise ContractValidationError(
            f"{field_name} must be one of: {supported}"
        ) from exc


def _require_table_fingerprint(value: object, *, field_name: str) -> TableFingerprint:
    if not isinstance(value, TableFingerprint):
        raise ContractValidationError(f"{field_name} must be TableFingerprint")
    return value


def _table_fingerprint_from_payload(
    payload: Mapping[str, object],
) -> TableFingerprint:
    try:
        return table_fingerprint_from_payload(payload)
    except PhosPyInputError as exc:
        raise ContractValidationError(
            "enrichment derived quantitative set provenance "
            f"source_result_fingerprint is invalid: {exc}"
        ) from exc


def _require_finite_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ContractValidationError(f"{field_name} must be a finite number")
    normalised = float(value)
    if not math.isfinite(normalised):
        raise ContractValidationError(f"{field_name} must be a finite number")
    return normalised


def _require_mapping_payload(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ContractValidationError(f"{field_name} must be a mapping")
    result: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str):
            raise ContractValidationError(
                f"{field_name} keys must be strings; got {type(key).__name__}"
            )
        if key in result:
            raise ContractValidationError(
                f"{field_name} contains duplicate key {key!r}"
            )
        result[key] = item
    return result


def _reject_unexpected_payload_keys(
    payload: Mapping[str, object],
    *,
    expected_keys: frozenset[str],
    field_name: str,
) -> None:
    unexpected = sorted(set(payload) - set(expected_keys))
    if unexpected:
        raise ContractValidationError(
            f"{field_name} contains unsupported keys: "
            + ", ".join(repr(key) for key in unexpected)
        )


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise ContractValidationError(f"{field_name} must be a non-empty string")
    return value.strip()


def _normalise_optional_text(value: object | None, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_non_empty_text(value, field_name=field_name)


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractValidationError(f"{field_name} must be an int")
    if value < 0:
        raise ContractValidationError(
            f"{field_name} must be greater than or equal to 0"
        )
    return value


__all__ = [
    "EnrichmentDerivedQuantitativeSetProvenance",
    "EnrichmentDerivedSetMissingValueRule",
    "EnrichmentDerivedSetSourceResultKind",
    "EnrichmentDerivedSetThresholdDirection",
    "EnrichmentDerivedSetValueMeaning",
    "EnrichmentDerivedSetValueScale",
    "EnrichmentIdentifierSetProvenance",
    "EnrichmentIdentifierSetSourceType",
]
