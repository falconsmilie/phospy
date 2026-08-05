"""Transformation-domain provenance value models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import cast

from phospy.errors.transformations import InvalidTransformationStateError
from phospy.provenance.immutability import (
    freeze_json_mapping_with_error_type,
    thaw_json_mapping,
)
from phospy.provenance.models import TableFingerprint
from phospy.science.transformations._payloads import (
    _canonical_table_fingerprint_payload,
    _canonical_table_fingerprint_payload_tuple,
    _normalize_diagnostic_caveat_codes,
    _normalize_quantitative_meaning_evidence_mode,
    _normalize_quantitative_meaning_operation_id,
    _normalize_required_provenance_text,
    _optional_payload_str,
    _require_payload_int,
    _require_payload_mapping,
    _require_payload_sequence,
    _require_payload_str,
)
from phospy.science.transformations.policy import (
    QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1,
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentSource,
    IntensityScaleEvidenceLevel,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    normalize_intensity_scale_evidence_level,
    normalize_optional_quantitative_meaning,
    normalize_required_quantitative_meaning,
)


def _default_provenance_parameters() -> dict[str, object]:
    return {}


@dataclass(frozen=True, slots=True)
class QuantitativeMeaningTransitionProvenance:
    """Structured provenance for quantitative-meaning establishment/transition."""

    source_quantity: QuantitativeMeaning | str | None
    target_quantity: QuantitativeMeaning | str
    operation_id: str
    producer_id: str
    evidence_mode: QuantitativeMeaningEvidenceMode | str
    parameters: Mapping[str, object] = field(
        default_factory=_default_provenance_parameters
    )
    input_table_fingerprints: tuple[TableFingerprint | Mapping[str, object], ...] = ()
    output_table_fingerprint: TableFingerprint | Mapping[str, object] | None = None
    trace_id: str | None = None
    diagnostic_caveat_codes: tuple[str, ...] = ()
    schema_version: int = QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1

    def __post_init__(self) -> None:
        if (
            int(self.schema_version)
            != QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1
        ):
            raise InvalidTransformationStateError(
                "quantitative meaning provenance schema_version must be "
                f"{QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1}"
            )
        object.__setattr__(
            self,
            "source_quantity",
            normalize_optional_quantitative_meaning(self.source_quantity),
        )
        target_quantity = normalize_required_quantitative_meaning(
            self.target_quantity,
            field_name="quantitative_meaning_provenance.target_quantity",
        )
        object.__setattr__(self, "target_quantity", target_quantity)
        evidence_mode = _normalize_quantitative_meaning_evidence_mode(
            self.evidence_mode
        )
        object.__setattr__(self, "evidence_mode", evidence_mode)
        object.__setattr__(
            self,
            "operation_id",
            _normalize_quantitative_meaning_operation_id(self.operation_id),
        )
        object.__setattr__(
            self,
            "producer_id",
            _normalize_required_provenance_text(
                self.producer_id,
                field_name="quantitative_meaning_provenance.producer_id",
            ),
        )
        object.__setattr__(
            self,
            "parameters",
            freeze_json_mapping_with_error_type(
                self.parameters,
                field_name="quantitative_meaning_provenance.parameters",
                error_type=InvalidTransformationStateError,
            ),
        )
        input_fingerprints = _canonical_table_fingerprint_payload_tuple(
            self.input_table_fingerprints,
            field_name="quantitative_meaning_provenance.input_table_fingerprints",
        )
        object.__setattr__(self, "input_table_fingerprints", input_fingerprints)
        output_fingerprint = (
            None
            if self.output_table_fingerprint is None
            else _canonical_table_fingerprint_payload(
                self.output_table_fingerprint,
                field_name=("quantitative_meaning_provenance.output_table_fingerprint"),
            )
        )
        object.__setattr__(self, "output_table_fingerprint", output_fingerprint)
        trace_id = None if self.trace_id is None else str(self.trace_id).strip() or None
        object.__setattr__(self, "trace_id", trace_id)
        object.__setattr__(
            self,
            "diagnostic_caveat_codes",
            _normalize_diagnostic_caveat_codes(self.diagnostic_caveat_codes),
        )
        if (
            evidence_mode is QuantitativeMeaningEvidenceMode.DERIVED_BY_PHOSPY_OPERATION
            and self.source_quantity is None
        ):
            raise InvalidTransformationStateError(
                "derived quantitative meaning transitions require a source "
                "quantitative meaning"
            )
        if (
            evidence_mode is QuantitativeMeaningEvidenceMode.DERIVED_BY_PHOSPY_OPERATION
            and (not input_fingerprints or output_fingerprint is None)
        ):
            raise InvalidTransformationStateError(
                "derived quantitative meaning transitions require input table "
                "fingerprints and an output table fingerprint"
            )

    def to_payload(self) -> dict[str, object]:
        """Return JSON-safe quantitative-meaning provenance payload."""

        source_quantity = cast(QuantitativeMeaning | None, self.source_quantity)
        target_quantity = cast(QuantitativeMeaning, self.target_quantity)
        evidence_mode = cast(QuantitativeMeaningEvidenceMode, self.evidence_mode)
        return {
            "schema_version": int(self.schema_version),
            "source_quantity": (
                None if source_quantity is None else source_quantity.value
            ),
            "target_quantity": target_quantity.value,
            "operation_id": self.operation_id,
            "producer_id": self.producer_id,
            "evidence_mode": evidence_mode.value,
            "parameters": thaw_json_mapping(
                self.parameters,
                field_name="quantitative_meaning_provenance.parameters",
            ),
            "input_table_fingerprints": [
                thaw_json_mapping(
                    fingerprint,
                    field_name=(
                        "quantitative_meaning_provenance.input_table_fingerprints[]"
                    ),
                )
                for fingerprint in self.input_table_fingerprints
            ],
            "output_table_fingerprint": (
                None
                if self.output_table_fingerprint is None
                else thaw_json_mapping(
                    self.output_table_fingerprint,
                    field_name=(
                        "quantitative_meaning_provenance.output_table_fingerprint"
                    ),
                )
            ),
            "trace_id": self.trace_id,
            "diagnostic_caveat_codes": list(self.diagnostic_caveat_codes),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> QuantitativeMeaningTransitionProvenance:
        """Deserialize quantitative-meaning provenance from a payload."""

        payload = _require_payload_mapping(
            payload,
            field_name="quantitative_meaning_provenance",
        )
        input_fingerprints = tuple(
            _require_payload_mapping(
                item,
                field_name=(
                    "quantitative_meaning_provenance.input_table_fingerprints"
                    f"[{position}]"
                ),
            )
            for position, item in enumerate(
                _require_payload_sequence(
                    payload.get("input_table_fingerprints"),
                    field_name=(
                        "quantitative_meaning_provenance.input_table_fingerprints"
                    ),
                )
            )
        )
        output_raw = payload.get("output_table_fingerprint")
        return cls(
            schema_version=_require_payload_int(
                payload.get("schema_version"),
                field_name="quantitative_meaning_provenance.schema_version",
            ),
            source_quantity=_optional_payload_str(
                payload.get("source_quantity"),
                field_name="quantitative_meaning_provenance.source_quantity",
            ),
            target_quantity=_require_payload_str(
                payload.get("target_quantity"),
                field_name="quantitative_meaning_provenance.target_quantity",
            ),
            operation_id=_require_payload_str(
                payload.get("operation_id"),
                field_name="quantitative_meaning_provenance.operation_id",
            ),
            producer_id=_require_payload_str(
                payload.get("producer_id"),
                field_name="quantitative_meaning_provenance.producer_id",
            ),
            evidence_mode=_require_payload_str(
                payload.get("evidence_mode"),
                field_name="quantitative_meaning_provenance.evidence_mode",
            ),
            parameters=_require_payload_mapping(
                payload.get("parameters"),
                field_name="quantitative_meaning_provenance.parameters",
            ),
            input_table_fingerprints=input_fingerprints,
            output_table_fingerprint=(
                None
                if output_raw is None
                else _require_payload_mapping(
                    output_raw,
                    field_name=(
                        "quantitative_meaning_provenance.output_table_fingerprint"
                    ),
                )
            ),
            trace_id=_optional_payload_str(
                payload.get("trace_id"),
                field_name="quantitative_meaning_provenance.trace_id",
            ),
            diagnostic_caveat_codes=tuple(
                _require_payload_str(
                    item,
                    field_name=(
                        "quantitative_meaning_provenance.diagnostic_caveat_codes"
                        f"[{position}]"
                    ),
                )
                for position, item in enumerate(
                    _require_payload_sequence(
                        payload.get("diagnostic_caveat_codes"),
                        field_name=(
                            "quantitative_meaning_provenance.diagnostic_caveat_codes"
                        ),
                    )
                )
            ),
        )


@dataclass(frozen=True, slots=True)
class IntensityScaleEstablishmentProvenance:
    """Structured provenance for intensity-scale establishment."""

    scale: str
    mode: IntensityScaleEstablishmentMode
    source: IntensityScaleEstablishmentSource
    evidence_level: IntensityScaleEvidenceLevel = IntensityScaleEvidenceLevel.UNKNOWN
    transformer_name: str | None = None
    input_declaration_source: str | None = None
    parameters: Mapping[str, object] = field(
        default_factory=_default_provenance_parameters
    )
    trace_id: str | None = None
    diagnostic_warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        scale = str(self.scale).strip()
        if not scale:
            raise InvalidTransformationStateError(
                "intensity-scale establishment provenance requires non-empty scale"
            )
        object.__setattr__(self, "scale", scale)
        if not isinstance(cast(object, self.mode), IntensityScaleEstablishmentMode):
            try:
                mode = IntensityScaleEstablishmentMode(str(self.mode))
            except ValueError as exc:
                supported = ", ".join(
                    item.value for item in IntensityScaleEstablishmentMode
                )
                raise InvalidTransformationStateError(
                    "unsupported intensity-scale establishment mode "
                    f"{self.mode!r}; supported: {supported}"
                ) from exc
            object.__setattr__(self, "mode", mode)
        if not isinstance(cast(object, self.source), IntensityScaleEstablishmentSource):
            try:
                source = IntensityScaleEstablishmentSource(str(self.source))
            except ValueError as exc:
                supported = ", ".join(
                    item.value for item in IntensityScaleEstablishmentSource
                )
                raise InvalidTransformationStateError(
                    "unsupported intensity-scale establishment source "
                    f"{self.source!r}; supported: {supported}"
                ) from exc
            object.__setattr__(self, "source", source)
        evidence_level = normalize_intensity_scale_evidence_level(self.evidence_level)
        object.__setattr__(self, "evidence_level", evidence_level)
        transformer_name = (
            None
            if self.transformer_name is None
            else str(self.transformer_name).strip() or None
        )
        object.__setattr__(self, "transformer_name", transformer_name)
        declaration_source = (
            None
            if self.input_declaration_source is None
            else str(self.input_declaration_source).strip() or None
        )
        object.__setattr__(self, "input_declaration_source", declaration_source)
        trace_id = None if self.trace_id is None else str(self.trace_id).strip() or None
        object.__setattr__(self, "trace_id", trace_id)
        object.__setattr__(
            self,
            "parameters",
            freeze_json_mapping_with_error_type(
                self.parameters,
                field_name="intensity_scale_establishment.parameters",
                error_type=InvalidTransformationStateError,
            ),
        )
        object.__setattr__(
            self,
            "diagnostic_warnings",
            tuple(
                str(item).strip()
                for item in self.diagnostic_warnings
                if str(item).strip()
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return JSON-safe payload for provenance surfaces."""

        return {
            "scale": self.scale,
            "establishment_mode": self.mode.value,
            "establishment_source": self.source.value,
            "evidence_level": self.evidence_level.value,
            "transformer_name": self.transformer_name,
            "input_declaration_source": self.input_declaration_source,
            "parameters": thaw_json_mapping(
                self.parameters,
                field_name="intensity_scale_establishment.parameters",
            ),
            "trace_id": self.trace_id,
            "diagnostic_warnings": list(self.diagnostic_warnings),
        }


__all__ = [
    "IntensityScaleEstablishmentProvenance",
    "QuantitativeMeaningTransitionProvenance",
]
