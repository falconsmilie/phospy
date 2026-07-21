"""Transformation-domain intensity scale models."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from typing import Final, cast

from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import InvalidTransformationStateError
from phospy.provenance.immutability import (
    freeze_json_mapping_with_error_type,
    thaw_json_mapping,
)
from phospy.provenance.models import TableFingerprint
from phospy.provenance.serialization import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)
from phospy.science.transformations._authority import (
    EstablishmentAuthority,
    QuantitativeMeaningTransitionAuthority,
    resolve_establishment_authority_source,
    resolve_quantitative_meaning_transition_authority_source,
)

IDENTITY_INTENSITY_SCALE_ESTABLISHER: Final[str] = (
    "phospy.science.transformations.transformers.identity"
)
QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1: Final[int] = 1
QUANTITATIVE_MEANING_OPERATION_CALLER_DECLARATION: Final[str] = (
    "phospy.dataset_builder.quantitative_meaning.declaration"
)
QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE: Final[str] = (
    "phospy.dataset_builder.quantitative_meaning.infer_from_scale_contract"
)
QUANTITATIVE_MEANING_OPERATION_TOTAL_PROTEIN_SUBTRACT_LOG_TOTAL: Final[str] = (
    "phospy.dataset_preprocessing.total_protein_correction.subtract_log_total"
)
QUANTITATIVE_MEANING_OPERATION_LEGACY_BUNDLE_MIGRATION: Final[str] = (
    "phospy.bundle.legacy_intensity_scale_state_quantitative_meaning_migration"
)
QUANTITATIVE_MEANING_USER_DECLARED_CAVEAT_CODE: Final[str] = (
    "quantitative_meaning_user_declared"
)
QUANTITATIVE_MEANING_LEGACY_UNVERIFIED_CAVEAT_CODE: Final[str] = (
    "quantitative_meaning_legacy_unverified"
)
_ESTABLISHED_INTENSITY_SCALE_STATE_MARKER: Final[object] = object()
_QUANTITATIVE_MEANING_OPERATION_PATTERN: Final[re.Pattern[str]] = re.compile(
    r"^[a-z0-9][a-z0-9_]*(?:\.[a-z0-9][a-z0-9_]*)*$"
)
_DATASET_QUANTITATIVE_MEANING_AUTHORITY_SOURCE: Final[str] = (
    "phospy.science.datasets.preprocessing.state_builder"
)
_BUNDLE_QUANTITATIVE_MEANING_AUTHORITY_SOURCE: Final[str] = (
    "phospy.io.bundles._shared.intensity_scale_state"
)


def _default_provenance_parameters() -> dict[str, object]:
    return {}


class IntensityScaleKind(str, Enum):
    """Supported quantitative-intensity scales."""

    LINEAR = "linear"
    LOG2 = "log2"


class IntensityScaleEstablishmentMode(str, Enum):
    """How intensity-scale state was established."""

    DECLARED = "declared"
    TRANSFORMED = "transformed"
    IDENTITY = "identity"
    DERIVED = "derived"


class DeclaredIntensityScaleDiagnosticPolicy(str, Enum):
    """Policy for suspicious declared input intensity-scale diagnostics."""

    WARN = "warn"
    ERROR = "error"


class IntensityScaleEstablishmentSource(str, Enum):
    """Provenance source for how intensity-scale truth was established."""

    TRANSFORMED_BY_PHOSPY = "transformed_by_phospy"
    DECLARED_BY_USER = "declared_by_user"
    RESTORED_FROM_TRUSTED_PROVENANCE = "restored_from_trusted_provenance"


class IntensityScaleEvidenceLevel(str, Enum):
    """Evidence level supporting an intensity-scale transition event."""

    OBSERVED_TRANSFORMATION = "observed_transformation"
    DECLARED_BY_USER = "declared_by_user"
    INFERRED_FROM_METADATA = "inferred_from_metadata"
    UNKNOWN = "unknown"


class QuantitativeMeaning(str, Enum):
    """Scientific interpretation of phospho matrix values."""

    PHOSPHOSITE_ABUNDANCE = "phosphosite_abundance"
    PHOSPHOSITE_LOG_ABUNDANCE = "phosphosite_log_abundance"
    PHOSPHO_TOTAL_LOG_RATIO = "phospho_total_log_ratio"
    CONTRAST_LOG2_FOLD_CHANGE = "contrast_log2_fold_change"
    DIFFERENTIAL_EFFECT_SIZE = "differential_effect_size"
    ACTIVITY_SCORE = "activity_score"
    MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE = (
        "mixed_phospho_total_log_ratio_and_phosphosite_log_abundance"
    )
    UNKNOWN = "unknown"


class QuantitativeMeaningEvidenceMode(str, Enum):
    """Evidence mode supporting quantitative-meaning state."""

    DERIVED_BY_PHOSPY_OPERATION = "derived_by_phospy_operation"
    DECLARED_BY_CALLER = "declared_by_caller"
    RESTORED_FROM_TRUSTED_SERIALIZED_PROVENANCE = (
        "restored_from_trusted_serialized_provenance"
    )
    INFERRED_FROM_SCALE_CONTRACT = "inferred_from_scale_contract"
    LEGACY_UNVERIFIED = "legacy_unverified"


CALLER_DECLARABLE_QUANTITATIVE_MEANINGS: Final[frozenset[QuantitativeMeaning]] = (
    frozenset(
        {
            QuantitativeMeaning.UNKNOWN,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        }
    )
)


@dataclass(frozen=True, slots=True)
class QuantitativeMeaningScaleRule:
    """Allowed scale and semantic role for one quantitative meaning."""

    meaning: QuantitativeMeaning
    allowed_scales: frozenset[IntensityScaleKind]
    semantic_role: str


_QUANTITATIVE_MEANING_SCALE_RULES: Final[tuple[QuantitativeMeaningScaleRule, ...]] = (
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        allowed_scales=frozenset({IntensityScaleKind.LINEAR}),
        semantic_role="phosphosite_abundance_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="phosphosite_abundance_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="total_corrected_log_ratio_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="contrast_or_effect_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="contrast_or_effect_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.ACTIVITY_SCORE,
        allowed_scales=frozenset({IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2}),
        semantic_role="activity_score_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=(
            QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE
        ),
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="mixed_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.UNKNOWN,
        allowed_scales=frozenset({IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2}),
        semantic_role="unknown_matrix",
    ),
)
_QUANTITATIVE_MEANING_SCALE_RULE_BY_MEANING: Final[
    dict[QuantitativeMeaning, QuantitativeMeaningScaleRule]
] = {rule.meaning: rule for rule in _QUANTITATIVE_MEANING_SCALE_RULES}
if set(_QUANTITATIVE_MEANING_SCALE_RULE_BY_MEANING) != set(QuantitativeMeaning):
    raise RuntimeError(
        "QuantitativeMeaning scale rules must cover every QuantitativeMeaning member"
    )


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
            _normalize_optional_quantitative_meaning(self.source_quantity),
        )
        target_quantity = _normalize_required_quantitative_meaning(
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
        evidence_level = _normalize_intensity_scale_evidence_level(self.evidence_level)
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


@dataclass(frozen=True, slots=True)
class MatrixIntensityScaleState:
    """Established intensity scale for one quantitative matrix."""

    kind: IntensityScaleKind
    transformed: bool
    established_by: str = IDENTITY_INTENSITY_SCALE_ESTABLISHER

    def __post_init__(self) -> None:
        if self.kind is IntensityScaleKind.LINEAR and self.transformed:
            raise InvalidTransformationStateError(
                "linear matrix state cannot be marked as transformed"
            )
        if self.kind is IntensityScaleKind.LOG2 and not self.transformed:
            raise InvalidTransformationStateError(
                "log2 matrix state must be marked as transformed"
            )
        if not self.established_by.strip():
            raise InvalidTransformationStateError(
                "matrix intensity scale state requires a non-empty established_by value"
            )

    @classmethod
    def linear(
        cls,
        *,
        established_by: str = IDENTITY_INTENSITY_SCALE_ESTABLISHER,
    ) -> MatrixIntensityScaleState:
        """Construct an established linear (untransformed) matrix state."""

        return cls(
            kind=IntensityScaleKind.LINEAR,
            transformed=False,
            established_by=established_by,
        )

    @classmethod
    def log2(
        cls,
        *,
        established_by: str = "phospy.science.transformations.transformers.log2",
    ) -> MatrixIntensityScaleState:
        """Construct an established log2-transformed matrix state."""

        return cls(
            kind=IntensityScaleKind.LOG2,
            transformed=True,
            established_by=established_by,
        )


@dataclass(frozen=True, slots=True)
class IntensityTransformationEvent:
    """Typed scientific evidence for an intensity-scale transition."""

    transformer_name: str
    input_scale: MatrixIntensityScaleState
    output_scale: MatrixIntensityScaleState
    evidence_level: IntensityScaleEvidenceLevel
    transformation_kind: str
    pseudocount: float | None = None
    input_fingerprint: str | None = None
    output_fingerprint: str | None = None

    def __post_init__(self) -> None:
        transformer_name = _normalize_required_text(
            self.transformer_name,
            field_name="transformer_name",
        )
        object.__setattr__(self, "transformer_name", transformer_name)

        evidence_level = _normalize_intensity_scale_evidence_level(self.evidence_level)
        object.__setattr__(self, "evidence_level", evidence_level)

        if not isinstance(cast(object, self.input_scale), MatrixIntensityScaleState):
            raise InvalidTransformationStateError(
                "intensity transformation event input_scale must be a "
                "MatrixIntensityScaleState"
            )
        if not isinstance(cast(object, self.output_scale), MatrixIntensityScaleState):
            if evidence_level is IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION:
                raise InvalidTransformationStateError(
                    "observed intensity transformation event requires a known "
                    "output scale"
                )
            raise InvalidTransformationStateError(
                "intensity transformation event output_scale must be a "
                "MatrixIntensityScaleState"
            )

        transformation_kind = _normalize_required_text(
            self.transformation_kind,
            field_name="transformation_kind",
        )
        object.__setattr__(self, "transformation_kind", transformation_kind)

        pseudocount = _normalize_intensity_transformation_pseudocount(self.pseudocount)
        object.__setattr__(self, "pseudocount", pseudocount)
        object.__setattr__(
            self,
            "input_fingerprint",
            _normalize_optional_text(self.input_fingerprint),
        )
        object.__setattr__(
            self,
            "output_fingerprint",
            _normalize_optional_text(self.output_fingerprint),
        )

        _validate_intensity_transformation_event_transition(
            input_scale=self.input_scale,
            output_scale=self.output_scale,
            evidence_level=evidence_level,
            transformation_kind=transformation_kind,
        )

    def to_payload(self) -> dict[str, object]:
        """Return JSON-safe payload for reporting/provenance surfaces."""

        return {
            "transformer_name": self.transformer_name,
            "input_scale": _matrix_intensity_scale_payload(self.input_scale),
            "output_scale": _matrix_intensity_scale_payload(self.output_scale),
            "evidence_level": self.evidence_level.value,
            "transformation_kind": self.transformation_kind,
            "pseudocount": self.pseudocount,
            "input_fingerprint": self.input_fingerprint,
            "output_fingerprint": self.output_fingerprint,
        }


@dataclass(frozen=True, slots=True)
class IntensityScaleState:
    """Validated intensity-scale metadata for dataset quantitative matrices."""

    phospho: MatrixIntensityScaleState = field(
        default_factory=MatrixIntensityScaleState.linear
    )
    total: MatrixIntensityScaleState | None = None
    quantity: QuantitativeMeaning | None = None
    _established_via: str | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )
    _establishment_authority_source: str | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )
    _establishment_marker: object | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )
    _establishment_mode: IntensityScaleEstablishmentMode | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )
    _establishment_provenance: IntensityScaleEstablishmentProvenance | None = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )
    _quantitative_meaning_provenance: QuantitativeMeaningTransitionProvenance | None = (
        field(
            default=None,
            init=False,
            repr=False,
            compare=True,
        )
    )

    def __post_init__(self) -> None:
        if not isinstance(cast(object, self.phospho), MatrixIntensityScaleState):
            raise InvalidTransformationStateError(
                "intensity_scale_state.phospho must be a MatrixIntensityScaleState"
            )
        if self.total is not None and not isinstance(
            cast(object, self.total), MatrixIntensityScaleState
        ):
            raise InvalidTransformationStateError(
                "intensity_scale_state.total must be a MatrixIntensityScaleState or None"
            )
        quantity = _normalize_quantitative_meaning(self.quantity)
        if quantity is not None:
            _validate_quantitative_meaning_kind_coherence(
                quantity=quantity,
                kind=self.phospho.kind,
            )
        object.__setattr__(self, "quantity", quantity)

    @property
    def kind(self) -> IntensityScaleKind:
        """Primary intensity scale for phosphosite intensities."""

        return self.phospho.kind

    @property
    def label(self) -> str:
        """Derived human-readable intensity-scale label."""

        if self.total is None or self.total.kind is self.phospho.kind:
            return self.phospho.kind.value
        return "mixed"

    @property
    def established_via(self) -> str | None:
        """Supported establishment source for this state, when available."""

        if self.is_established:
            return self._established_via
        return None

    @property
    def is_established(self) -> bool:
        """Whether this state was established through a supported PhosPy path."""

        return (
            self._establishment_marker is _ESTABLISHED_INTENSITY_SCALE_STATE_MARKER
            and self._established_via is not None
            and self._establishment_authority_source is not None
            and self._establishment_mode is not None
            and self._establishment_provenance is not None
        )

    @property
    def establishment_mode(self) -> IntensityScaleEstablishmentMode | None:
        """Intensity-scale establishment mode for established states."""

        if not self.is_established:
            return None
        return self._establishment_mode

    @property
    def establishment_provenance(self) -> IntensityScaleEstablishmentProvenance | None:
        """Structured establishment provenance, when state is established."""

        if not self.is_established:
            return None
        return self._establishment_provenance

    @property
    def establishment_authority_source(self) -> str | None:
        """Internal authority lane that established this state, when available."""

        if not self.is_established:
            return None
        return self._establishment_authority_source

    @property
    def quantitative_meaning_provenance(
        self,
    ) -> QuantitativeMeaningTransitionProvenance | None:
        """Structured quantitative-meaning provenance, when established."""

        return self._quantitative_meaning_provenance

    @classmethod
    def raw(cls, *, has_total_matrix: bool = False) -> IntensityScaleState:
        """Create a declared linear state (not yet established)."""

        if has_total_matrix:
            return cls(
                phospho=MatrixIntensityScaleState.linear(),
                total=MatrixIntensityScaleState.linear(),
            )
        return cls(phospho=MatrixIntensityScaleState.linear(), total=None)

    @classmethod
    def established_raw(
        cls,
        *,
        has_total_matrix: bool = False,
        established_via: str = IDENTITY_INTENSITY_SCALE_ESTABLISHER,
        establishment_mode: IntensityScaleEstablishmentMode = (
            IntensityScaleEstablishmentMode.IDENTITY
        ),
        evidence_level: IntensityScaleEvidenceLevel = (
            IntensityScaleEvidenceLevel.UNKNOWN
        ),
        transformer_name: str | None = IDENTITY_INTENSITY_SCALE_ESTABLISHER,
        input_declaration_source: str | None = None,
        parameters: Mapping[str, object] | None = None,
        trace_id: str | None = None,
        diagnostic_warnings: tuple[str, ...] = (),
        _authority: EstablishmentAuthority | None = None,
    ) -> IntensityScaleState:
        """Create a linear state through an approved internal authority."""

        return cls.raw(has_total_matrix=has_total_matrix).with_establishment(
            established_via=established_via,
            authority=_authority,
            establishment_mode=establishment_mode,
            evidence_level=evidence_level,
            transformer_name=transformer_name,
            input_declaration_source=input_declaration_source,
            parameters={} if parameters is None else parameters,
            trace_id=trace_id,
            diagnostic_warnings=diagnostic_warnings,
        )

    def with_establishment(
        self,
        *,
        established_via: str,
        authority: EstablishmentAuthority | None,
        establishment_mode: IntensityScaleEstablishmentMode,
        evidence_level: IntensityScaleEvidenceLevel | str,
        transformer_name: str | None,
        input_declaration_source: str | None,
        parameters: Mapping[str, object],
        trace_id: str | None,
        diagnostic_warnings: tuple[str, ...],
    ) -> IntensityScaleState:
        source = established_via.strip()
        if not source:
            raise InvalidTransformationStateError(
                "intensity scale state establishment source must be a non-empty string"
            )
        if isinstance(
            cast(object, establishment_mode), IntensityScaleEstablishmentMode
        ):
            resolved_mode = establishment_mode
        else:
            try:
                resolved_mode = IntensityScaleEstablishmentMode(str(establishment_mode))
            except ValueError as exc:
                supported = ", ".join(
                    item.value for item in IntensityScaleEstablishmentMode
                )
                raise InvalidTransformationStateError(
                    "unsupported intensity-scale establishment mode "
                    f"{establishment_mode!r}; supported: {supported}"
                ) from exc
        authority_source = resolve_establishment_authority_source(authority)
        resolved_evidence_level = _normalize_intensity_scale_evidence_level(
            evidence_level
        )
        established = IntensityScaleState(
            phospho=self.phospho,
            total=self.total,
            quantity=self.quantity,
        )
        object.__setattr__(established, "_established_via", source)
        object.__setattr__(
            established,
            "_establishment_authority_source",
            authority_source,
        )
        object.__setattr__(
            established,
            "_establishment_marker",
            _ESTABLISHED_INTENSITY_SCALE_STATE_MARKER,
        )
        object.__setattr__(established, "_establishment_mode", resolved_mode)
        object.__setattr__(
            established,
            "_establishment_provenance",
            IntensityScaleEstablishmentProvenance(
                scale=established.label,
                mode=resolved_mode,
                source=_resolve_establishment_source(
                    authority_source=authority_source,
                    establishment_mode=resolved_mode,
                ),
                evidence_level=resolved_evidence_level,
                transformer_name=transformer_name,
                input_declaration_source=input_declaration_source,
                parameters=parameters,
                trace_id=trace_id,
                diagnostic_warnings=diagnostic_warnings,
            ),
        )
        object.__setattr__(
            established,
            "_quantitative_meaning_provenance",
            self._quantitative_meaning_provenance,
        )
        return established

    def with_quantitative_meaning(
        self,
        quantity: QuantitativeMeaning,
    ) -> IntensityScaleState:
        """Deprecated unrestricted relabel path; use authority-gated transition."""

        raise InvalidTransformationStateError(
            "IntensityScaleState.with_quantitative_meaning() is no longer supported "
            "because changing quantitative meaning requires explicit semantic "
            "provenance and transition authority. Use "
            "transition_quantitative_meaning(..., provenance=..., authority=...)."
        )

    def transition_quantitative_meaning(
        self,
        *,
        target_quantity: QuantitativeMeaning | str,
        provenance: QuantitativeMeaningTransitionProvenance,
        authority: QuantitativeMeaningTransitionAuthority | None,
    ) -> IntensityScaleState:
        """Return state with authority-gated quantitative-meaning provenance."""

        authority_source = resolve_quantitative_meaning_transition_authority_source(
            authority
        )
        if authority_source != _DATASET_QUANTITATIVE_MEANING_AUTHORITY_SOURCE:
            raise InvalidTransformationStateError(
                "bundle restoration authority cannot mint a new quantitative "
                "meaning transition"
            )
        target = _normalize_required_quantitative_meaning(
            target_quantity,
            field_name="quantitative meaning transition target_quantity",
        )
        if not isinstance(
            cast(object, provenance), QuantitativeMeaningTransitionProvenance
        ):
            raise InvalidTransformationStateError(
                "quantitative meaning transition requires "
                "QuantitativeMeaningTransitionProvenance"
            )
        if not self.is_established:
            raise InvalidTransformationStateError(
                "quantitative meaning transition requires an established intensity "
                "scale state"
            )
        _validate_quantitative_meaning_transition_contract(
            current_quantity=self.quantity,
            target_quantity=target,
            provenance=provenance,
            kind=self.phospho.kind,
        )
        if target is self.quantity:
            if self._quantitative_meaning_provenance == provenance:
                return self
            raise InvalidTransformationStateError(
                "quantitative meaning transition target equals the current meaning "
                "but supplied provenance differs from the existing semantic "
                "provenance"
            )
        provenance_source = cast(QuantitativeMeaning | None, provenance.source_quantity)
        if provenance_source != self.quantity:
            current = None if self.quantity is None else self.quantity.value
            source = None if provenance_source is None else provenance_source.value
            raise InvalidTransformationStateError(
                "quantitative meaning transition provenance source_quantity must "
                f"match current state quantity; current={current!r}, "
                f"provenance source={source!r}"
            )
        updated = IntensityScaleState(
            phospho=self.phospho,
            total=self.total,
            quantity=target,
        )
        _copy_intensity_scale_establishment(source=self, target=updated)
        object.__setattr__(
            updated,
            "_quantitative_meaning_provenance",
            provenance,
        )
        return updated

    def restore_quantitative_meaning_provenance(
        self,
        *,
        provenance: QuantitativeMeaningTransitionProvenance,
        authority: QuantitativeMeaningTransitionAuthority | None,
    ) -> IntensityScaleState:
        """Restore validated serialized quantitative-meaning provenance."""

        authority_source = resolve_quantitative_meaning_transition_authority_source(
            authority
        )
        if authority_source != _BUNDLE_QUANTITATIVE_MEANING_AUTHORITY_SOURCE:
            raise InvalidTransformationStateError(
                "quantitative meaning provenance restoration requires bundle "
                "reconstruction authority"
            )
        if not isinstance(
            cast(object, provenance), QuantitativeMeaningTransitionProvenance
        ):
            raise InvalidTransformationStateError(
                "quantitative meaning provenance restoration requires "
                "QuantitativeMeaningTransitionProvenance"
            )
        if not self.is_established:
            raise InvalidTransformationStateError(
                "quantitative meaning provenance restoration requires an "
                "established intensity scale state"
            )
        source_quantity = cast(QuantitativeMeaning | None, provenance.source_quantity)
        target_quantity = cast(QuantitativeMeaning, provenance.target_quantity)
        _validate_quantitative_meaning_source_target_scale_contract(
            source_quantity=source_quantity,
            target_quantity=target_quantity,
            kind=self.phospho.kind,
        )
        updated = IntensityScaleState(
            phospho=self.phospho,
            total=self.total,
            quantity=target_quantity,
        )
        _copy_intensity_scale_establishment(source=self, target=updated)
        object.__setattr__(
            updated,
            "_quantitative_meaning_provenance",
            provenance,
        )
        return updated


def establish_intensity_scale_state(
    state: IntensityScaleState,
    *,
    established_via: str,
    establishment_mode: IntensityScaleEstablishmentMode = (
        IntensityScaleEstablishmentMode.DERIVED
    ),
    evidence_level: IntensityScaleEvidenceLevel = IntensityScaleEvidenceLevel.UNKNOWN,
    transformer_name: str | None = None,
    input_declaration_source: str | None = None,
    parameters: Mapping[str, object] | None = None,
    trace_id: str | None = None,
    diagnostic_warnings: tuple[str, ...] = (),
    _authority: EstablishmentAuthority | None = None,
) -> IntensityScaleState:
    """Return state marked as established through approved internal authority."""

    if not isinstance(cast(object, state), IntensityScaleState):
        raise InvalidTransformationStateError(
            "intensity scale state establishment requires an IntensityScaleState value"
        )
    return state.with_establishment(
        established_via=established_via,
        authority=_authority,
        establishment_mode=establishment_mode,
        evidence_level=evidence_level,
        transformer_name=transformer_name,
        input_declaration_source=input_declaration_source,
        parameters={} if parameters is None else parameters,
        trace_id=trace_id,
        diagnostic_warnings=diagnostic_warnings,
    )


def _copy_intensity_scale_establishment(
    *,
    source: IntensityScaleState,
    target: IntensityScaleState,
) -> None:
    if not source.is_established:
        return
    established_via = source.established_via
    authority_source = source.establishment_authority_source
    establishment_mode = source.establishment_mode
    establishment_provenance = source.establishment_provenance
    if (
        established_via is None
        or authority_source is None
        or establishment_mode is None
        or establishment_provenance is None
    ):
        return
    object.__setattr__(target, "_established_via", established_via)
    object.__setattr__(target, "_establishment_authority_source", authority_source)
    object.__setattr__(
        target,
        "_establishment_marker",
        _ESTABLISHED_INTENSITY_SCALE_STATE_MARKER,
    )
    object.__setattr__(target, "_establishment_mode", establishment_mode)
    object.__setattr__(target, "_establishment_provenance", establishment_provenance)


def _validate_quantitative_meaning_transition_contract(
    *,
    current_quantity: QuantitativeMeaning | None,
    target_quantity: QuantitativeMeaning,
    provenance: QuantitativeMeaningTransitionProvenance,
    kind: IntensityScaleKind,
) -> None:
    source_quantity = cast(QuantitativeMeaning | None, provenance.source_quantity)
    provenance_target = cast(QuantitativeMeaning, provenance.target_quantity)
    _validate_quantitative_meaning_source_target_scale_contract(
        source_quantity=source_quantity,
        target_quantity=target_quantity,
        kind=kind,
    )
    if provenance_target is not target_quantity:
        raise InvalidTransformationStateError(
            "quantitative meaning transition provenance target_quantity must match "
            "the requested target_quantity"
        )
    if (
        current_quantity is not None
        and source_quantity is None
        and provenance.evidence_mode
        is not QuantitativeMeaningEvidenceMode.LEGACY_UNVERIFIED
    ):
        raise InvalidTransformationStateError(
            "quantitative meaning transition provenance source_quantity may be "
            "None only for initial establishment"
        )


def _validate_quantitative_meaning_source_target_scale_contract(
    *,
    source_quantity: QuantitativeMeaning | None,
    target_quantity: QuantitativeMeaning,
    kind: IntensityScaleKind,
) -> None:
    if source_quantity is not None:
        _validate_quantitative_meaning_kind_coherence(
            quantity=source_quantity,
            kind=kind,
        )
    _validate_quantitative_meaning_kind_coherence(
        quantity=target_quantity,
        kind=kind,
    )


def _normalize_required_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise InvalidTransformationStateError(
            f"intensity transformation event {field_name} must be a non-empty string"
        )
    normalized = value.strip()
    if not normalized:
        raise InvalidTransformationStateError(
            f"intensity transformation event {field_name} must be a non-empty string"
        )
    return normalized


def _normalize_required_provenance_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise InvalidTransformationStateError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise InvalidTransformationStateError(
            f"{field_name} must be a non-empty string"
        )
    return normalized


def _normalize_optional_text(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        return None
    return normalized


def _normalize_intensity_scale_evidence_level(
    evidence_level: IntensityScaleEvidenceLevel | str,
) -> IntensityScaleEvidenceLevel:
    raw_evidence_level = cast(object, evidence_level)
    if isinstance(raw_evidence_level, IntensityScaleEvidenceLevel):
        return raw_evidence_level
    try:
        return IntensityScaleEvidenceLevel(str(raw_evidence_level).strip())
    except ValueError as exc:
        supported = ", ".join(member.value for member in IntensityScaleEvidenceLevel)
        raise InvalidTransformationStateError(
            "unsupported intensity-scale evidence level "
            f"{raw_evidence_level!r}; supported: {supported}"
        ) from exc


def _normalize_quantitative_meaning_evidence_mode(
    evidence_mode: QuantitativeMeaningEvidenceMode | str,
) -> QuantitativeMeaningEvidenceMode:
    raw_evidence_mode = cast(object, evidence_mode)
    if isinstance(raw_evidence_mode, QuantitativeMeaningEvidenceMode):
        return raw_evidence_mode
    try:
        return QuantitativeMeaningEvidenceMode(str(raw_evidence_mode).strip())
    except ValueError as exc:
        supported = ", ".join(
            member.value for member in QuantitativeMeaningEvidenceMode
        )
        raise InvalidTransformationStateError(
            "unsupported quantitative meaning evidence mode "
            f"{raw_evidence_mode!r}; supported: {supported}"
        ) from exc


def _normalize_quantitative_meaning_operation_id(value: object) -> str:
    raw = _normalize_required_provenance_text(
        value,
        field_name="quantitative_meaning_provenance.operation_id",
    )
    normalized = raw.strip().lower().replace("-", "_")
    normalized = re.sub(r"[\s_]+", "_", normalized)
    normalized = re.sub(r"\.+", ".", normalized)
    if not _QUANTITATIVE_MEANING_OPERATION_PATTERN.fullmatch(normalized):
        raise InvalidTransformationStateError(
            "quantitative_meaning_provenance.operation_id must be a stable "
            "dot-separated machine identifier using lowercase ASCII letters, "
            "digits, and underscores"
        )
    return normalized


def _normalize_diagnostic_caveat_codes(values: object) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise InvalidTransformationStateError(
            "quantitative_meaning_provenance.diagnostic_caveat_codes must be a "
            "sequence of strings"
        )
    if not isinstance(values, Sequence):
        raise InvalidTransformationStateError(
            "quantitative_meaning_provenance.diagnostic_caveat_codes must be a "
            "sequence of strings"
        )
    raw_values = tuple(cast(Sequence[object], values))
    codes: list[str] = []
    for code in raw_values:
        normalized = str(code).strip().lower().replace("-", "_").replace(" ", "_")
        normalized = re.sub(r"_+", "_", normalized)
        if not normalized:
            continue
        if not re.fullmatch(r"[a-z0-9][a-z0-9_]*", normalized):
            raise InvalidTransformationStateError(
                "quantitative_meaning_provenance.diagnostic_caveat_codes must "
                "contain stable lowercase ASCII caveat codes"
            )
        if normalized not in codes:
            codes.append(normalized)
    return tuple(codes)


def _canonical_table_fingerprint_payload_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[Mapping[str, object], ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise InvalidTransformationStateError(f"{field_name} must be a sequence")
    if not isinstance(values, Sequence):
        raise InvalidTransformationStateError(f"{field_name} must be a sequence")
    raw_values = tuple(cast(Sequence[object], values))
    return tuple(
        _canonical_table_fingerprint_payload(
            value,
            field_name=f"{field_name}[{position}]",
        )
        for position, value in enumerate(raw_values)
    )


def _canonical_table_fingerprint_payload(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    try:
        if isinstance(value, TableFingerprint):
            payload = table_fingerprint_to_payload(value)
        elif isinstance(value, Mapping):
            fingerprint = table_fingerprint_from_payload(
                cast(Mapping[str, object], value)
            )
            payload = table_fingerprint_to_payload(fingerprint)
        else:
            raise InvalidTransformationStateError(
                f"{field_name} must be a TableFingerprint or table fingerprint payload"
            )
    except PhosPyInputError as exc:
        raise InvalidTransformationStateError(
            f"{field_name} must be a valid table fingerprint payload: {exc}"
        ) from exc
    return cast(
        Mapping[str, object],
        freeze_json_mapping_with_error_type(
            payload,
            field_name=field_name,
            error_type=InvalidTransformationStateError,
        ),
    )


def _require_payload_mapping(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        raw_mapping = cast(Mapping[object, object], value)
        result: dict[str, object] = {}
        for key, item in raw_mapping.items():
            if not isinstance(key, str):
                raise InvalidTransformationStateError(
                    f"{field_name} JSON object keys must be strings; "
                    f"got {type(key).__name__}"
                )
            if key in result:
                raise InvalidTransformationStateError(
                    f"{field_name} contains duplicate JSON object key {key!r}"
                )
            result[key] = item
        return result
    raise InvalidTransformationStateError(f"{field_name} must be an object")


def _require_payload_sequence(value: object, *, field_name: str) -> list[object]:
    if isinstance(value, (str, bytes, bytearray)):
        raise InvalidTransformationStateError(f"{field_name} must be an array")
    if isinstance(value, Sequence):
        return list(cast(Sequence[object], value))
    raise InvalidTransformationStateError(f"{field_name} must be an array")


def _require_payload_str(value: object, *, field_name: str) -> str:
    return _normalize_required_provenance_text(value, field_name=field_name)


def _optional_payload_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_payload_str(value, field_name=field_name)


def _require_payload_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidTransformationStateError(f"{field_name} must be an int")
    return int(value)


def _normalize_intensity_transformation_pseudocount(
    pseudocount: float | None,
) -> float | None:
    if pseudocount is None:
        return None
    try:
        resolved = float(pseudocount)
    except (TypeError, ValueError) as exc:
        raise InvalidTransformationStateError(
            "intensity transformation event pseudocount must be a finite number"
        ) from exc
    if not math.isfinite(resolved):
        raise InvalidTransformationStateError(
            "intensity transformation event pseudocount must be finite"
        )
    if resolved < 0:
        raise InvalidTransformationStateError(
            "intensity transformation event pseudocount must be greater than or "
            "equal to 0"
        )
    return resolved


def _validate_intensity_transformation_event_transition(
    *,
    input_scale: MatrixIntensityScaleState,
    output_scale: MatrixIntensityScaleState,
    evidence_level: IntensityScaleEvidenceLevel,
    transformation_kind: str,
) -> None:
    normalized_kind = transformation_kind.lower().replace("-", "_")
    if normalized_kind in {
        "identity",
        "passthrough",
        "pass_through",
        "declaration",
        "declared",
        "declared_by_user",
    }:
        if not _same_matrix_intensity_scale(input_scale, output_scale):
            raise InvalidTransformationStateError(
                "inconsistent intensity transformation event scale transition: "
                f"{transformation_kind!r} requires matching input and output scales"
            )
        return

    if normalized_kind in {
        "log2",
        "log2_transform",
        "log2_transformation",
        "linear_to_log2",
    }:
        if (
            input_scale.kind is not IntensityScaleKind.LINEAR
            or input_scale.transformed
            or output_scale.kind is not IntensityScaleKind.LOG2
            or not output_scale.transformed
        ):
            raise InvalidTransformationStateError(
                "inconsistent intensity transformation event scale transition: "
                "log2 requires linear input scale and log2 output scale"
            )
        return

    if _same_matrix_intensity_scale(input_scale, output_scale):
        return

    evidence = evidence_level.value
    raise InvalidTransformationStateError(
        "inconsistent intensity transformation event scale transition: "
        f"{transformation_kind!r} with evidence level {evidence!r} does not define "
        f"a supported transition from {input_scale.kind.value!r} to "
        f"{output_scale.kind.value!r}"
    )


def _same_matrix_intensity_scale(
    left: MatrixIntensityScaleState,
    right: MatrixIntensityScaleState,
) -> bool:
    return left.kind is right.kind and left.transformed == right.transformed


def _matrix_intensity_scale_payload(
    state: MatrixIntensityScaleState,
) -> dict[str, object]:
    return {
        "kind": state.kind.value,
        "transformed": bool(state.transformed),
        "established_by": state.established_by,
    }


def _normalize_quantitative_meaning(
    quantity: QuantitativeMeaning | str | None,
) -> QuantitativeMeaning | None:
    if quantity is None:
        return None
    if isinstance(quantity, QuantitativeMeaning):
        return quantity
    try:
        return QuantitativeMeaning(str(quantity))
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise InvalidTransformationStateError(
            "unsupported intensity-scale quantitative meaning "
            f"'{quantity}'; supported: {supported}"
        ) from exc


def _normalize_optional_quantitative_meaning(
    quantity: QuantitativeMeaning | str | None,
) -> QuantitativeMeaning | None:
    return _normalize_quantitative_meaning(quantity)


def _normalize_required_quantitative_meaning(
    quantity: QuantitativeMeaning | str,
    *,
    field_name: str,
) -> QuantitativeMeaning:
    normalized = _normalize_quantitative_meaning(quantity)
    if normalized is None:
        raise InvalidTransformationStateError(f"{field_name} must not be None")
    return normalized


def default_quantitative_meaning_for_scale_kind(
    kind: IntensityScaleKind,
) -> QuantitativeMeaning:
    """Return the base quantitative meaning implied by a scale contract."""

    if kind is IntensityScaleKind.LINEAR:
        return QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE
    if kind is IntensityScaleKind.LOG2:
        return QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE
    return QuantitativeMeaning.UNKNOWN


def is_caller_declarable_quantitative_meaning(
    meaning: QuantitativeMeaning | str,
) -> bool:
    """Return whether a public caller may declare this direct input meaning."""

    normalized = _normalize_required_quantitative_meaning(
        meaning,
        field_name="quantitative_meaning",
    )
    return normalized in CALLER_DECLARABLE_QUANTITATIVE_MEANINGS


def caller_declarable_quantitative_meaning_values() -> tuple[str, ...]:
    """Return stable public caller-declarable quantitative meaning values."""

    return tuple(
        meaning.value
        for meaning in (
            QuantitativeMeaning.UNKNOWN,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        )
    )


def _validate_quantitative_meaning_kind_coherence(
    *,
    quantity: QuantitativeMeaning,
    kind: IntensityScaleKind,
) -> None:
    rule = _QUANTITATIVE_MEANING_SCALE_RULE_BY_MEANING[quantity]
    if kind in rule.allowed_scales:
        return
    allowed = ", ".join(sorted(scale.value for scale in rule.allowed_scales))
    if len(rule.allowed_scales) == 1:
        allowed = f"{allowed} intensity scale"
    else:
        allowed = f"one of these intensity scales: {allowed}"
    if quantity is QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE:
        raise InvalidTransformationStateError(
            "quantitative meaning 'phosphosite_abundance' requires linear "
            "intensity scale"
        )
    raise InvalidTransformationStateError(
        f"quantitative meaning '{quantity.value}' requires {allowed}"
    )


def _resolve_establishment_source(
    *,
    authority_source: str,
    establishment_mode: IntensityScaleEstablishmentMode,
) -> IntensityScaleEstablishmentSource:
    if authority_source == "phospy.io.bundles._shared.intensity_scale_state":
        return IntensityScaleEstablishmentSource.RESTORED_FROM_TRUSTED_PROVENANCE
    if establishment_mode is IntensityScaleEstablishmentMode.DECLARED:
        return IntensityScaleEstablishmentSource.DECLARED_BY_USER
    return IntensityScaleEstablishmentSource.TRANSFORMED_BY_PHOSPY
