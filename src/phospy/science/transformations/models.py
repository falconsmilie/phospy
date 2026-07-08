"""Transformation-domain intensity scale models."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import Final, cast

from phospy.errors.transformations import InvalidTransformationStateError
from phospy.science.transformations._authority import (
    EstablishmentAuthority,
    resolve_establishment_authority_source,
)

IDENTITY_INTENSITY_SCALE_ESTABLISHER: Final[str] = (
    "phospy.science.transformations.transformers.identity"
)
_ESTABLISHED_INTENSITY_SCALE_STATE_MARKER: Final[object] = object()


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
    MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE = (
        "mixed_phospho_total_log_ratio_and_phosphosite_log_abundance"
    )
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class IntensityScaleEstablishmentProvenance:
    """Structured provenance for intensity-scale establishment."""

    scale: str
    mode: IntensityScaleEstablishmentMode
    source: IntensityScaleEstablishmentSource
    transformer_name: str | None = None
    input_declaration_source: str | None = None
    parameters: dict[str, object] = field(
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
            cast(
                Mapping[str, object],
                {str(key): value for key, value in self.parameters.items()},
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
            "transformer_name": self.transformer_name,
            "input_declaration_source": self.input_declaration_source,
            "parameters": dict(self.parameters),
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
        quantity = _normalize_quantitative_meaning(
            self.quantity,
            default_kind=self.phospho.kind,
        )
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
                transformer_name=transformer_name,
                input_declaration_source=input_declaration_source,
                parameters=dict(parameters),
                trace_id=trace_id,
                diagnostic_warnings=diagnostic_warnings,
            ),
        )
        return established

    def with_quantitative_meaning(
        self,
        quantity: QuantitativeMeaning,
    ) -> IntensityScaleState:
        """Clone this state with a different quantitative meaning."""

        normalized_quantity = _normalize_quantitative_meaning(
            quantity,
            default_kind=self.phospho.kind,
        )
        if normalized_quantity is self.quantity:
            return self
        updated = IntensityScaleState(
            phospho=self.phospho,
            total=self.total,
            quantity=normalized_quantity,
        )
        if not self.is_established:
            return updated
        object.__setattr__(updated, "_established_via", self._established_via)
        object.__setattr__(
            updated,
            "_establishment_authority_source",
            self._establishment_authority_source,
        )
        object.__setattr__(
            updated,
            "_establishment_marker",
            _ESTABLISHED_INTENSITY_SCALE_STATE_MARKER,
        )
        object.__setattr__(updated, "_establishment_mode", self._establishment_mode)
        object.__setattr__(
            updated,
            "_establishment_provenance",
            self._establishment_provenance,
        )
        return updated


def establish_intensity_scale_state(
    state: IntensityScaleState,
    *,
    established_via: str,
    establishment_mode: IntensityScaleEstablishmentMode = (
        IntensityScaleEstablishmentMode.DERIVED
    ),
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
        transformer_name=transformer_name,
        input_declaration_source=input_declaration_source,
        parameters={} if parameters is None else parameters,
        trace_id=trace_id,
        diagnostic_warnings=diagnostic_warnings,
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
    *,
    default_kind: IntensityScaleKind,
) -> QuantitativeMeaning:
    if quantity is None:
        if default_kind is IntensityScaleKind.LINEAR:
            return QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE
        if default_kind is IntensityScaleKind.LOG2:
            return QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE
        return QuantitativeMeaning.UNKNOWN
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


def _validate_quantitative_meaning_kind_coherence(
    *,
    quantity: QuantitativeMeaning,
    kind: IntensityScaleKind,
) -> None:
    if quantity is QuantitativeMeaning.UNKNOWN:
        return
    if (
        quantity is QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE
        and kind is not IntensityScaleKind.LINEAR
    ):
        raise InvalidTransformationStateError(
            "quantitative meaning 'phosphosite_abundance' requires linear "
            "intensity scale"
        )
    if (
        quantity
        in {
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
            QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE,
        }
        and kind is not IntensityScaleKind.LOG2
    ):
        raise InvalidTransformationStateError(
            f"quantitative meaning '{quantity.value}' requires log2 intensity scale"
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
