"""Transformation-domain intensity scale models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import Final

from phospy.errors.transformations import InvalidTransformationStateError
from phospy.transformations._authority import (
    _EstablishmentAuthority,
    _resolve_establishment_authority_source,
)

IDENTITY_INTENSITY_SCALE_ESTABLISHER: Final[str] = (
    "phospy.transformations.transformers.identity"
)
_ESTABLISHED_INTENSITY_SCALE_STATE_MARKER: Final[object] = object()


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
    transformer_name: str | None = None
    input_declaration_source: str | None = None
    parameters: Mapping[str, object] = field(default_factory=dict)
    trace_id: str | None = None
    diagnostic_warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        scale = str(self.scale).strip()
        if not scale:
            raise InvalidTransformationStateError(
                "intensity-scale establishment provenance requires non-empty scale"
            )
        object.__setattr__(self, "scale", scale)
        if not isinstance(self.mode, IntensityScaleEstablishmentMode):
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
            {str(key): value for key, value in self.parameters.items()},
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
        established_by: str = "phospy.transformations.transformers.log2",
    ) -> MatrixIntensityScaleState:
        """Construct an established log2-transformed matrix state."""

        return cls(
            kind=IntensityScaleKind.LOG2,
            transformed=True,
            established_by=established_by,
        )


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
        if not isinstance(self.phospho, MatrixIntensityScaleState):
            raise InvalidTransformationStateError(
                "intensity_scale_state.phospho must be a MatrixIntensityScaleState"
            )
        if self.total is not None and not isinstance(
            self.total, MatrixIntensityScaleState
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

    @classmethod
    def raw(cls, *, has_total_matrix: bool = False) -> IntensityScaleState:
        """Create a canonical declared linear state (not yet established)."""

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
        _authority: _EstablishmentAuthority | None = None,
    ) -> IntensityScaleState:
        """Create a canonical linear state through an approved internal authority."""

        return cls.raw(has_total_matrix=has_total_matrix)._with_establishment(
            established_via=established_via,
            authority=_authority,
            establishment_mode=establishment_mode,
            transformer_name=transformer_name,
            input_declaration_source=input_declaration_source,
            parameters={} if parameters is None else parameters,
            trace_id=trace_id,
            diagnostic_warnings=diagnostic_warnings,
        )

    def _with_establishment(
        self,
        *,
        established_via: str,
        authority: _EstablishmentAuthority | None,
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
        if isinstance(establishment_mode, IntensityScaleEstablishmentMode):
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
        authority_source = _resolve_establishment_authority_source(authority)
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
                transformer_name=transformer_name,
                input_declaration_source=input_declaration_source,
                parameters=parameters,
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
    _authority: _EstablishmentAuthority | None = None,
) -> IntensityScaleState:
    """Return state marked as established through approved internal authority."""

    if not isinstance(state, IntensityScaleState):
        raise InvalidTransformationStateError(
            "intensity scale state establishment requires an IntensityScaleState value"
        )
    return state._with_establishment(
        established_via=established_via,
        authority=_authority,
        establishment_mode=establishment_mode,
        transformer_name=transformer_name,
        input_declaration_source=input_declaration_source,
        parameters={} if parameters is None else parameters,
        trace_id=trace_id,
        diagnostic_warnings=diagnostic_warnings,
    )


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
