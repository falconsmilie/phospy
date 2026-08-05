"""Transformation-domain scale state establishment and transitions."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Final, cast

from phospy.errors.transformations import InvalidTransformationStateError
from phospy.science.transformations._authority import (
    EstablishmentAuthority,
    QuantitativeMeaningTransitionAuthority,
    resolve_establishment_authority_source,
    resolve_quantitative_meaning_transition_authority_source,
)
from phospy.science.transformations.policy import (
    BUNDLE_QUANTITATIVE_MEANING_AUTHORITY_SOURCE,
    DATASET_QUANTITATIVE_MEANING_AUTHORITY_SOURCE,
    IDENTITY_INTENSITY_SCALE_ESTABLISHER,
    IntensityScaleEstablishmentMode,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    normalize_intensity_scale_evidence_level,
    normalize_quantitative_meaning,
    normalize_required_quantitative_meaning,
    resolve_establishment_source,
    validate_quantitative_meaning_kind_coherence,
)
from phospy.science.transformations.provenance import (
    IntensityScaleEstablishmentProvenance,
    QuantitativeMeaningTransitionProvenance,
)
from phospy.science.transformations.scale_values import MatrixIntensityScaleState

_ESTABLISHED_INTENSITY_SCALE_STATE_MARKER: Final[object] = object()


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
        quantity = normalize_quantitative_meaning(self.quantity)
        if quantity is not None:
            validate_quantitative_meaning_kind_coherence(
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
        resolved_evidence_level = normalize_intensity_scale_evidence_level(
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
                source=resolve_establishment_source(
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
        if authority_source != DATASET_QUANTITATIVE_MEANING_AUTHORITY_SOURCE:
            raise InvalidTransformationStateError(
                "bundle restoration authority cannot mint a new quantitative "
                "meaning transition"
            )
        target = normalize_required_quantitative_meaning(
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
        if authority_source != BUNDLE_QUANTITATIVE_MEANING_AUTHORITY_SOURCE:
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
        validate_quantitative_meaning_kind_coherence(
            quantity=source_quantity,
            kind=kind,
        )
    validate_quantitative_meaning_kind_coherence(
        quantity=target_quantity,
        kind=kind,
    )


__all__ = [
    "IntensityScaleState",
    "establish_intensity_scale_state",
]
