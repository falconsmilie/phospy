"""Transformation domain models."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Final

from phospy.errors.transformations import InvalidTransformationStateError
from phospy.transformations._authority import (
    _EstablishmentAuthority,
    _resolve_establishment_authority_source,
)

IDENTITY_TRANSFORMATION_ESTABLISHER: Final[str] = (
    "phospy.transformations.transformers.identity"
)
_ESTABLISHED_TRANSFORMATION_STATE_MARKER: Final[object] = object()


class TransformationKind(str, Enum):
    """Supported transformation-state kinds."""

    LINEAR = "linear"
    LOG2 = "log2"


@dataclass(frozen=True, slots=True)
class MatrixTransformationState:
    """Established transformation state for one quantitative matrix."""

    kind: TransformationKind
    transformed: bool
    established_by: str = IDENTITY_TRANSFORMATION_ESTABLISHER

    def __post_init__(self) -> None:
        if self.kind is TransformationKind.LINEAR and self.transformed:
            raise InvalidTransformationStateError(
                "linear matrix state cannot be marked as transformed"
            )
        if self.kind is TransformationKind.LOG2 and not self.transformed:
            raise InvalidTransformationStateError(
                "log2 matrix state must be marked as transformed"
            )
        if not self.established_by.strip():
            raise InvalidTransformationStateError(
                "matrix transformation state requires a non-empty established_by value"
            )

    @classmethod
    def linear(
        cls,
        *,
        established_by: str = IDENTITY_TRANSFORMATION_ESTABLISHER,
    ) -> MatrixTransformationState:
        """Construct an established linear (untransformed) matrix state."""

        return cls(
            kind=TransformationKind.LINEAR,
            transformed=False,
            established_by=established_by,
        )

    @classmethod
    def log2(
        cls,
        *,
        established_by: str = "phospy.transformations.transformers.log2",
    ) -> MatrixTransformationState:
        """Construct an established log2-transformed matrix state."""

        return cls(
            kind=TransformationKind.LOG2,
            transformed=True,
            established_by=established_by,
        )


@dataclass(frozen=True, slots=True)
class TransformationState:
    """Validated transformation metadata for dataset inputs."""

    phospho: MatrixTransformationState = field(
        default_factory=MatrixTransformationState.linear
    )
    total: MatrixTransformationState | None = None
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

    def __post_init__(self) -> None:
        if not isinstance(self.phospho, MatrixTransformationState):
            raise InvalidTransformationStateError(
                "transformation_state.phospho must be a MatrixTransformationState"
            )
        if self.total is not None and not isinstance(
            self.total, MatrixTransformationState
        ):
            raise InvalidTransformationStateError(
                "transformation_state.total must be a MatrixTransformationState or None"
            )

    @property
    def kind(self) -> TransformationKind:
        """Primary transformation kind for phosphosite intensities."""

        return self.phospho.kind

    @property
    def label(self) -> str:
        """Derived human-readable transformation label."""

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
            self._establishment_marker is _ESTABLISHED_TRANSFORMATION_STATE_MARKER
            and self._established_via is not None
            and self._establishment_authority_source is not None
        )

    @classmethod
    def raw(cls, *, has_total_matrix: bool = False) -> TransformationState:
        """Create a canonical declared linear state (not yet established)."""

        if has_total_matrix:
            return cls(
                phospho=MatrixTransformationState.linear(),
                total=MatrixTransformationState.linear(),
            )
        return cls(phospho=MatrixTransformationState.linear(), total=None)

    @classmethod
    def established_raw(
        cls,
        *,
        has_total_matrix: bool = False,
        established_via: str = IDENTITY_TRANSFORMATION_ESTABLISHER,
        _authority: _EstablishmentAuthority | None = None,
    ) -> TransformationState:
        """Create a canonical linear state through an approved internal authority."""

        return cls.raw(has_total_matrix=has_total_matrix)._with_establishment(
            established_via=established_via,
            authority=_authority,
        )

    def _with_establishment(
        self,
        *,
        established_via: str,
        authority: _EstablishmentAuthority | None,
    ) -> TransformationState:
        source = established_via.strip()
        if not source:
            raise InvalidTransformationStateError(
                "transformation state establishment source must be a non-empty string"
            )
        authority_source = _resolve_establishment_authority_source(authority)
        established = TransformationState(
            phospho=self.phospho,
            total=self.total,
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
            _ESTABLISHED_TRANSFORMATION_STATE_MARKER,
        )
        return established


def establish_transformation_state(
    state: TransformationState,
    *,
    established_via: str,
    _authority: _EstablishmentAuthority | None = None,
) -> TransformationState:
    """Return state marked as established through approved internal authority."""

    if not isinstance(state, TransformationState):
        raise InvalidTransformationStateError(
            "transformation state establishment requires a TransformationState value"
        )
    return state._with_establishment(
        established_via=established_via,
        authority=_authority,
    )
