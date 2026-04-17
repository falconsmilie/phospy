"""Transformation domain models."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from phospy.errors.transformations import InvalidTransformationStateError


class TransformationKind(str, Enum):
    """Supported transformation-state kinds."""

    LINEAR = "linear"
    LOG2 = "log2"


@dataclass(frozen=True, slots=True)
class MatrixTransformationState:
    """Established transformation state for one quantitative matrix."""

    kind: TransformationKind
    transformed: bool
    established_by: str = "phospy.transformations.transformers.identity"

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
        established_by: str = "phospy.transformations.transformers.identity",
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

    @classmethod
    def raw(cls, *, has_total_matrix: bool = False) -> TransformationState:
        """Create the canonical untransformed state used by the identity path."""

        if has_total_matrix:
            return cls(
                phospho=MatrixTransformationState.linear(),
                total=MatrixTransformationState.linear(),
            )
        return cls(phospho=MatrixTransformationState.linear(), total=None)
