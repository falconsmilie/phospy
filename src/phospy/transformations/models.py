"""Transformation domain models."""

from __future__ import annotations

import inspect
from dataclasses import dataclass, field
from enum import Enum
from typing import Final

from phospy.errors.transformations import InvalidTransformationStateError

IDENTITY_TRANSFORMATION_ESTABLISHER: Final[str] = (
    "phospy.transformations.transformers.identity"
)
_SUPPORTED_ESTABLISHMENT_CALLER_PREFIXES: Final[tuple[str, ...]] = (
    "phospy.datasets.builders.transformation_resolver",
    "phospy.io.bundles._shared.transformation_state",
    "phospy.transformations.transformers",
)
_SUPPORTED_ESTABLISHMENT_MARKER: Final[object] = object()


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
    _supported_establishment_marker: object | None = field(
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
            self._supported_establishment_marker is _SUPPORTED_ESTABLISHMENT_MARKER
            and self._established_via is not None
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
    ) -> TransformationState:
        """Create a canonical linear state through a supported establisher lane."""

        _require_supported_establishment_caller()
        return cls.raw(has_total_matrix=has_total_matrix)._with_establishment(
            established_via=established_via,
            supported_establishment_marker=_SUPPORTED_ESTABLISHMENT_MARKER,
        )

    def _with_establishment(
        self,
        *,
        established_via: str,
        supported_establishment_marker: object,
    ) -> TransformationState:
        source = established_via.strip()
        if not source:
            raise InvalidTransformationStateError(
                "transformation state establishment source must be a non-empty string"
            )
        if supported_establishment_marker is not _SUPPORTED_ESTABLISHMENT_MARKER:
            raise InvalidTransformationStateError(
                "transformation state establishment requires a supported-path marker"
            )
        established = TransformationState(
            phospho=self.phospho,
            total=self.total,
        )
        object.__setattr__(established, "_established_via", source)
        object.__setattr__(
            established,
            "_supported_establishment_marker",
            supported_establishment_marker,
        )
        return established


def establish_transformation_state(
    state: TransformationState,
    *,
    established_via: str,
) -> TransformationState:
    """Return a state marked as established through a supported source lane."""

    if not isinstance(state, TransformationState):
        raise InvalidTransformationStateError(
            "transformation state establishment requires a TransformationState value"
        )
    _require_supported_establishment_caller()
    return state._with_establishment(
        established_via=established_via,
        supported_establishment_marker=_SUPPORTED_ESTABLISHMENT_MARKER,
    )


def _require_supported_establishment_caller() -> None:
    caller_module = _caller_module_name()
    if _is_supported_establishment_caller(caller_module):
        return
    raise InvalidTransformationStateError(
        "transformation state can be established only through supported PhosPy "
        "builder/transformer or bundle reconstruction paths"
    )


def _caller_module_name() -> str | None:
    frame = inspect.currentframe()
    if frame is None:
        return None
    try:
        target_frame = frame
        for _ in range(3):
            target_frame = target_frame.f_back
            if target_frame is None:
                return None
        module_name = target_frame.f_globals.get("__name__")
        return module_name if isinstance(module_name, str) else None
    finally:
        del frame


def _is_supported_establishment_caller(module_name: str | None) -> bool:
    if module_name is None:
        return False
    return any(
        module_name == prefix or module_name.startswith(f"{prefix}.")
        for prefix in _SUPPORTED_ESTABLISHMENT_CALLER_PREFIXES
    )
