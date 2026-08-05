"""Transformation-domain immutable scale value models."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import cast

from phospy.errors.transformations import InvalidTransformationStateError
from phospy.science.transformations.policy import (
    IDENTITY_INTENSITY_SCALE_ESTABLISHER,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    normalize_intensity_scale_evidence_level,
)


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

        evidence_level = normalize_intensity_scale_evidence_level(self.evidence_level)
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


__all__ = [
    "IntensityTransformationEvent",
    "MatrixIntensityScaleState",
]
