"""Private preprocessing transformation-event sequence validation."""

from __future__ import annotations

from phospy.errors.build import DatasetBuildError
from phospy.science.transformations.models import (
    IntensityTransformationEvent,
    MatrixIntensityScaleState,
)


class _TransformationEventSequenceValidator:
    """Validate continuity between stage-emitted intensity transformation events."""

    def __init__(self) -> None:
        self._previous_event: IntensityTransformationEvent | None = None
        self._previous_stage: str | None = None

    def run(
        self,
        *,
        stage_key: str,
        event: object | None,
    ) -> IntensityTransformationEvent | None:
        current_event = _normalize_intensity_transformation_event(
            stage_key=stage_key,
            event=event,
        )
        if current_event is None:
            return None
        self._validate_current_event(stage_key=stage_key, event=current_event)
        self._previous_event = current_event
        self._previous_stage = stage_key
        return current_event

    def _validate_current_event(
        self,
        *,
        stage_key: str,
        event: IntensityTransformationEvent,
    ) -> None:
        previous_event = self._previous_event
        if previous_event is None:
            return
        if _same_matrix_intensity_scale(previous_event.output_scale, event.input_scale):
            return
        previous_label = (
            "unknown" if self._previous_stage is None else self._previous_stage
        )
        raise DatasetBuildError(
            "dataset preprocessing intensity transformation event conflict: "
            f"stage {stage_key!r} declares input scale "
            f"{_format_matrix_intensity_scale(event.input_scale)} but previous "
            f"event from stage {previous_label!r} produced output scale "
            f"{_format_matrix_intensity_scale(previous_event.output_scale)}"
        )


def _normalize_intensity_transformation_event(
    *,
    stage_key: str,
    event: object | None,
) -> IntensityTransformationEvent | None:
    if event is None:
        return None
    if not isinstance(event, IntensityTransformationEvent):
        raise DatasetBuildError(
            "dataset preprocessing intensity transformation event parse error: "
            f"stage={stage_key!r}, expected IntensityTransformationEvent or None, "
            f"got {event!r} ({type(event).__name__})"
        )
    return event


def _same_matrix_intensity_scale(
    left: MatrixIntensityScaleState,
    right: MatrixIntensityScaleState,
) -> bool:
    return left.kind is right.kind and left.transformed == right.transformed


def _format_matrix_intensity_scale(state: MatrixIntensityScaleState) -> str:
    transformed_label = "transformed" if state.transformed else "untransformed"
    return f"{state.kind.value} ({transformed_label})"
