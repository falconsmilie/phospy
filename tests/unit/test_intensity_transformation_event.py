from __future__ import annotations

import pytest

from phospy.errors.transformations import InvalidTransformationStateError
from phospy.science.transformations import (
    IntensityScaleEvidenceLevel,
    IntensityTransformationEvent,
    MatrixIntensityScaleState,
)


def _linear_scale() -> MatrixIntensityScaleState:
    return MatrixIntensityScaleState.linear(established_by="tests.input")


def _log2_scale() -> MatrixIntensityScaleState:
    return MatrixIntensityScaleState.log2(established_by="tests.log2")


def test_valid_log2_transformation_event() -> None:
    event = IntensityTransformationEvent(
        transformer_name="tests.Log2Transformer",
        input_scale=_linear_scale(),
        output_scale=_log2_scale(),
        evidence_level=IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION,
        transformation_kind="log2",
        pseudocount=1,
        input_fingerprint=" input-fingerprint ",
        output_fingerprint="output-fingerprint",
    )

    assert event.evidence_level is IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION
    assert event.pseudocount == 1.0
    assert event.input_fingerprint == "input-fingerprint"
    assert event.output_scale.kind.value == "log2"
    assert event.to_payload()["evidence_level"] == "observed_transformation"


def test_valid_declared_by_user_event() -> None:
    scale = _log2_scale()

    event = IntensityTransformationEvent(
        transformer_name="dataset_build_request.input_intensity_scale",
        input_scale=scale,
        output_scale=scale,
        evidence_level=IntensityScaleEvidenceLevel.DECLARED_BY_USER,
        transformation_kind="declared_by_user",
    )

    assert event.evidence_level is IntensityScaleEvidenceLevel.DECLARED_BY_USER
    assert event.input_scale is scale
    assert event.output_scale is scale


def test_intensity_transformation_event_rejects_missing_transformer_name() -> None:
    with pytest.raises(InvalidTransformationStateError, match="transformer_name"):
        IntensityTransformationEvent(
            transformer_name=" ",
            input_scale=_linear_scale(),
            output_scale=_log2_scale(),
            evidence_level=IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION,
            transformation_kind="log2",
        )


def test_observed_intensity_transformation_event_rejects_unknown_output_scale() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="observed intensity transformation event requires a known output scale",
    ):
        IntensityTransformationEvent(
            transformer_name="tests.Log2Transformer",
            input_scale=_linear_scale(),
            output_scale=None,  # type: ignore[arg-type]
            evidence_level=IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION,
            transformation_kind="log2",
        )


def test_intensity_transformation_event_rejects_inconsistent_transition() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="log2 requires linear input scale and log2 output scale",
    ):
        IntensityTransformationEvent(
            transformer_name="tests.Log2Transformer",
            input_scale=_log2_scale(),
            output_scale=_linear_scale(),
            evidence_level=IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION,
            transformation_kind="log2",
        )


@pytest.mark.parametrize("pseudocount", [-1.0, float("nan"), float("inf")])
def test_intensity_transformation_event_rejects_invalid_pseudocount(
    pseudocount: float,
) -> None:
    with pytest.raises(InvalidTransformationStateError, match="pseudocount"):
        IntensityTransformationEvent(
            transformer_name="tests.Log2Transformer",
            input_scale=_linear_scale(),
            output_scale=_log2_scale(),
            evidence_level=IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION,
            transformation_kind="log2",
            pseudocount=pseudocount,
        )
