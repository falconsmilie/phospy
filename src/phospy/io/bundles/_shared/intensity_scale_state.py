"""Intensity-scale-state payload serialization for bundle manifests."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import (
    require_bool,
    require_mapping,
    require_str,
)
from phospy.transformations._authority import (
    _bundle_reconstruction_establishment_authority,
)
from phospy.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    establish_intensity_scale_state,
)


def intensity_scale_state_to_payload(state: IntensityScaleState) -> dict[str, object]:
    """Serialize intensity scale state to manifest payload."""

    return {
        "phospho": _matrix_state_to_payload(state.phospho),
        "total": None if state.total is None else _matrix_state_to_payload(state.total),
    }


def intensity_scale_state_from_payload(
    payload: Mapping[str, object],
) -> IntensityScaleState:
    """Deserialize intensity scale state from manifest payload."""

    phospho_payload = require_mapping(
        payload.get("phospho"),
        field_name="dataset.metadata.intensity_scale_state.phospho",
    )
    total_raw = payload.get("total")
    if total_raw is None:
        total_state = None
    else:
        total_state = _matrix_state_from_payload(
            require_mapping(
                total_raw,
                field_name="dataset.metadata.intensity_scale_state.total",
            )
        )
    state = IntensityScaleState(
        phospho=_matrix_state_from_payload(phospho_payload),
        total=total_state,
    )
    return establish_intensity_scale_state(
        state,
        established_via="phospy.io.bundles._shared.intensity_scale_state",
        _authority=_bundle_reconstruction_establishment_authority(),
    )


def _matrix_state_to_payload(state: MatrixIntensityScaleState) -> dict[str, object]:
    return {
        "kind": state.kind.value,
        "transformed": state.transformed,
        "established_by": state.established_by,
    }


def _matrix_state_from_payload(
    payload: Mapping[str, object],
) -> MatrixIntensityScaleState:
    kind_token = require_str(
        payload.get("kind"), field_name="matrix_intensity_scale_state.kind"
    )
    try:
        kind = IntensityScaleKind(kind_token)
    except ValueError as exc:
        supported = ", ".join(member.value for member in IntensityScaleKind)
        raise PhosPyInputError(
            f"unsupported intensity scale kind '{kind_token}'; supported: {supported}"
        ) from exc
    return MatrixIntensityScaleState(
        kind=kind,
        transformed=require_bool(
            payload.get("transformed"),
            field_name="matrix_intensity_scale_state.transformed",
        ),
        established_by=require_str(
            payload.get("established_by"),
            field_name="matrix_intensity_scale_state.established_by",
        ),
    )
