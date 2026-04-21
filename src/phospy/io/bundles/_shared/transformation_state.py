"""Transformation-state payload serialization for bundle manifests."""

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
    MatrixTransformationState,
    TransformationKind,
    TransformationState,
    establish_transformation_state,
)


def transformation_state_to_payload(state: TransformationState) -> dict[str, object]:
    """Serialize transformation state to manifest payload."""

    return {
        "phospho": _matrix_state_to_payload(state.phospho),
        "total": None if state.total is None else _matrix_state_to_payload(state.total),
    }


def transformation_state_from_payload(
    payload: Mapping[str, object],
) -> TransformationState:
    """Deserialize transformation state from manifest payload."""

    phospho_payload = require_mapping(
        payload.get("phospho"),
        field_name="dataset.metadata.transformation_state.phospho",
    )
    total_raw = payload.get("total")
    if total_raw is None:
        total_state = None
    else:
        total_state = _matrix_state_from_payload(
            require_mapping(
                total_raw,
                field_name="dataset.metadata.transformation_state.total",
            )
        )
    state = TransformationState(
        phospho=_matrix_state_from_payload(phospho_payload),
        total=total_state,
    )
    return establish_transformation_state(
        state,
        established_via="phospy.io.bundles._shared.transformation_state",
        _authority=_bundle_reconstruction_establishment_authority(),
    )


def _matrix_state_to_payload(state: MatrixTransformationState) -> dict[str, object]:
    return {
        "kind": state.kind.value,
        "transformed": state.transformed,
        "established_by": state.established_by,
    }


def _matrix_state_from_payload(
    payload: Mapping[str, object],
) -> MatrixTransformationState:
    kind_token = require_str(
        payload.get("kind"), field_name="matrix_transformation_state.kind"
    )
    try:
        kind = TransformationKind(kind_token)
    except ValueError as exc:
        supported = ", ".join(member.value for member in TransformationKind)
        raise PhosPyInputError(
            f"unsupported transformation kind '{kind_token}'; supported: {supported}"
        ) from exc
    return MatrixTransformationState(
        kind=kind,
        transformed=require_bool(
            payload.get("transformed"),
            field_name="matrix_transformation_state.transformed",
        ),
        established_by=require_str(
            payload.get("established_by"),
            field_name="matrix_transformation_state.established_by",
        ),
    )
