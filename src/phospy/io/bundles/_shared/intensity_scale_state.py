"""Intensity-scale-state payload serialization for bundle manifests."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import InvalidTransformationStateError
from phospy.io.bundles._shared.primitives import (
    require_bool,
    require_mapping,
    require_str,
)
from phospy.science.transformations._authority import (
    _bundle_reconstruction_establishment_authority,
    bundle_quantitative_meaning_restoration_authority,
)
from phospy.science.transformations.models import (
    QUANTITATIVE_MEANING_LEGACY_UNVERIFIED_CAVEAT_CODE,
    QUANTITATIVE_MEANING_OPERATION_LEGACY_BUNDLE_MIGRATION,
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    QuantitativeMeaningTransitionProvenance,
    establish_intensity_scale_state,
)


def intensity_scale_state_to_payload(state: IntensityScaleState) -> dict[str, object]:
    """Serialize intensity scale state to manifest payload."""

    quantity = state.quantity
    if quantity is None:
        raise PhosPyInputError(
            "intensity_scale_state.quantity must be established before serialization"
        )
    quantitative_meaning_provenance = state.quantitative_meaning_provenance
    if quantitative_meaning_provenance is None:
        raise PhosPyInputError(
            "intensity_scale_state.quantitative_meaning_provenance must be "
            "established before serialization"
        )
    if quantitative_meaning_provenance.target_quantity is not quantity:
        raise PhosPyInputError(
            "intensity_scale_state.quantitative_meaning_provenance target must "
            "match intensity_scale_state.quantity"
        )
    return {
        "phospho": _matrix_state_to_payload(state.phospho),
        "total": None if state.total is None else _matrix_state_to_payload(state.total),
        "quantity": quantity.value,
        "quantitative_meaning_provenance": (
            quantitative_meaning_provenance.to_payload()
        ),
    }


def intensity_scale_state_from_payload(
    payload: Mapping[str, object],
    *,
    legacy_quantitative_meaning_policy: str = "reject",
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
    quantity = _quantitative_meaning_from_payload(payload)
    state = IntensityScaleState(
        phospho=_matrix_state_from_payload(phospho_payload),
        total=total_state,
    )
    established = establish_intensity_scale_state(
        state,
        established_via="phospy.io.bundles._shared.intensity_scale_state",
        _authority=_bundle_reconstruction_establishment_authority(),
    )
    provenance = _quantitative_meaning_provenance_from_payload(
        payload,
        quantity=quantity,
        legacy_policy=legacy_quantitative_meaning_policy,
    )
    return established.restore_quantitative_meaning_provenance(
        provenance=provenance,
        authority=bundle_quantitative_meaning_restoration_authority(),
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


def _quantitative_meaning_from_payload(
    payload: Mapping[str, object],
) -> QuantitativeMeaning:
    token = require_str(
        payload.get("quantity"),
        field_name="dataset.metadata.intensity_scale_state.quantity",
    )
    try:
        return QuantitativeMeaning(token)
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise PhosPyInputError(
            "unsupported intensity scale quantitative meaning "
            f"'{token}'; supported: {supported}"
        ) from exc


def _quantitative_meaning_provenance_from_payload(
    payload: Mapping[str, object],
    *,
    quantity: QuantitativeMeaning,
    legacy_policy: str,
) -> QuantitativeMeaningTransitionProvenance:
    raw = payload.get("quantitative_meaning_provenance")
    if raw is None:
        if str(legacy_policy).strip() != "migrate_unverified":
            raise PhosPyInputError(
                "dataset.metadata.intensity_scale_state."
                "quantitative_meaning_provenance is required; pass "
                "legacy_quantitative_meaning_policy='migrate_unverified' only "
                "when loading a historical payload without semantic provenance"
            )
        return QuantitativeMeaningTransitionProvenance(
            source_quantity=None,
            target_quantity=quantity,
            operation_id=QUANTITATIVE_MEANING_OPERATION_LEGACY_BUNDLE_MIGRATION,
            producer_id="phospy.io.bundles._shared.intensity_scale_state",
            evidence_mode=QuantitativeMeaningEvidenceMode.LEGACY_UNVERIFIED,
            parameters={
                "legacy_payload_missing_quantitative_meaning_provenance": True,
            },
            diagnostic_caveat_codes=(
                QUANTITATIVE_MEANING_LEGACY_UNVERIFIED_CAVEAT_CODE,
            ),
        )
    if str(legacy_policy).strip() not in {"reject", "migrate_unverified"}:
        raise PhosPyInputError(
            "legacy_quantitative_meaning_policy must be 'reject' or "
            "'migrate_unverified'"
        )
    provenance_payload = require_mapping(
        raw,
        field_name=(
            "dataset.metadata.intensity_scale_state.quantitative_meaning_provenance"
        ),
    )
    try:
        provenance = QuantitativeMeaningTransitionProvenance.from_payload(
            provenance_payload
        )
    except InvalidTransformationStateError as exc:
        raise PhosPyInputError(
            "dataset.metadata.intensity_scale_state."
            "quantitative_meaning_provenance is invalid: "
            f"{exc}"
        ) from exc
    if provenance.target_quantity is not quantity:
        raise PhosPyInputError(
            "dataset.metadata.intensity_scale_state."
            "quantitative_meaning_provenance.target_quantity must match "
            "dataset.metadata.intensity_scale_state.quantity"
        )
    return provenance
