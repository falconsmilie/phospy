"""Shared workflow provenance for input intensity-scale evidence."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

from phospy.contracts.result_caveats import ResultCaveat
from phospy.provenance.models import InputIntensityScaleEvidence
from phospy.science.transformations.models import (
    IntensityScaleEvidenceLevel,
    IntensityScaleState,
)

if TYPE_CHECKING:
    from phospy.science.datasets.models import AnalysisReadyPhosphoDataset

INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE = "input_intensity_scale_declared_by_user"
INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_MESSAGE = (
    "Input intensity scale was accepted from user-declared metadata rather than "
    "an observed transformation event. Downstream quantitative results depend "
    "on that declaration being correct."
)


def input_intensity_scale_evidence_from_dataset(
    dataset: AnalysisReadyPhosphoDataset,
) -> InputIntensityScaleEvidence:
    """Build workflow-visible intensity-scale evidence from a dataset."""

    return input_intensity_scale_evidence_from_state(dataset.intensity_scale_state)


def input_intensity_scale_evidence_from_state(
    state: IntensityScaleState,
) -> InputIntensityScaleEvidence:
    """Build workflow-visible intensity-scale evidence from scale state."""

    provenance = state.establishment_provenance
    if provenance is None:
        return InputIntensityScaleEvidence(
            input_intensity_scale=state.label,
            input_intensity_scale_evidence_level=(
                IntensityScaleEvidenceLevel.UNKNOWN.value
            ),
            input_intensity_scale_source="unknown",
        )
    return InputIntensityScaleEvidence(
        input_intensity_scale=provenance.scale,
        input_intensity_scale_evidence_level=provenance.evidence_level.value,
        input_intensity_scale_source=provenance.source.value,
        input_intensity_scale_source_detail=_source_detail(state),
    )


def input_intensity_scale_evidence_payload(
    dataset: AnalysisReadyPhosphoDataset,
) -> dict[str, object]:
    """Return the JSON-safe workflow evidence payload for a dataset."""

    return input_intensity_scale_evidence_from_dataset(dataset).to_payload()


def with_input_intensity_scale_evidence(
    parameters: Mapping[str, object] | None,
    *,
    dataset: AnalysisReadyPhosphoDataset,
) -> dict[str, object]:
    """Merge required intensity-scale evidence into workflow parameters."""

    payload = (
        {}
        if parameters is None
        else {str(key): value for key, value in parameters.items()}
    )
    payload.update(input_intensity_scale_evidence_payload(dataset))
    return payload


def build_declared_input_intensity_scale_caveat(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    workflow_scope: str,
) -> ResultCaveat | None:
    """Return the shared caveat for declared-only input scale evidence."""

    evidence = input_intensity_scale_evidence_from_dataset(dataset)
    if (
        evidence.input_intensity_scale_evidence_level
        != IntensityScaleEvidenceLevel.DECLARED_BY_USER.value
    ):
        return None
    details = evidence.to_payload()
    details["workflow_scope"] = str(workflow_scope)
    return ResultCaveat(
        code=INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE,
        severity="warning",
        message=INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_MESSAGE,
        details=details,
    )


def _source_detail(state: IntensityScaleState) -> str | None:
    provenance = state.establishment_provenance
    if provenance is None:
        return None
    if provenance.input_declaration_source is not None:
        return provenance.input_declaration_source
    if provenance.transformer_name is not None:
        return provenance.transformer_name
    return state.establishment_authority_source


__all__ = [
    "INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE",
    "INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_MESSAGE",
    "build_declared_input_intensity_scale_caveat",
    "input_intensity_scale_evidence_from_dataset",
    "input_intensity_scale_evidence_from_state",
    "input_intensity_scale_evidence_payload",
    "with_input_intensity_scale_evidence",
]
