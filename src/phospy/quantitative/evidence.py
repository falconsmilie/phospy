"""Input intensity-scale evidence and state-only caveat policy."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.provenance.models._shared import (
    _optional_provenance_text,
    _required_provenance_text,
)
from phospy.provenance.result_caveats import ResultCaveat
from phospy.quantitative.models import (
    IntensityScaleEvidenceLevel,
    IntensityScaleState,
)

INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE = "input_intensity_scale_declared_by_user"
INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_MESSAGE = (
    "Input intensity scale was accepted from user-declared metadata rather than "
    "an observed transformation event. Downstream quantitative results depend "
    "on that declaration being correct."
)


@dataclass(frozen=True, slots=True)
class InputIntensityScaleEvidence:
    """Workflow-visible evidence for how input intensity scale was established."""

    input_intensity_scale: str
    input_intensity_scale_evidence_level: str
    input_intensity_scale_source: str
    input_intensity_scale_source_detail: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "input_intensity_scale",
            _required_provenance_text(
                self.input_intensity_scale,
                field_name="input_intensity_scale",
            ),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_evidence_level",
            _required_provenance_text(
                self.input_intensity_scale_evidence_level,
                field_name="input_intensity_scale_evidence_level",
            ),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_source",
            _required_provenance_text(
                self.input_intensity_scale_source,
                field_name="input_intensity_scale_source",
            ),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_source_detail",
            _optional_provenance_text(self.input_intensity_scale_source_detail),
        )

    def to_payload(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "input_intensity_scale": self.input_intensity_scale,
            "input_intensity_scale_evidence_level": (
                self.input_intensity_scale_evidence_level
            ),
            "input_intensity_scale_source": self.input_intensity_scale_source,
        }
        if self.input_intensity_scale_source_detail is not None:
            payload["input_intensity_scale_source_detail"] = (
                self.input_intensity_scale_source_detail
            )
        return payload


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


def build_declared_input_intensity_scale_evidence_caveat(
    *,
    evidence: InputIntensityScaleEvidence,
    workflow_scope: str,
    extra_details: Mapping[str, object] | None = None,
) -> ResultCaveat | None:
    """Return the shared caveat for declared-only input scale evidence."""

    if (
        evidence.input_intensity_scale_evidence_level
        != IntensityScaleEvidenceLevel.DECLARED_BY_USER.value
    ):
        return None
    details = evidence.to_payload()
    details["workflow_scope"] = str(workflow_scope)
    if extra_details is not None:
        details.update({str(key): value for key, value in extra_details.items()})
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
    "InputIntensityScaleEvidence",
    "build_declared_input_intensity_scale_evidence_caveat",
    "input_intensity_scale_evidence_from_state",
]
