"""Typed enrichment identifier-set provenance contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum

from phospy.provenance.models import InputIntensityScaleEvidence


class EnrichmentIdentifierSetSourceType(str, Enum):
    """Source category for an enrichment identifier set."""

    MANUAL = "manual"
    RAW_IDENTIFIER_LIST = "raw_identifier_list"
    PHOSPY_DERIVED_QUANTITATIVE = "phospy_derived_quantitative"


@dataclass(frozen=True, slots=True)
class EnrichmentIdentifierSetProvenance:
    """Typed provenance for selected or background enrichment identifiers."""

    source_type: EnrichmentIdentifierSetSourceType
    source_label: str
    identifier_count: int
    upstream_workflow_id: str | None = None
    upstream_result_id: str | None = None
    input_intensity_scale_evidence: InputIntensityScaleEvidence | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "source_type",
            _coerce_enrichment_identifier_set_source_type(self.source_type),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_evidence",
            _coerce_input_intensity_scale_evidence(
                self.input_intensity_scale_evidence,
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a compact JSON-compatible provenance payload."""

        payload: dict[str, object] = {
            "source_type": self.source_type.value,
            "source_label": self.source_label,
            "identifier_count": int(self.identifier_count),
        }
        if self.upstream_workflow_id is not None:
            payload["upstream_workflow_id"] = self.upstream_workflow_id
        if self.upstream_result_id is not None:
            payload["upstream_result_id"] = self.upstream_result_id
        if self.input_intensity_scale_evidence is not None:
            payload["input_intensity_scale_evidence"] = (
                self.input_intensity_scale_evidence.to_payload()
            )
        return payload


def _coerce_enrichment_identifier_set_source_type(
    value: object,
) -> EnrichmentIdentifierSetSourceType:
    if isinstance(value, EnrichmentIdentifierSetSourceType):
        return value
    try:
        return EnrichmentIdentifierSetSourceType(str(value).strip())
    except ValueError as exc:
        supported = ", ".join(item.value for item in EnrichmentIdentifierSetSourceType)
        raise ValueError(
            "enrichment identifier-set provenance source_type must be one of: "
            + supported
        ) from exc


def _coerce_input_intensity_scale_evidence(
    value: object | None,
) -> InputIntensityScaleEvidence | None:
    if value is None or isinstance(value, InputIntensityScaleEvidence):
        return value
    if isinstance(value, Mapping):
        return InputIntensityScaleEvidence(
            input_intensity_scale=value.get("input_intensity_scale"),
            input_intensity_scale_evidence_level=value.get(
                "input_intensity_scale_evidence_level"
            ),
            input_intensity_scale_source=value.get("input_intensity_scale_source"),
            input_intensity_scale_source_detail=value.get(
                "input_intensity_scale_source_detail"
            ),
        )
    raise TypeError(
        "input_intensity_scale_evidence must be InputIntensityScaleEvidence, "
        "mapping, or None"
    )


__all__ = [
    "EnrichmentIdentifierSetProvenance",
    "EnrichmentIdentifierSetSourceType",
]
