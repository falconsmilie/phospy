"""Immutable planning boundary for validated SPS discovery requests."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from phospy.preprocessing.control_sites.discovery_contracts import (
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsReferenceDataset,
)
from phospy.provenance.models import TableFingerprint


def _empty_intensity_fingerprints() -> Mapping[str, TableFingerprint]:
    return MappingProxyType({})


@dataclass(frozen=True, slots=True)
class SpsDiscoveryPlan:
    """Execution-ready references and configuration after request validation."""

    reference_datasets: tuple[SpsReferenceDataset, ...]
    config: SpsDiscoveryConfig
    intensity_fingerprints_by_dataset: Mapping[str, TableFingerprint] = field(
        default_factory=_empty_intensity_fingerprints
    )


class SpsDiscoveryRequestInterpreter:
    """Resolve an accepted passive request without repeating selection science."""

    def __init__(
        self,
        *,
        intensity_fingerprints_by_dataset: Mapping[str, TableFingerprint] | None = None,
    ) -> None:
        self._intensity_fingerprints_by_dataset = MappingProxyType(
            dict(intensity_fingerprints_by_dataset or {})
        )

    def run(self, request: SpsDiscoveryRequest) -> SpsDiscoveryPlan:
        if not isinstance(request, SpsDiscoveryRequest):
            raise TypeError("request must be an SpsDiscoveryRequest")
        return SpsDiscoveryPlan(
            reference_datasets=tuple(request.reference_datasets),
            config=request.config,
            intensity_fingerprints_by_dataset=(self._intensity_fingerprints_by_dataset),
        )


__all__: list[str] = []
