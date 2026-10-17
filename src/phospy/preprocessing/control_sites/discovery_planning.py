"""Immutable planning boundary for validated SPS discovery requests."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.preprocessing.control_sites.discovery_contracts import (
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsReferenceDataset,
)


@dataclass(frozen=True, slots=True)
class SpsDiscoveryPlan:
    """Execution-ready references and configuration after request validation."""

    reference_datasets: tuple[SpsReferenceDataset, ...]
    config: SpsDiscoveryConfig


class SpsDiscoveryRequestInterpreter:
    """Resolve an accepted passive request without repeating selection science."""

    def run(self, request: SpsDiscoveryRequest) -> SpsDiscoveryPlan:
        if not isinstance(request, SpsDiscoveryRequest):
            raise TypeError("request must be an SpsDiscoveryRequest")
        return SpsDiscoveryPlan(
            reference_datasets=tuple(request.reference_datasets),
            config=request.config,
        )


__all__: list[str] = []
