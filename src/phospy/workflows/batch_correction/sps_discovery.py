"""Public request and validation boundary for future SPS discovery execution."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from phospy.science.batch_correction.sps_discovery import (
    SpsDiscoveryConfig,
    SpsDiscoveryProvenance,
    SpsDiscoveryValidationError,
    SpsDiscoveryValidationIssue,
    SpsDiscoveryValidationResult,
    SpsReferenceDataset,
    SpsSelectionBoundaryCounts,
)
from phospy.validation.workflows.batch_correction.sps_discovery import (
    SpsDiscoveryRequestValidator,
)


@dataclass(frozen=True, slots=True, eq=False)
class SpsDiscoveryRequest:
    """Passive caller-intent payload for later SPS discovery execution."""

    reference_datasets: Sequence[SpsReferenceDataset]
    config: SpsDiscoveryConfig


class SpsDiscoveryWorkflow:
    """Validate SPS inputs and assemble provenance before future execution."""

    def validate(self, request: object) -> SpsDiscoveryValidationResult:
        """Return structured validation without executing ranking science."""

        if not isinstance(request, SpsDiscoveryRequest):
            return SpsDiscoveryValidationResult(
                issues=(
                    SpsDiscoveryValidationIssue(
                        code="invalid_discovery_request",
                        message="request must be an SpsDiscoveryRequest",
                        field_name="request",
                    ),
                ),
                reference_dataset_count=0,
                potential_overlap_site_count=0,
            )
        return SpsDiscoveryRequestValidator().run(
            reference_datasets=request.reference_datasets,
            config=request.config,
        )

    def require_valid(self, request: object) -> SpsDiscoveryValidationResult:
        """Validate a complete request or raise its structured domain error."""

        validation = self.validate(request)
        if not validation.valid:
            raise SpsDiscoveryValidationError(validation)
        return validation

    def assemble_provenance(
        self,
        request: SpsDiscoveryRequest,
        *,
        actual_control_count: int,
        selection_boundaries: SpsSelectionBoundaryCounts,
    ) -> SpsDiscoveryProvenance:
        """Build provenance only after the complete request passes validation."""

        if not isinstance(request, SpsDiscoveryRequest):
            self.require_valid(request)
            raise AssertionError("unreachable after structured request validation")
        try:
            references = tuple(request.reference_datasets)
        except Exception as exc:
            raise SpsDiscoveryValidationError(
                SpsDiscoveryValidationResult(
                    issues=(
                        SpsDiscoveryValidationIssue(
                            code="invalid_reference_collection",
                            message=(
                                "reference_datasets could not be read as a stable "
                                f"sequence: {exc}"
                            ),
                            field_name="reference_datasets",
                        ),
                    ),
                    reference_dataset_count=0,
                    potential_overlap_site_count=0,
                )
            ) from exc
        validation = SpsDiscoveryRequestValidator().run(
            reference_datasets=references,
            config=request.config,
        )
        if not validation.valid:
            raise SpsDiscoveryValidationError(validation)
        return SpsDiscoveryProvenance(
            source_datasets=tuple(
                reference._to_provenance() for reference in references
            ),
            config=request.config,
            requested_control_count=request.config.top_n,
            actual_control_count=actual_control_count,
            selection_boundaries=selection_boundaries,
        )


__all__ = ["SpsDiscoveryRequest", "SpsDiscoveryWorkflow"]
