"""Request and orchestration boundary for SPS discovery."""

from __future__ import annotations

from phospy.preprocessing.control_sites.discovery_contracts import (
    SpsDiscoveryRequest,
    SpsDiscoveryValidationError,
    SpsDiscoveryValidationIssue,
    SpsDiscoveryValidationResult,
)
from phospy.preprocessing.control_sites.discovery_execution import (
    SpsDiscoveryExecutor,
)
from phospy.preprocessing.control_sites.discovery_planning import (
    SpsDiscoveryRequestInterpreter,
)
from phospy.preprocessing.control_sites.discovery_provenance import (
    SpsDiscoveryProvenance,
    SpsSelectionBoundaryCounts,
    reference_provenance,
)
from phospy.preprocessing.control_sites.discovery_results import (
    SpsDiscoveryResult,
)
from phospy.preprocessing.control_sites.discovery_validation import (
    SpsDiscoveryRequestValidator,
)


class SpsDiscoveryWorkflow:
    """Validate explicit reference inputs and execute SPS discovery."""

    def run(self, request: object) -> SpsDiscoveryResult:
        """Validate references and return their deterministic consensus ranking."""

        if not isinstance(request, SpsDiscoveryRequest):
            self.require_valid(request)
            raise AssertionError("unreachable after structured request validation")
        validator = SpsDiscoveryRequestValidator()
        validation = validator.run(
            reference_datasets=request.reference_datasets,
            config=request.config,
        )
        if not validation.valid:
            raise SpsDiscoveryValidationError(validation)
        plan = SpsDiscoveryRequestInterpreter(
            intensity_fingerprints_by_dataset=(
                validator.intensity_fingerprints_by_dataset
            ),
        ).run(request)
        return SpsDiscoveryExecutor().run(plan)

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
        validator = SpsDiscoveryRequestValidator()
        validation = validator.run(
            reference_datasets=references,
            config=request.config,
        )
        if not validation.valid:
            raise SpsDiscoveryValidationError(validation)
        ordered_references = tuple(
            sorted(references, key=lambda reference: reference.dataset_id)
        )
        return SpsDiscoveryProvenance(
            source_datasets=tuple(
                reference_provenance(
                    reference,
                    intensity_fingerprint=(
                        validator.intensity_fingerprints_by_dataset.get(
                            reference.dataset_id
                        )
                    ),
                )
                for reference in ordered_references
            ),
            config=request.config,
            requested_control_count=request.config.top_n,
            actual_control_count=actual_control_count,
            selection_boundaries=selection_boundaries,
        )


__all__ = ["SpsDiscoveryRequest", "SpsDiscoveryWorkflow"]
