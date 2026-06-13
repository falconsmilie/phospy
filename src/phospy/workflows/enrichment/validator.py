"""Internal validator for enrichment workflow requests."""

from __future__ import annotations

from phospy.validation.workflows.enrichment import (
    EnrichmentWorkflowValidator as EnrichmentContractValidator,
)
from phospy.workflows.enrichment.models import ValidatedEnrichmentWorkflowRequest


class EnrichmentWorkflowValidator:
    """Validate `EnrichmentWorkflowRequest` before interpretation."""

    def __init__(
        self,
        *,
        contract_validator: EnrichmentContractValidator | None = None,
    ) -> None:
        self._contract_validator = contract_validator or EnrichmentContractValidator()

    def run(self, request: object) -> ValidatedEnrichmentWorkflowRequest:
        validated = self._contract_validator.run(request)
        return ValidatedEnrichmentWorkflowRequest(
            request=validated.request,
            identifier_column=validated.identifier_column,
            identifier_kind=validated.identifier_kind,
            set_collection=validated.set_collection,
            background_universe=validated.background_universe,
            selected_identifiers=validated.selected_identifiers,
            config=validated.config,
            selected_identifier_source=validated.selected_identifier_source,
        )


__all__ = ["EnrichmentWorkflowValidator"]
