"""Internal validator for enrichment workflow requests."""

from __future__ import annotations

from collections.abc import Sequence

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
        selected_identifier_input_count = _selected_identifier_input_count(
            request=validated.request,
            identifier_column=validated.identifier_column,
        )
        background_identifier_input_count = _sequence_input_count(
            validated.request.background_universe,
            fallback_count=len(validated.background_universe),
        )
        return ValidatedEnrichmentWorkflowRequest(
            request=validated.request,
            identifier_column=validated.identifier_column,
            identifier_kind=validated.identifier_kind,
            set_collection=validated.set_collection,
            background_universe=validated.background_universe,
            selected_identifiers=validated.selected_identifiers,
            config=validated.config,
            selected_identifier_source=validated.selected_identifier_source,
            selected_identifier_input_count=selected_identifier_input_count,
            background_identifier_input_count=background_identifier_input_count,
        )


def _selected_identifier_input_count(
    *,
    request: object,
    identifier_column: str,
) -> int:
    selected_identifiers = getattr(request, "selected_identifiers", None)
    if selected_identifiers is not None:
        return _sequence_input_count(selected_identifiers, fallback_count=0)
    input_table = getattr(request, "input_table", None)
    if input_table is not None and hasattr(input_table, "shape"):
        return int(input_table.shape[0])
    _ = identifier_column
    return 0


def _sequence_input_count(value: object, *, fallback_count: int) -> int:
    if isinstance(value, str | bytes | bytearray):
        return int(fallback_count)
    if isinstance(value, Sequence):
        return int(len(value))
    return int(fallback_count)


__all__ = ["EnrichmentWorkflowValidator"]
