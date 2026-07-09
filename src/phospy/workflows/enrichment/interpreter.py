"""Internal interpreter for enrichment workflow requests."""

from __future__ import annotations

from typing import cast

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.enrichment.models import ENRICHMENT_COLLECTION_KIND_GENE_SET
from phospy.science.enrichment.ora import (
    ORA_OUTSIDE_BACKGROUND_POLICY_DROP,
    ORA_STATISTICAL_TEST_HYPERGEOMETRIC,
    OraConfig,
)
from phospy.workflows.enrichment.models import (
    EnrichmentAnalysisLevel,
    EnrichmentIdentifierSemantics,
    InterpretedEnrichmentWorkflowRequest,
    ValidatedEnrichmentWorkflowRequest,
)


class EnrichmentWorkflowInterpreter:
    """Resolve a validated enrichment request into execution-ready inputs."""

    def run(
        self, request: ValidatedEnrichmentWorkflowRequest
    ) -> InterpretedEnrichmentWorkflowRequest:
        if not isinstance(cast(object, request), ValidatedEnrichmentWorkflowRequest):
            raise WorkflowBoundaryError(
                seam="enrichment.interpreter.validated_request_type",
                next_action=(
                    "pass validator output into EnrichmentWorkflowInterpreter.run"
                ),
                message_prefix="enrichment workflow boundary validation failed",
            )

        identifier_semantics = EnrichmentIdentifierSemantics(
            identifier_column=request.identifier_column,
            identifier_kind=request.identifier_kind,
            collection_kind=request.set_collection.collection_kind,
            analysis_level=_analysis_level_for_collection(
                request.set_collection.collection_kind
            ),
        )
        method_config = OraConfig(
            statistical_test=ORA_STATISTICAL_TEST_HYPERGEOMETRIC,
            selected_outside_background_policy=ORA_OUTSIDE_BACKGROUND_POLICY_DROP,
            set_outside_background_policy=ORA_OUTSIDE_BACKGROUND_POLICY_DROP,
            multiple_testing_correction=request.config.multiple_testing_correction,
        )
        method_metadata: dict[str, object] = {
            "method": request.config.method,
            "statistical_test": method_config.statistical_test,
            "multiple_testing_correction": method_config.multiple_testing_correction,
            "selected_outside_background_policy": (
                method_config.selected_outside_background_policy
            ),
            "set_outside_background_policy": method_config.set_outside_background_policy,
        }
        background_summary: dict[str, object] = {
            "source": "explicit",
            "provided_identifier_count": request.background_identifier_input_count,
            "universe_size": len(request.background_universe),
            "selected_identifier_count": len(request.selected_identifiers),
            "selected_identifier_input_count": request.selected_identifier_input_count,
            "selected_identifier_source": request.selected_identifier_source,
        }
        set_collection_summary = _summarise_set_collection(
            request=request,
        )
        diagnostics: dict[str, object] = {
            "interpreter": {
                "selected_identifiers_prepared": len(request.selected_identifiers),
                "background_universe_prepared": len(request.background_universe),
                "set_count_prepared": len(request.set_collection.enrichment_sets),
            }
        }
        return InterpretedEnrichmentWorkflowRequest(
            selected_identifiers=request.selected_identifiers,
            background_universe=request.background_universe,
            set_collection=request.set_collection,
            method_config=method_config,
            identifier_semantics=identifier_semantics,
            config=request.config,
            selected_identifier_source=request.selected_identifier_source,
            method_metadata=method_metadata,
            background_summary=background_summary,
            set_collection_summary=set_collection_summary,
            diagnostics=diagnostics,
            selected_identifier_input_count=request.selected_identifier_input_count,
            background_identifier_input_count=request.background_identifier_input_count,
        )


def _analysis_level_for_collection(
    collection_kind: str,
) -> EnrichmentAnalysisLevel:
    if collection_kind == ENRICHMENT_COLLECTION_KIND_GENE_SET:
        return "gene"
    return "ptm"


def _summarise_set_collection(
    *,
    request: ValidatedEnrichmentWorkflowRequest,
) -> dict[str, object]:
    enrichment_sets = request.set_collection.enrichment_sets
    all_members = tuple(
        identifier
        for enrichment_set in enrichment_sets
        for identifier in enrichment_set.identifiers
    )
    return {
        "collection_kind": request.set_collection.collection_kind,
        "identifier_kind": request.set_collection.identifier_kind,
        "set_count": len(enrichment_sets),
        "member_count": len(all_members),
        "distinct_member_count": len(frozenset(all_members)),
        "source_name": request.set_collection.source_name,
        "source_version": request.set_collection.source_version,
    }


__all__ = ["EnrichmentWorkflowInterpreter"]
