"""Structured caveats for enrichment workflow results."""

from __future__ import annotations

from phospy.contracts.result_caveats import ResultCaveat
from phospy.science.enrichment.models import ENRICHMENT_METHOD_OVER_REPRESENTATION
from phospy.workflows.enrichment.models import InterpretedEnrichmentWorkflowRequest
from phospy.workflows.result_caveat_helpers import deduplicate_caveats

ENRICHMENT_OFFLINE_ORA_ONLY_SCOPE_CAVEAT_CODE = "enrichment_offline_ora_only_scope"
ENRICHMENT_BACKGROUND_UNIVERSE_ASSUMPTION_CAVEAT_CODE = (
    "enrichment_background_universe_assumption"
)
ENRICHMENT_IDENTIFIER_KIND_ASSUMPTION_CAVEAT_CODE = (
    "enrichment_identifier_kind_assumption"
)
ENRICHMENT_NO_RANK_BASED_SEMANTICS_CAVEAT_CODE = "enrichment_no_rank_based_semantics"


def build_enrichment_result_caveats(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    background_summary: dict[str, object],
    set_collection_summary: dict[str, object],
) -> tuple[ResultCaveat, ...]:
    """Build compact machine-readable caveats for enrichment workflow results."""

    return deduplicate_caveats(
        (
            _offline_ora_only_scope_caveat(
                request=request,
                set_collection_summary=set_collection_summary,
            ),
            _background_universe_assumption_caveat(
                request=request,
                background_summary=background_summary,
            ),
            _identifier_kind_assumption_caveat(request),
            _no_rank_based_semantics_caveat(request),
        )
    )


def _offline_ora_only_scope_caveat(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    set_collection_summary: dict[str, object],
) -> ResultCaveat:
    return ResultCaveat(
        code=ENRICHMENT_OFFLINE_ORA_ONLY_SCOPE_CAVEAT_CODE,
        severity="info",
        message=(
            "Enrichment workflow runs offline over-representation analysis against "
            "the caller-supplied set collection only; pathway database resources "
            "are not bundled, fetched, or expanded."
        ),
        details={
            "method": ENRICHMENT_METHOD_OVER_REPRESENTATION,
            "workflow_scope": "offline_ora",
            "online_resources_used": False,
            "set_collection_scope": "caller_supplied_only",
            "pathway_database_coverage": "caller_supplied_sets_only",
            "collection_kind": request.set_collection.collection_kind,
            "set_count": set_collection_summary.get(
                "set_count",
                len(request.set_collection.enrichment_sets),
            ),
            "source_name": request.set_collection.source_name,
            "source_version": request.set_collection.source_version,
        },
    )


def _background_universe_assumption_caveat(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    background_summary: dict[str, object],
) -> ResultCaveat:
    return ResultCaveat(
        code=ENRICHMENT_BACKGROUND_UNIVERSE_ASSUMPTION_CAVEAT_CODE,
        severity="warning",
        message=(
            "ORA denominators and foreground/set filtering depend on the "
            "caller-supplied background_universe; the workflow does not infer a "
            "detectable universe."
        ),
        details={
            "background_universe_source": "caller_supplied",
            "background_universe_required": True,
            "background_universe_inferred": False,
            "background_universe_size": len(request.background_universe),
            "selected_identifier_count": len(request.selected_identifiers),
            "selected_in_background_count": background_summary.get(
                "selected_in_background_count"
            ),
            "dropped_selected_count": background_summary.get("dropped_selected_count"),
            "foreground_identifiers_missing_from_background_count": (
                background_summary.get(
                    "foreground_identifiers_missing_from_background_count"
                )
            ),
        },
    )


def _identifier_kind_assumption_caveat(
    request: InterpretedEnrichmentWorkflowRequest,
) -> ResultCaveat:
    semantics = request.identifier_semantics
    return ResultCaveat(
        code=ENRICHMENT_IDENTIFIER_KIND_ASSUMPTION_CAVEAT_CODE,
        severity="info",
        message=(
            "Enrichment identifiers are treated as the caller-declared "
            "identifier_kind; the workflow does not map identifiers or verify "
            "equivalence across identifier namespaces."
        ),
        details={
            "identifier_column": semantics.identifier_column,
            "identifier_kind": semantics.identifier_kind,
            "collection_identifier_kind": request.set_collection.identifier_kind,
            "collection_kind": semantics.collection_kind,
            "analysis_level": semantics.analysis_level,
            "identifier_mapping_performed": False,
        },
    )


def _no_rank_based_semantics_caveat(
    request: InterpretedEnrichmentWorkflowRequest,
) -> ResultCaveat:
    return ResultCaveat(
        code=ENRICHMENT_NO_RANK_BASED_SEMANTICS_CAVEAT_CODE,
        severity="info",
        message=(
            "Enrichment results contain ORA overlap counts and p-values only; "
            "GSEA, ssGSEA, PTM-SEA, ranking, leading-edge, and enrichment-score "
            "semantics are not computed by this workflow."
        ),
        details={
            "computed_method": ENRICHMENT_METHOD_OVER_REPRESENTATION,
            "ranked_input_consumed": False,
            "rank_based_enrichment_supported": False,
            "unsupported_rank_based_methods": ("GSEA", "ssGSEA", "PTM-SEA"),
            "analysis_level": request.identifier_semantics.analysis_level,
        },
    )


__all__ = [
    "ENRICHMENT_BACKGROUND_UNIVERSE_ASSUMPTION_CAVEAT_CODE",
    "ENRICHMENT_IDENTIFIER_KIND_ASSUMPTION_CAVEAT_CODE",
    "ENRICHMENT_NO_RANK_BASED_SEMANTICS_CAVEAT_CODE",
    "ENRICHMENT_OFFLINE_ORA_ONLY_SCOPE_CAVEAT_CODE",
    "build_enrichment_result_caveats",
]
