from __future__ import annotations

from phospy.api import (
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
    ResultCaveat,
)
from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
from phospy.workflows.enrichment.caveats import (
    ENRICHMENT_BACKGROUND_UNIVERSE_ASSUMPTION_CAVEAT_CODE,
    ENRICHMENT_IDENTIFIER_KIND_ASSUMPTION_CAVEAT_CODE,
    ENRICHMENT_NO_RANK_BASED_SEMANTICS_CAVEAT_CODE,
    ENRICHMENT_OFFLINE_ORA_ONLY_SCOPE_CAVEAT_CODE,
)


def _request() -> EnrichmentWorkflowRequest:
    return EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={
                "MAPK_PATHWAY": ("AKT1", "MAPK1"),
                "MTOR_SIGNALING": ("MTOR",),
            },
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            term_names={
                "MAPK_PATHWAY": "MAPK pathway",
                "MTOR_SIGNALING": "MTOR signaling",
            },
            source_name="unit_test",
        ),
        selected_identifiers=("AKT1", "MAPK1", "OUTSIDE_BACKGROUND"),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(selected_outside_background_policy="drop"),
    )


def _caveat_by_code() -> dict[str, ResultCaveat]:
    result = EnrichmentWorkflow().run(_request())
    return {caveat.code: caveat for caveat in result.caveats}


def test_enrichment_result_caveats_include_ora_only_scope() -> None:
    caveats = _caveat_by_code()

    caveat = caveats[ENRICHMENT_OFFLINE_ORA_ONLY_SCOPE_CAVEAT_CODE]

    assert caveat.severity == "info"
    assert caveat.details["workflow_scope"] == "offline_ora"
    assert caveat.details["online_resources_used"] is False
    assert caveat.details["set_collection_scope"] == "caller_supplied_only"
    assert caveat.details["pathway_database_coverage"] == ("caller_supplied_sets_only")


def test_enrichment_result_caveats_include_background_universe_assumption() -> None:
    caveats = _caveat_by_code()

    caveat = caveats[ENRICHMENT_BACKGROUND_UNIVERSE_ASSUMPTION_CAVEAT_CODE]

    assert caveat.severity == "warning"
    assert caveat.details["background_universe_source"] == "caller_supplied"
    assert caveat.details["background_universe_required"] is True
    assert caveat.details["background_universe_inferred"] is False
    assert caveat.details["background_universe_size"] == 3
    assert caveat.details["foreground_identifiers_missing_from_background_count"] == 1


def test_enrichment_result_caveats_include_identifier_kind_assumption() -> None:
    caveats = _caveat_by_code()

    caveat = caveats[ENRICHMENT_IDENTIFIER_KIND_ASSUMPTION_CAVEAT_CODE]

    assert caveat.details["identifier_column"] == "gene_symbol"
    assert caveat.details["identifier_kind"] == ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    assert caveat.details["collection_identifier_kind"] == (
        ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    )
    assert caveat.details["identifier_mapping_performed"] is False


def test_enrichment_result_caveats_distinguish_ora_from_rank_based_methods() -> None:
    caveats = _caveat_by_code()

    caveat = caveats[ENRICHMENT_NO_RANK_BASED_SEMANTICS_CAVEAT_CODE]

    assert caveat.details["computed_method"] == "over_representation"
    assert caveat.details["ranked_input_consumed"] is False
    assert caveat.details["rank_based_enrichment_supported"] is False
    assert caveat.details["unsupported_rank_based_methods"] == (
        "GSEA",
        "ssGSEA",
        "PTM-SEA",
    )
