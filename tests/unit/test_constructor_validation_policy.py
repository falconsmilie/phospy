from __future__ import annotations

import pytest

from phospy.api import (
    ContractValidationError,
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
    EnrichmentSet,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
    KinaseScoringConfig,
    WorkflowValidationError,
)
from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
from phospy.api.results import ResultCaveat


def test_contract_value_object_constructor_failures_use_contract_error() -> None:
    with pytest.raises(ContractValidationError, match="scoring_config.min_substrates"):
        KinaseScoringConfig(reliability_profile="custom", min_substrates=1)

    with pytest.raises(ContractValidationError, match="result_caveat.severity"):
        ResultCaveat(code="low_support", severity="bad", message="Low support.")

    with pytest.raises(ContractValidationError, match="identifiers must not be empty"):
        EnrichmentSet(
            set_id="EMPTY",
            name="Empty set",
            identifiers=(),
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        )


def test_request_contextual_invalidity_reaches_workflow_validation() -> None:
    provenance = EnrichmentIdentifierSetProvenance(
        source_type=EnrichmentIdentifierSetSourceType.MANUAL,
        source_label="manual list",
        identifier_count=1,
    )
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={"mapk_pathway": ("AKT1", "MAPK1")},
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        ),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        selected_identifier_provenance=provenance,
    )

    with pytest.raises(
        WorkflowValidationError,
        match="Selected identifier-set provenance count mismatch",
    ):
        EnrichmentWorkflow().run(request)
