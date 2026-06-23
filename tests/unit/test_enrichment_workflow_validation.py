from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    EnrichmentConfig,
    EnrichmentSetCollection,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
    PtmSetCollection,
    WorkflowValidationError,
)
from phospy.validation.workflows.enrichment import EnrichmentWorkflowValidator


def _gene_collection() -> GeneSetCollection:
    return GeneSetCollection(
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
    )


def _site_collection() -> PtmSetCollection:
    return PtmSetCollection(
        sets={
            "MOTIF_SITES": (
                "rat|P12345|S10",
                "rat|P12345|T20",
            )
        },
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        source_name="unit_test",
    )


def _valid_gene_request() -> EnrichmentWorkflowRequest:
    return EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(),
    )


def test_enrichment_validation_accepts_valid_gene_level_request() -> None:
    validated = EnrichmentWorkflowValidator().run(_valid_gene_request())

    assert validated.identifier_kind == ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    assert validated.selected_identifier_source == "selected_identifiers"
    assert validated.selected_identifiers == ("AKT1", "MAPK1")
    assert validated.background_universe == ("AKT1", "MAPK1", "MTOR")
    assert validated.set_collection.collection_kind == "gene_set"


def test_enrichment_validation_accepts_valid_site_level_request() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="site_key",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        set_collection=_site_collection(),
        input_table=pd.DataFrame(
            {
                "site_key": [
                    "rat|P12345|S10",
                    " rat|P12345|T20 ",
                ]
            }
        ),
        background_universe=(
            "rat|P12345|S10",
            "rat|P12345|T20",
            "rat|P12345|Y30",
        ),
    )

    validated = EnrichmentWorkflowValidator().run(request)

    assert validated.identifier_kind == ENRICHMENT_IDENTIFIER_KIND_SITE_KEY
    assert validated.selected_identifier_source == "input_table"
    assert validated.selected_identifiers == (
        "rat|P12345|S10",
        "rat|P12345|T20",
    )
    assert validated.set_collection.collection_kind == "ptm_set"


def test_enrichment_validation_rejects_missing_background() -> None:
    request = _valid_gene_request()
    object.__setattr__(request, "background_universe", ())

    with pytest.raises(WorkflowValidationError, match="background_universe"):
        EnrichmentWorkflowValidator().run(request)


def test_enrichment_validation_rejects_identifier_kind_mismatch() -> None:
    request = _valid_gene_request()
    object.__setattr__(request, "identifier_kind", ENRICHMENT_IDENTIFIER_KIND_SITE_KEY)

    with pytest.raises(WorkflowValidationError, match="identifier_kind.*match"):
        EnrichmentWorkflowValidator().run(request)


def test_enrichment_validation_rejects_missing_identifier_column() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        input_table=pd.DataFrame({"protein_id": ["P31749"]}),
        background_universe=("AKT1", "MAPK1", "MTOR"),
    )

    with pytest.raises(WorkflowValidationError, match="missing required columns"):
        EnrichmentWorkflowValidator().run(request)


def test_enrichment_validation_rejects_empty_set_collection() -> None:
    empty_collection = object.__new__(EnrichmentSetCollection)
    object.__setattr__(empty_collection, "enrichment_sets", ())
    object.__setattr__(
        empty_collection,
        "identifier_kind",
        ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    )
    object.__setattr__(empty_collection, "collection_kind", "gene_set")
    object.__setattr__(empty_collection, "source_name", "unit_test")
    object.__setattr__(empty_collection, "source_version", None)
    request = _valid_gene_request()
    object.__setattr__(request, "set_collection", empty_collection)

    with pytest.raises(WorkflowValidationError, match="at least one set"):
        EnrichmentWorkflowValidator().run(request)


def test_enrichment_validation_rejects_unsupported_method() -> None:
    request = _valid_gene_request()
    config = EnrichmentConfig()
    object.__setattr__(config, "method", "competitive")
    object.__setattr__(request, "config", config)

    with pytest.raises(WorkflowValidationError, match="config.method"):
        EnrichmentWorkflowValidator().run(request)


def test_enrichment_validation_rejects_unsupported_correction_method() -> None:
    request = _valid_gene_request()
    config = EnrichmentConfig()
    object.__setattr__(config, "multiple_testing_correction", "storey")
    object.__setattr__(request, "config", config)

    with pytest.raises(
        WorkflowValidationError,
        match="config.multiple_testing_correction",
    ):
        EnrichmentWorkflowValidator().run(request)


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        ({"min_set_size": 0}, "enrichment.min_set_size"),
        ({"max_set_size": 0}, "enrichment.max_set_size"),
        ({"min_set_size": True}, "enrichment.min_set_size"),
        ({"min_set_size": 4, "max_set_size": 3}, "min_set_size"),
    ],
)
def test_enrichment_config_rejects_invalid_set_size_filters(
    kwargs: dict[str, object],
    pattern: str,
) -> None:
    with pytest.raises(WorkflowValidationError, match=pattern):
        EnrichmentConfig(**kwargs)  # type: ignore[arg-type]


def test_enrichment_validation_rejects_mutated_invalid_set_size_filters() -> None:
    request = _valid_gene_request()
    config = EnrichmentConfig()
    object.__setattr__(config, "min_set_size", 4)
    object.__setattr__(config, "max_set_size", 3)
    object.__setattr__(request, "config", config)

    with pytest.raises(WorkflowValidationError, match="config.min_set_size"):
        EnrichmentWorkflowValidator().run(request)


def test_enrichment_validation_preserves_selected_identifiers_outside_background() -> (
    None
):
    request = _valid_gene_request()
    object.__setattr__(request, "selected_identifiers", ("AKT1", "UNKNOWN"))

    validated = EnrichmentWorkflowValidator().run(request)

    assert validated.selected_identifiers == ("AKT1", "UNKNOWN")
    assert validated.background_universe == ("AKT1", "MAPK1", "MTOR")


def test_enrichment_validation_does_not_run_ora(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.science.enrichment.ora as ora

    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("ORA must not run during enrichment validation")

    monkeypatch.setattr(ora, "run", fail_if_called)

    validated = EnrichmentWorkflowValidator().run(_valid_gene_request())

    assert validated.selected_identifiers == ("AKT1", "MAPK1")
