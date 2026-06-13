from __future__ import annotations

import pytest

from phospy.api import (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    EnrichmentSet,
    EnrichmentSetCollection,
    GeneSetCollection,
    WorkflowValidationError,
)


def test_enrichment_collection_constructs_in_memory() -> None:
    collection = EnrichmentSetCollection(
        sets=(
            EnrichmentSet(
                set_id="MAPK_PATHWAY",
                name="MAPK pathway",
                identifiers=("AKT1", "MAPK1"),
                identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
                source_name="curated",
                source_version="2026.06",
                description="Local curated gene set",
            ),
        ),
        source_name="curated",
        source_version="2026.06",
    )

    assert collection.identifier_kind == ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    assert collection.collection_kind == "gene_set"
    assert collection.set_ids == ("MAPK_PATHWAY",)
    assert collection.members_by_set_id == {"MAPK_PATHWAY": ("AKT1", "MAPK1")}
    assert collection.term_names == {"MAPK_PATHWAY": "MAPK pathway"}
    assert collection.set_by_id["MAPK_PATHWAY"].source_name == "curated"
    assert collection.set_by_id["MAPK_PATHWAY"].source_version == "2026.06"
    assert collection.set_by_id["MAPK_PATHWAY"].description == (
        "Local curated gene set"
    )


def test_enrichment_collection_duplicate_identifiers_within_set_are_deduplicated() -> (
    None
):
    enrichment_set = EnrichmentSet(
        set_id="DUPLICATES",
        name="Duplicate members",
        identifiers=(" AKT1 ", "MAPK1", "AKT1", "MTOR", "MAPK1"),
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    )

    assert enrichment_set.identifiers == ("AKT1", "MAPK1", "MTOR")


def test_enrichment_collection_mixed_identifier_kind_rejected() -> None:
    with pytest.raises(WorkflowValidationError, match="cannot mix identifier_kind"):
        EnrichmentSetCollection(
            sets=(
                EnrichmentSet(
                    set_id="GENES",
                    name="Genes",
                    identifiers=("AKT1",),
                    identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
                ),
                EnrichmentSet(
                    set_id="SITES",
                    name="Sites",
                    identifiers=("rat|P12345|S10",),
                    identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
                ),
            )
        )


def test_enrichment_collection_empty_set_rejected() -> None:
    with pytest.raises(WorkflowValidationError, match="identifiers must not be empty"):
        EnrichmentSet(
            set_id="EMPTY",
            name="Empty set",
            identifiers=(),
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        )


def test_enrichment_collection_legacy_gene_wrapper_preserves_mapping_view() -> None:
    collection = GeneSetCollection(
        sets={"MAPK_PATHWAY": ("AKT1", "AKT1", "MAPK1")},
        term_names={"MAPK_PATHWAY": "MAPK pathway"},
        source_name="legacy_user",
    )

    assert collection.sets == {"MAPK_PATHWAY": ("AKT1", "MAPK1")}
    assert collection.term_names == {"MAPK_PATHWAY": "MAPK pathway"}
    assert collection.set_by_id["MAPK_PATHWAY"].source_name == "legacy_user"
