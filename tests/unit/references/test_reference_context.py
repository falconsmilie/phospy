from __future__ import annotations

import pytest

from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.models import ReferenceProvenance
from phospy.science.references.models import (
    ReferenceContext,
    reference_context_from_provenance,
)


def _context(**overrides: object) -> ReferenceContext:
    values = {
        "organism": "rat",
        "protein_namespace": "display_site_id_gene_symbol_residue",
        "source_name": "unit reference",
        "source_version": "2026-06",
        "proteome_version": None,
        "reference_table_sha256": "a" * 64,
    }
    values.update(overrides)
    return ReferenceContext(**values)


def test_reference_context_id_is_deterministic() -> None:
    first = _context()
    second = _context()

    assert first.reference_context_id == second.reference_context_id
    assert first.reference_context_id.startswith("reference-context-v1:")
    assert len(first.reference_context_id.removeprefix("reference-context-v1:")) == 64


def test_different_source_versions_produce_different_ids() -> None:
    first = _context(source_version="2026-06")
    second = _context(source_version="2026-07")

    assert first.reference_context_id != second.reference_context_id
    assert first != second


def test_different_protein_namespaces_produce_different_ids() -> None:
    first = _context(protein_namespace="gene_symbol")
    second = _context(protein_namespace="uniprot_accession")

    assert first.reference_context_id != second.reference_context_id
    assert first != second


def test_same_fields_produce_same_id_and_compare_equal() -> None:
    first = _context(
        organism=" RAT ",
        reference_table_sha256=("A" * 64),
    )
    second = _context(
        organism="rat",
        reference_table_sha256=("a" * 64),
    )

    assert first.reference_context_id == second.reference_context_id
    assert first == second


def test_invalid_empty_organism_rejected() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="reference_context.organism must be non-empty",
    ):
        _context(organism=" ")


def test_invalid_empty_protein_namespace_rejected() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="reference_context.protein_namespace must be non-empty",
    ):
        _context(protein_namespace="")


def test_reference_context_can_be_resolved_from_complete_provenance() -> None:
    provenance = ReferenceProvenance(
        source_type="explicit",
        organism="rat",
        bundle_id=None,
        source_name="unit reference",
        source_version="2026-06",
        identifier_namespace="gene_symbol",
        table_fingerprints=(),
    )

    context = reference_context_from_provenance(provenance)

    assert context is not None
    assert context.organism == "rat"
    assert context.protein_namespace == "gene_symbol"
