from __future__ import annotations

from pathlib import Path

import pytest

from phospy._deprecations import PhosPyDeprecationWarning
from phospy.errors import PhosPyInputError
from phospy.io.readers.enrichment_sets import (
    load_enrichment_sets_csv,
    load_enrichment_sets_gmt,
    load_enrichment_sets_table,
    load_enrichment_sets_tsv,
    read_enrichment_sets_csv,
    read_enrichment_sets_gmt,
    read_enrichment_sets_table,
    read_enrichment_sets_tsv,
)


def _fixture_path(name: str) -> Path:
    return Path(__file__).resolve().parents[1] / "fixtures" / "enrichment" / name


def test_reader_enrichment_gmt_loading() -> None:
    collection = read_enrichment_sets_gmt(
        _fixture_path("gene_sets.gmt"),
        identifier_kind="gene_symbol",
        source_name="local_gmt",
        source_version="v1",
    )

    assert collection.identifier_kind == "gene_symbol"
    assert collection.collection_kind == "gene_set"
    assert collection.members_by_set_id["MAPK_PATHWAY"] == (
        "AKT1",
        "MAPK1",
        "MTOR",
    )
    assert collection.set_by_id["MAPK_PATHWAY"].name == "MAPK_PATHWAY"
    assert collection.set_by_id["MAPK_PATHWAY"].description == (
        "MAPK signaling pathway"
    )
    assert collection.set_by_id["MAPK_PATHWAY"].source_name == "local_gmt"
    assert collection.source_version == "v1"


def test_reader_enrichment_csv_loading_preserves_source_metadata() -> None:
    collection = read_enrichment_sets_table(_fixture_path("gene_sets.csv"))

    assert collection.identifier_kind == "gene_symbol"
    assert collection.source_name == "curated_unit"
    assert collection.source_version == "2026.06"
    assert collection.members_by_set_id["KINASE_RESPONSE"] == ("AKT1", "MAPK1")
    assert collection.set_by_id["KINASE_RESPONSE"].name == "Kinase response"
    assert collection.set_by_id["KINASE_RESPONSE"].description == "Response genes"
    assert collection.set_by_id["KINASE_RESPONSE"].source_name == "curated_unit"
    assert collection.set_by_id["CELL_CYCLE"].source_version == "2026.06"


def test_reader_enrichment_tsv_loading_requires_explicit_identifier_kind() -> None:
    path = _fixture_path("ptm_sets.tsv")

    with pytest.raises(PhosPyInputError, match="identifier_kind must be provided"):
        read_enrichment_sets_table(path)

    collection = read_enrichment_sets_table(path, identifier_kind="site_key")

    assert collection.identifier_kind == "site_key"
    assert collection.collection_kind == "ptm_set"
    assert collection.members_by_set_id["MOTIF_A"] == (
        "rat|P12345|S10",
        "rat|P12345|T20",
    )
    assert collection.set_by_id["MOTIF_A"].source_name == "local_ptm"


def test_reader_enrichment_mixed_identifier_kind_rejected(tmp_path: Path) -> None:
    path = tmp_path / "mixed.csv"
    path.write_text(
        "\n".join(
            (
                "set_id,name,identifier,identifier_kind",
                "GENES,Genes,AKT1,gene_symbol",
                "SITES,Sites,rat|P12345|S10,site_key",
            )
        ),
        encoding="utf-8",
    )

    with pytest.raises(PhosPyInputError, match="cannot mix identifier_kind"):
        read_enrichment_sets_table(path)


def test_reader_enrichment_empty_set_handling(tmp_path: Path) -> None:
    path = tmp_path / "empty.gmt"
    path.write_text("EMPTY\tEmpty set\n", encoding="utf-8")

    with pytest.raises(PhosPyInputError, match="identifiers must not be empty"):
        read_enrichment_sets_gmt(path, identifier_kind="gene_symbol")


@pytest.mark.parametrize(
    ("alias", "reader", "fixture_name", "kwargs", "replacement_name"),
    [
        (
            load_enrichment_sets_gmt,
            read_enrichment_sets_gmt,
            "gene_sets.gmt",
            {"identifier_kind": "gene_symbol"},
            "read_enrichment_sets_gmt",
        ),
        (
            load_enrichment_sets_table,
            read_enrichment_sets_table,
            "gene_sets.csv",
            {},
            "read_enrichment_sets_table",
        ),
        (
            load_enrichment_sets_csv,
            read_enrichment_sets_csv,
            "gene_sets.csv",
            {},
            "read_enrichment_sets_csv",
        ),
        (
            load_enrichment_sets_tsv,
            read_enrichment_sets_tsv,
            "ptm_sets.tsv",
            {"identifier_kind": "site_key"},
            "read_enrichment_sets_tsv",
        ),
    ],
)
def test_reader_enrichment_load_aliases_warn_and_forward(
    alias,
    reader,
    fixture_name: str,
    kwargs: dict[str, str],
    replacement_name: str,
) -> None:
    path = _fixture_path(fixture_name)
    expected = reader(path, **kwargs)

    with pytest.warns(
        PhosPyDeprecationWarning,
        match=rf"{alias.__name__}.*{replacement_name}",
    ):
        observed = alias(path, **kwargs)

    assert observed == expected
