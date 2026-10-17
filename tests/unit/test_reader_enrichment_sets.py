from __future__ import annotations

from pathlib import Path

import pytest

from phospy.errors import PhosPyInputError
from phospy.io.readers.enrichment_sets import (
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


def test_reader_enrichment_malformed_gmt_error_is_exact(tmp_path: Path) -> None:
    path = tmp_path / "malformed.gmt"
    path.write_text("ONLY_A_SET_ID\n", encoding="utf-8")

    with pytest.raises(PhosPyInputError) as exc_info:
        read_enrichment_sets_gmt(path, identifier_kind="gene_symbol")

    assert str(exc_info.value) == (
        f"failed to parse GMT enrichment set file '{path}': line 1 must contain "
        "set_id and description columns"
    )


def test_reader_enrichment_malformed_csv_error_is_exact(tmp_path: Path) -> None:
    path = tmp_path / "malformed.csv"
    path.write_text("set_id,name\nA,Alpha\n", encoding="utf-8")

    with pytest.raises(PhosPyInputError) as exc_info:
        read_enrichment_sets_csv(path, identifier_kind="gene_symbol")

    assert str(exc_info.value) == (
        f"enrichment set table '{path}' is missing required columns: identifier"
    )


def test_reader_enrichment_malformed_tsv_error_is_exact(tmp_path: Path) -> None:
    path = tmp_path / "malformed.tsv"
    path.write_text(
        "set_id\tname\tidentifier\tidentifier_kind\nMOTIF_A\tMotif A\t\tsite_key\n",
        encoding="utf-8",
    )

    with pytest.raises(PhosPyInputError) as exc_info:
        read_enrichment_sets_tsv(path)

    assert str(exc_info.value) == (
        f"enrichment set table '{path}' row 2 column 'identifier' must be non-empty"
    )


def test_reader_enrichment_csv_and_tsv_preserve_metadata_precedence(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "metadata.csv"
    csv_path.write_text(
        "set_id,name,identifier,identifier_kind,source_name,source_version\n"
        "A,Alpha,AKT1,gene_symbol,file_source,file-v1\n",
        encoding="utf-8",
    )
    tsv_path = tmp_path / "metadata.tsv"
    tsv_path.write_text(
        "set_id\tname\tidentifier\tidentifier_kind\tsource_name\tsource_version\n"
        "M\tMotif\trat|P1|S1\tsite_key\t\t\n",
        encoding="utf-8",
    )

    csv_collection = read_enrichment_sets_csv(
        csv_path,
        source_name="caller_source",
        source_version="caller-v1",
    )
    tsv_collection = read_enrichment_sets_tsv(
        tsv_path,
        source_name="caller_source",
        source_version="caller-v1",
    )

    assert csv_collection.source_name == "caller_source"
    assert csv_collection.source_version == "caller-v1"
    assert csv_collection.set_by_id["A"].source_name == "file_source"
    assert csv_collection.set_by_id["A"].source_version == "file-v1"
    assert tsv_collection.source_name == "caller_source"
    assert tsv_collection.source_version == "caller-v1"
    assert tsv_collection.set_by_id["M"].source_name == "caller_source"
    assert tsv_collection.set_by_id["M"].source_version == "caller-v1"
