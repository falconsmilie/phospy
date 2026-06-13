from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from phospy.science.datasets.preprocessing.protein_mapping import (
    ProteinMappingConfig,
    ProteinMappingResolver,
    ProteinMappingStatus,
)


def _resolve(
    *,
    site_metadata: pd.DataFrame,
    phospho_index: list[str],
    total_index: list[str],
    config: ProteinMappingConfig | None = None,
):
    return ProteinMappingResolver().run(
        site_metadata=site_metadata,
        phospho_matrix_index=pd.Index(phospho_index, name="site_id"),
        total_protein_matrix_index=pd.Index(total_index, name="protein_id"),
        config=config,
    )


def test_protein_mapping_exact_one_site_to_one_protein_mapping() -> None:
    site_metadata = pd.DataFrame(
        {"protein_accession": ["P53778"]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    original_metadata = site_metadata.copy(deep=True)

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["MAPK14;Y182;"],
        total_index=["P53778"],
    )

    assert result.status_by_site == {"MAPK14;Y182;": ProteinMappingStatus.MATCHED}
    assert result.site_to_protein_identifier == {"MAPK14;Y182;": "P53778"}
    assert result.site_to_total_protein_row_key == {"MAPK14;Y182;": "P53778"}
    assert result.records[0].status is ProteinMappingStatus.MATCHED
    pdt.assert_frame_equal(site_metadata, original_metadata)


def test_protein_mapping_uses_protein_identifier_from_site_metadata() -> None:
    site_metadata = pd.DataFrame(
        {
            "protein_accession": ["P31749"],
            "display_id": ["AKT1;T308;"],
        },
        index=pd.Index(["AKT1;T308;"], name="site_id"),
    )

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["AKT1;T308;"],
        total_index=["AKT1;T308;", "P31749"],
        config=ProteinMappingConfig(allow_display_label_fallback=True),
    )

    record = result.records[0]
    assert record.status is ProteinMappingStatus.MATCHED
    assert record.protein_identifier == "P31749"
    assert record.total_protein_row_key == "P31749"
    assert record.protein_identifier_source == "protein_identifier:protein_accession"


def test_protein_mapping_missing_protein_identifier_status() -> None:
    site_metadata = pd.DataFrame(
        {"protein_accession": [""]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["MAPK14;Y182;"],
        total_index=["MAPK14;Y182;"],
    )

    record = result.records[0]
    assert record.status is ProteinMappingStatus.MISSING_SITE_PROTEIN_IDENTIFIER
    assert record.protein_identifier is None
    assert result.site_to_protein_identifier == {}
    assert result.site_to_total_protein_row_key == {}


def test_protein_mapping_missing_total_protein_row_status() -> None:
    site_metadata = pd.DataFrame(
        {"protein_accession": ["P53778"]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["MAPK14;Y182;"],
        total_index=["P31749"],
    )

    record = result.records[0]
    assert record.status is ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW
    assert record.protein_identifier == "P53778"
    assert record.total_protein_row_key is None
    assert result.site_to_protein_identifier == {"MAPK14;Y182;": "P53778"}
    assert result.site_to_total_protein_row_key == {}


def test_protein_mapping_ambiguous_site_to_protein_mapping_status() -> None:
    site_metadata = pd.DataFrame(
        {"protein_accession": ["P53778", "Q16539"]},
        index=pd.Index(["MAPK14;Y182;", "MAPK14;Y182;"], name="site_id"),
    )

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["MAPK14;Y182;"],
        total_index=["P53778", "Q16539"],
    )

    record = result.records[0]
    assert record.status is ProteinMappingStatus.AMBIGUOUS_SITE_PROTEIN_MAPPING
    assert record.protein_identifier is None
    assert record.candidate_protein_identifiers == ("P53778", "Q16539")
    assert result.site_to_protein_identifier == {}
    assert result.site_to_total_protein_row_key == {}


def test_protein_mapping_ambiguous_total_protein_row_mapping_status() -> None:
    site_metadata = pd.DataFrame(
        {"protein_accession": ["P53778"]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["MAPK14;Y182;"],
        total_index=["P53778", "P53778"],
    )

    record = result.records[0]
    assert record.status is ProteinMappingStatus.AMBIGUOUS_TOTAL_PROTEIN_MAPPING
    assert record.protein_identifier == "P53778"
    assert record.total_protein_row_key is None
    assert record.candidate_total_protein_row_keys == ("P53778", "P53778")
    assert result.site_to_protein_identifier == {"MAPK14;Y182;": "P53778"}
    assert result.site_to_total_protein_row_key == {}


def test_protein_mapping_output_ordering_is_deterministic() -> None:
    site_metadata = pd.DataFrame(
        {
            "protein_accession": ["P_C", "P_B", "P_A"],
        },
        index=pd.Index(["site_c", "site_b", "site_a"], name="site_id"),
    )

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["site_b", "site_a", "site_c"],
        total_index=["P_A", "P_B", "P_C"],
    )

    assert tuple(record.site_key for record in result.records) == (
        "site_b",
        "site_a",
        "site_c",
    )
    assert list(result.site_to_total_protein_row_key) == [
        "site_b",
        "site_a",
        "site_c",
    ]


def test_protein_mapping_does_not_use_gene_symbol_without_policy() -> None:
    site_metadata = pd.DataFrame(
        {"gene_symbol": ["MAPK14"]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["MAPK14;Y182;"],
        total_index=["MAPK14"],
    )

    assert (
        result.records[0].status is ProteinMappingStatus.MISSING_SITE_PROTEIN_IDENTIFIER
    )


def test_protein_mapping_uses_gene_symbol_when_policy_permits_it() -> None:
    site_metadata = pd.DataFrame(
        {"gene_symbol": ["MAPK14"]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    result = _resolve(
        site_metadata=site_metadata,
        phospho_index=["MAPK14;Y182;"],
        total_index=["MAPK14"],
        config=ProteinMappingConfig(allow_gene_symbol_matching=True),
    )

    record = result.records[0]
    assert record.status is ProteinMappingStatus.MATCHED
    assert record.protein_identifier == "MAPK14"
    assert record.protein_identifier_source == "gene_symbol:gene_symbol"
