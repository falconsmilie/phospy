from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.science.evidence import PeptideEvidenceTable


def _base_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide_row_id": ["row_1", "row_2", "row_3"],
            "site_id": ["MAPK1;S10;", "MAPK1;S10;", "MAPK1;T12;"],
            "unique_feature_id": ["feat_1", "feat_2", "feat_3"],
            "gene_symbol": ["MAPK1", "MAPK1", "MAPK1"],
            "protein_accession": ["P28482", "P28482", "P28482"],
            "site_string": ["S10", "S10", "T12"],
            "sample_a": [10.0, 12.0, 9.0],
            "sample_b": [11.0, 13.0, 8.0],
            "peptide_sequence": ["AAAAA", "BBBBB", "CCCCC"],
            "modified_peptide_sequence": ["AA[+80]AAA", "BB[+80]BBB", "CC[+80]CCC"],
            "site_sequence": ["QAAAAAAA", "QBBBBBBB", "QCCCCCCC"],
            "localisation_confidence": [0.95, 0.93, 0.99],
            "missingness_flags": ["", "", ""],
            "imputation_flags": ["", "", ""],
            "multi_site": [False, False, False],
            "provenance_source": ["maxquant", "maxquant", "maxquant"],
        }
    )


def test_evidence_constructs_from_peptide_table_and_preserves_duplicate_site_rows() -> (
    None
):
    evidence = PeptideEvidenceTable(
        frame=_base_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    frame = evidence.to_dataframe()
    assert frame.shape[0] == 3
    duplicate_site_count = int((frame.loc[:, "site_id"] == "MAPK1;S10;").sum())
    assert duplicate_site_count == 2


def test_evidence_preserves_duplicate_peptide_sequences_as_distinct_rows() -> None:
    frame = _base_frame()
    frame.loc[1, "peptide_sequence"] = frame.loc[0, "peptide_sequence"]
    frame.loc[1, "modified_peptide_sequence"] = frame.loc[
        0, "modified_peptide_sequence"
    ]

    evidence = PeptideEvidenceTable(
        frame=frame,
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    records = evidence.records()
    duplicated_sequences = [
        record.peptide_row_id
        for record in records
        if record.peptide_sequence == "AAAAA"
    ]
    assert set(duplicated_sequences) == {"row_1", "row_2"}


def test_evidence_allows_missing_site_ids_and_mapping_excludes_missing_rows() -> None:
    frame = _base_frame()
    frame.loc[2, "site_id"] = None

    evidence = PeptideEvidenceTable(
        frame=frame,
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    mapping = evidence.site_mapping.to_dataframe()
    assert mapping.shape[0] == 2
    assert set(mapping.loc[:, "peptide_row_id"]) == {"row_1", "row_2"}


def test_evidence_keeps_multi_site_rows_unexploded_by_default_and_validates_mapping() -> (
    None
):
    frame = _base_frame().iloc[[0]].copy(deep=True)
    frame.loc[:, "multi_site"] = [True]
    frame.loc[:, "site_string"] = ["S10;T12"]
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["row_1", "row_1"],
            "site_id": ["MAPK1;S10;", "MAPK1;T12;"],
        }
    )

    evidence = PeptideEvidenceTable(
        frame=frame,
        sample_intensity_columns=("sample_a", "sample_b"),
        site_mapping=mapping,
    )
    assert evidence.to_dataframe().shape[0] == 1
    mapped = evidence.site_mapping.to_dataframe()
    assert mapped.shape[0] == 2
    assert set(mapped.loc[:, "site_id"]) == {"MAPK1;S10;", "MAPK1;T12;"}


def test_evidence_mapping_preserves_optional_weight_metadata_columns() -> None:
    frame = _base_frame().iloc[[0]].copy(deep=True)
    frame.loc[:, "multi_site"] = [True]
    frame.loc[:, "site_string"] = ["S10;T12"]
    mapping = pd.DataFrame(
        {
            "peptide_row_id": ["row_1", "row_1"],
            "site_id": ["MAPK1;S10;", "MAPK1;T12;"],
            "mapping_weight": [0.75, 0.25],
            "mapping_uncertainty": [True, True],
            "multi_site_policy": ["split_equal_weight", "split_equal_weight"],
            "is_multi_site": [True, True],
        }
    )

    evidence = PeptideEvidenceTable(
        frame=frame,
        sample_intensity_columns=("sample_a", "sample_b"),
        site_mapping=mapping,
    )
    mapped = evidence.site_mapping.to_dataframe()
    assert set(mapped.columns) == {
        "peptide_row_id",
        "site_id",
        "mapping_weight",
        "mapping_uncertainty",
        "multi_site_policy",
        "is_multi_site",
    }
    assert mapped.loc[:, "mapping_weight"].tolist() == [0.75, 0.25]


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (
            lambda frame: frame.drop(columns=["peptide_row_id"]),
            "peptide_evidence_table is missing required columns: peptide_row_id",
        ),
        (
            lambda frame: frame.assign(
                site_id=["not_a_site_id", "MAPK1;S10;", "MAPK1;T12;"]
            ),
            "site identifiers must use 'GENE;SITE;' format",
        ),
        (
            lambda frame: frame.assign(peptide_row_id=["row_1", "row_1", "row_3"]),
            "peptide_evidence_table.peptide_row_id must be unique",
        ),
    ],
)
def test_evidence_rejects_malformed_tables_with_clear_errors(
    mutator, message: str
) -> None:
    bad = mutator(_base_frame())
    with pytest.raises(PhosPyInputError, match=message):
        PeptideEvidenceTable(
            frame=bad,
            sample_intensity_columns=("sample_a", "sample_b"),
        )
