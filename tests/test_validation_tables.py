from __future__ import annotations

import pandas as pd
import pytest

from phospy.io import (
    DEFAULT_TEXT_ENCODING,
    default_text_encoding,
    infer_text_encoding,
    load_phospho_table,
    load_pred_mat,
    load_total_table,
)
from phospy.validation.errors import TableSchemaError
from phospy.validation.tables import (
    PhosphoInputSchema,
    PredictionScoreMatrixSchema,
    PredMatSchema,
    SiteMatrixSchema,
    SiteMatrixSourceSchema,
    TotalInputSchema,
)


def test_total_input_schema_rejects_non_numeric_group_values() -> None:
    frame = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": ["bad"],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )

    with pytest.raises(TableSchemaError, match="non-numeric values"):
        TotalInputSchema.validate(frame)


def test_total_input_schema_uses_canonical_gene_column_with_custom_value_columns() -> (
    None
):
    frame = pd.DataFrame(
        {
            "total_gene": ["PRKACA"],
            "sample_1": [1.0],
            "sample_2": [1.0],
        }
    )

    with pytest.raises(TableSchemaError, match="missing required columns: genes"):
        TotalInputSchema.validate(frame, total_cols=["sample_1", "sample_2"])


def test_phospho_input_schema_keeps_canonical_metadata_columns_with_custom_value_columns() -> (
    None
):
    frame = pd.DataFrame(
        {
            "uid": ["u1"],
            "phospho_gene": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "sample_1": [1.0],
            "sample_2": [1.0],
        }
    )

    with pytest.raises(TableSchemaError, match="missing required columns: gene_names"):
        PhosphoInputSchema.validate(frame, phospho_cols=["sample_1", "sample_2"])


def test_phospho_input_schema_rejects_malformed_gene_p_site() -> None:
    frame = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    with pytest.raises(TableSchemaError, match="malformed gene_p_site"):
        PhosphoInputSchema.validate(frame)


def test_load_phospho_table_rejects_duplicate_cleaned_columns(tmp_path) -> None:
    phospho_path = tmp_path / "phospho.tsv"
    phospho_path.write_text(
        "uid\tGene Names\tgene-names\tgene_p_site\tlocalization_prob\tcentralized_sequence"
        "\tp_group1\tp_group2\tp_group3\tp_group4\tp_group5\tp_group6\n"
        "u1\tPRKACA\tPRKACA\tPRKACA_S339\t0.95\tAAAAAA\t1\t1\t1\t1\t1\t1\n"
    )

    with pytest.raises(TableSchemaError, match="duplicate column names"):
        load_phospho_table(phospho_path)


def test_load_pred_mat_rejects_out_of_range_scores(tmp_path) -> None:
    pred_path = tmp_path / "pred.csv"
    pd.DataFrame(
        {
            "PRKACA": [1.2],
            "BTK": [0.8],
        },
        index=["PRKACA;S339;"],
    ).to_csv(pred_path)

    with pytest.raises(TableSchemaError, match="outside the allowed range"):
        load_pred_mat(pred_path)


def test_pred_mat_schema_rejects_zero_column_frames() -> None:
    frame = pd.DataFrame(index=["SITE_1", "SITE_2"], dtype=float)

    with pytest.raises(
        TableSchemaError, match="pred_mat must contain at least one kinase column"
    ):
        PredMatSchema.validate(frame)


def test_site_matrix_schema_rejects_duplicate_index() -> None:
    frame = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
        },
        index=["SITE_1", "SITE_1"],
    )

    with pytest.raises(TableSchemaError, match="duplicate index entries"):
        SiteMatrixSchema.validate(frame)


def test_load_total_table_returns_numeric_frame(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    total_path.write_text(
        "genes\tgroup1\tgroup2\tgroup3\tgroup4\tgroup5\tgroup6\n"
        "PRKACA\t1\t2\t3\t4\t5\t6\n"
    )

    loaded = load_total_table(total_path)

    assert loaded["group1"].dtype.kind in {"f", "i"}
    assert loaded.loc[0, "genes"] == "PRKACA"


def test_load_phospho_table_uses_explicit_encoding_when_provided(tmp_path) -> None:
    phospho_path = tmp_path / "phospho-utf16.tsv"
    phospho_path.write_text(
        "uid\tgene_names\tgene_p_site\tlocalization_prob\tcentralized_sequence"
        "\tp_group1\tp_group2\tp_group3\tp_group4\tp_group5\tp_group6\n"
        "u1\tPRKACA\tPRKACA_S339\t0.95\tAAAAAA\t1\t1\t1\t1\t1\t1\n",
        encoding="utf-16",
    )

    loaded = load_phospho_table(phospho_path, encoding="utf-16")

    assert loaded.loc[0, "gene_names"] == "PRKACA"
    assert default_text_encoding(phospho_path) == DEFAULT_TEXT_ENCODING
    assert infer_text_encoding(phospho_path) == DEFAULT_TEXT_ENCODING


def test_site_matrix_source_schema_rejects_empty_or_extra_delimiter_gene_p_site() -> (
    None
):
    frame = pd.DataFrame(
        {
            "gene_p_site": ["GENE_", "_S123", "GENE__S1"],
            "centralized_sequence": ["AAAAAA", "BBBBBB", "CCCCCC"],
            "sample_1": [1.0, 2.0, 3.0],
        }
    )

    with pytest.raises(
        TableSchemaError,
        match=(
            "site-matrix source table contains malformed gene_p_site values "
            "that cannot be split into non-empty gene and site parts using a single underscore"
        ),
    ):
        SiteMatrixSourceSchema.validate(
            frame,
            gene_p_site_col="gene_p_site",
            sequence_col="centralized_sequence",
            value_cols=["sample_1"],
        )


def test_prediction_score_matrix_schema_rejects_non_finite_scores() -> None:
    frame = pd.DataFrame(
        {
            "KINASE_A": [0.8, float("nan")],
            "KINASE_B": [0.2, 0.4],
        },
        index=["SITE_1", "SITE_2"],
    )

    with pytest.raises(TableSchemaError, match="non-finite values"):
        PredictionScoreMatrixSchema.validate(frame)


def test_site_matrix_schema_rejects_non_finite_values() -> None:
    frame = pd.DataFrame(
        {
            "sample_1": [1.0, float("-inf")],
            "sample_2": [1.1, 2.1],
        },
        index=["SITE_1", "SITE_2"],
    )

    with pytest.raises(TableSchemaError, match="non-finite values"):
        SiteMatrixSchema.validate(frame)
