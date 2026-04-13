from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import TableSchemaError
from phospy.matrices import SiteMatrixPolicy, build_site_matrix


def test_build_site_matrix_creates_site_ids_and_deduplicates_by_mean() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "PRKACA_S339", "BTK_Y551"],
            "centralized_sequence": ["AAAAAA", "BBBBBB", "CCCCCC"],
            "phospho_corrected_1": [1.0, 10.0, 3.0],
            "phospho_corrected_2": [1.0, 10.0, 3.0],
            "phospho_corrected_3": [1.0, 10.0, 3.0],
            "phospho_corrected_4": [1.0, 10.0, 3.0],
            "phospho_corrected_5": [1.0, 10.0, 3.0],
            "phospho_corrected_6": [1.0, 10.0, 3.0],
        }
    )
    phosr_input, matrix, sequences = build_site_matrix(
        df=df,
        gene_p_site_col="gene_p_site",
        sequence_col="centralized_sequence",
        value_cols=[
            "phospho_corrected_1",
            "phospho_corrected_2",
            "phospho_corrected_3",
            "phospho_corrected_4",
            "phospho_corrected_5",
            "phospho_corrected_6",
        ],
    )
    assert "PRKACA;S339;" in matrix.index
    assert matrix.loc["PRKACA;S339;", "phospho_corrected_1"] == 10.0
    assert sequences.loc["PRKACA;S339;"] == "BBBBBB"
    assert phosr_input.shape[0] == 2


def test_build_site_matrix_exposes_row_drop_stats() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "BTK_Y551", "BTK_Y551"],
            "centralized_sequence": ["AAAAAA", None, "CCCCCC"],
            "phospho_corrected_1": [1.0, 2.0, 3.0],
            "phospho_corrected_2": [1.0, 2.0, None],
        }
    )
    phosr_input, matrix, sequences = build_site_matrix(
        df=df,
        gene_p_site_col="gene_p_site",
        sequence_col="centralized_sequence",
        value_cols=["phospho_corrected_1", "phospho_corrected_2"],
    )

    stats = phosr_input.attrs["row_drop_stats"]
    assert stats["input_rows"] == 3
    assert stats["dropped_missing_sequence"] == 1
    assert stats["dropped_incomplete_values"] == 1
    assert stats["retained_rows"] == 1
    assert matrix.attrs["row_drop_stats"] == stats
    assert sequences.attrs["row_drop_stats"] == stats


def test_build_site_matrix_raises_table_schema_error_for_malformed_gene_p_site() -> (
    None
):
    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA", "BTK_Y551"],
            "centralized_sequence": ["AAAAAA", "CCCCCC"],
            "phospho_corrected_1": [1.0, 3.0],
            "phospho_corrected_2": [1.0, 3.0],
        }
    )

    with pytest.raises(
        TableSchemaError,
        match=(
            "site-matrix source table contains malformed gene_p_site values "
            "that cannot be split into non-empty gene and site parts using a single underscore"
        ),
    ):
        build_site_matrix(
            df=df,
            gene_p_site_col="gene_p_site",
            sequence_col="centralized_sequence",
            value_cols=["phospho_corrected_1", "phospho_corrected_2"],
        )


def test_site_matrix_builder_reports_row_drop_diagnostics_when_all_rows_are_dropped() -> (
    None
):
    from phospy.preprocessing import SiteMatrixBuilder

    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "BTK_Y551", "SRC_Y419"],
            "centralized_sequence": [None, "CCCCCC", "DDDDDD"],
            "phospho_corrected_1": [1.0, None, 3.0],
            "phospho_corrected_2": [1.0, 2.0, None],
        }
    )

    builder = SiteMatrixBuilder(
        value_cols=["phospho_corrected_1", "phospho_corrected_2"]
    )

    with pytest.raises(TableSchemaError, match="row-drop diagnostics") as exc_info:
        builder.build(df)

    message = str(exc_info.value)
    assert "site matrix must contain at least one phosphosite row" in message
    assert "dropped_missing_sequence=1" in message
    assert "dropped_incomplete_values=2" in message
    assert "deduplicated_site_rows=0" in message
    assert "other_dropped_rows=0" in message
    assert "retained_rows=0" in message


def test_build_site_matrix_rejects_empty_or_extra_delimiter_gene_p_site_values() -> (
    None
):
    invalid_values = ["GENE_", "_S123", "GENE__S1"]

    for invalid_value in invalid_values:
        df = pd.DataFrame(
            {
                "gene_p_site": [invalid_value],
                "centralized_sequence": ["AAAAAA"],
                "phospho_corrected_1": [1.0],
                "phospho_corrected_2": [1.0],
            }
        )

        with pytest.raises(
            TableSchemaError,
            match=(
                "site-matrix source table contains malformed gene_p_site values "
                "that cannot be split into non-empty gene and site parts using a single underscore"
            ),
        ):
            build_site_matrix(
                df=df,
                gene_p_site_col="gene_p_site",
                sequence_col="centralized_sequence",
                value_cols=["phospho_corrected_1", "phospho_corrected_2"],
            )


def test_build_site_matrix_can_keep_first_duplicate_row() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "PRKACA_S339"],
            "centralized_sequence": ["AAAAAA", "BBBBBB"],
            "phospho_corrected_1": [1.0, 10.0],
            "phospho_corrected_2": [1.0, 10.0],
        }
    )

    phosr_input, matrix, _ = build_site_matrix(
        df=df,
        gene_p_site_col="gene_p_site",
        sequence_col="centralized_sequence",
        value_cols=["phospho_corrected_1", "phospho_corrected_2"],
        policy=SiteMatrixPolicy(duplicate_site_strategy="first"),
    )

    assert matrix.loc["PRKACA;S339;", "phospho_corrected_1"] == 1.0
    assert phosr_input.attrs["row_drop_stats"]["duplicate_site_strategy"] == "first"


def test_build_site_matrix_can_aggregate_duplicate_rows_by_mean() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "PRKACA_S339"],
            "centralized_sequence": ["AAAAAA", "BBBBBB"],
            "phospho_corrected_1": [1.0, 3.0],
            "phospho_corrected_2": [2.0, 4.0],
        }
    )

    _, matrix, _ = build_site_matrix(
        df=df,
        gene_p_site_col="gene_p_site",
        sequence_col="centralized_sequence",
        value_cols=["phospho_corrected_1", "phospho_corrected_2"],
        policy=SiteMatrixPolicy(duplicate_site_strategy="aggregate_mean"),
    )

    assert float(matrix.loc["PRKACA;S339;", "phospho_corrected_1"]) == pytest.approx(
        2.0
    )
    assert float(matrix.loc["PRKACA;S339;", "phospho_corrected_2"]) == pytest.approx(
        3.0
    )


def test_build_site_matrix_can_reject_duplicate_rows_explicitly() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "PRKACA_S339"],
            "centralized_sequence": ["AAAAAA", "BBBBBB"],
            "phospho_corrected_1": [1.0, 3.0],
            "phospho_corrected_2": [2.0, 4.0],
        }
    )

    with pytest.raises(TableSchemaError, match="duplicate_site_strategy='error'"):
        build_site_matrix(
            df=df,
            gene_p_site_col="gene_p_site",
            sequence_col="centralized_sequence",
            value_cols=["phospho_corrected_1", "phospho_corrected_2"],
            policy=SiteMatrixPolicy(duplicate_site_strategy="error"),
        )


def test_site_matrix_builder_build_owned_routes_through_no_copy_validation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from phospy.preprocessing import SiteMatrixBuilder
    from phospy.validation.schema.tables import SiteMatrixSourceSchema

    source_df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "BTK_Y551"],
            "centralized_sequence": ["AAAAAA", "BBBBBB"],
            "phospho_corrected_1": [1.0, 2.0],
            "phospho_corrected_2": [1.0, 2.0],
        }
    )
    builder = SiteMatrixBuilder(
        value_cols=["phospho_corrected_1", "phospho_corrected_2"]
    )

    seen_copy_frame: list[bool] = []
    original_validate = SiteMatrixSourceSchema.validate

    def counting_validate(*args, **kwargs):
        seen_copy_frame.append(bool(kwargs.get("copy_frame", True)))
        return original_validate(*args, **kwargs)

    monkeypatch.setattr(SiteMatrixSourceSchema, "validate", counting_validate)

    public_result = builder.build(source_df)
    owned_result = builder.build_owned(source_df.copy(deep=True))

    assert not public_result.matrix.empty
    assert not owned_result.matrix.empty
    assert seen_copy_frame == [True, False]
