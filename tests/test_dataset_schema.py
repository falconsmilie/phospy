from __future__ import annotations

import pandas as pd
import pytest

from phospy import DatasetSchema, PhosphoDataset
from phospy.validation.errors import InputCompatibilityError


def test_dataset_schema_validates_aligned_column_groups() -> None:
    with pytest.raises(
        InputCompatibilityError,
        match="same number of total and phospho",
    ):
        DatasetSchema(
            total_cols=("group1", "group2"),
            phospho_cols=("p_group1",),
        )


def test_dataset_stores_schema_object() -> None:
    schema = DatasetSchema(
        corrected_cols=(
            "sample_a",
            "sample_b",
            "sample_c",
            "sample_d",
            "sample_e",
            "sample_f",
        ),
    )
    dataset = PhosphoDataset(
        total_df=pd.DataFrame(
            {
                "genes": ["PRKACA"],
                "group1": [1.0],
                "group2": [1.0],
                "group3": [1.0],
                "group4": [1.0],
                "group5": [1.0],
                "group6": [1.0],
            }
        ),
        phospho_df=pd.DataFrame(
            {
                "uid": ["u1"],
                "gene_names": ["PRKACA"],
                "gene_p_site": ["PRKACA_S339"],
                "localization_prob": [0.95],
                "centralized_sequence": ["AAAAAA"],
                "p_group1": [8.0],
                "p_group2": [8.0],
                "p_group3": [8.0],
                "p_group4": [8.0],
                "p_group5": [8.0],
                "p_group6": [8.0],
            }
        ),
        schema=schema,
    )

    assert dataset.schema is schema
