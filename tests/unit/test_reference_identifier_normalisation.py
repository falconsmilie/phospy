from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.validation import ReferenceValidationError
from phospy.references.identifiers import (
    normalise_reference_kinase_id,
    normalise_reference_site_id,
)
from phospy.references.models import Organism, ReferenceBundle


def test_kinase_identifier_normaliser_normalises_mixed_case_values() -> None:
    records = [
        normalise_reference_kinase_id(
            "akt1",
            table_name="references.kinase_substrate_map",
            column_name="kinase",
            row_position=0,
        ),
        normalise_reference_kinase_id(
            "Akt1",
            table_name="references.kinase_substrate_map",
            column_name="kinase",
            row_position=1,
        ),
        normalise_reference_kinase_id(
            "AKT1",
            table_name="references.kinase_substrate_map",
            column_name="kinase",
            row_position=2,
        ),
    ]

    assert [record.normalised_value for record in records] == ["AKT1", "AKT1", "AKT1"]


def test_site_identifier_normaliser_normalises_supported_inputs() -> None:
    records = [
        normalise_reference_site_id(
            "mapk1;s123",
            table_name="references.kinase_substrate_map",
            column_name="substrate_site",
            row_position=0,
        ),
        normalise_reference_site_id(
            "Mapk1 ; s123",
            table_name="references.kinase_substrate_map",
            column_name="substrate_site",
            row_position=1,
        ),
        normalise_reference_site_id(
            "MAPK1;S123;",
            table_name="references.kinase_substrate_map",
            column_name="substrate_site",
            row_position=2,
        ),
    ]

    assert [record.normalised_value for record in records] == [
        "MAPK1;S123;",
        "MAPK1;S123;",
        "MAPK1;S123;",
    ]


def test_reference_bundle_rejects_duplicate_pairs_after_kinase_normalisation() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["akt1", "AKT1"],
                    "substrate_site": ["MAPK1;S123;", "MAPK1;S123;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK1;S123;"], name="site_id"),
            ),
        )
