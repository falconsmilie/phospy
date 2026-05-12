from __future__ import annotations

import pandas as pd
import pytest

from phospy.sites.identity import (
    validate_identity_optional_columns,
    validate_no_conflicting_identity_collisions,
)


def test_validate_identity_optional_columns_accepts_missing_and_strings() -> None:
    frame = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["A" * 31, "A" * 31],
            "protein_id": ["P28482", pd.NA],
            "protein_accession": ["P28482-1", ""],
            "source_namespace": ["uniprot", " "],
        },
        index=pd.Index(["MAPK14;Y182;", "MAPK14;Y182;"], name="site_id"),
    )

    validate_identity_optional_columns(
        site_metadata=frame,
        field_name="dataset.site_metadata",
        error_type=ValueError,
    )


def test_validate_identity_optional_columns_rejects_non_string_values() -> None:
    frame = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["A" * 31],
            "protein_accession": [123],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    with pytest.raises(ValueError, match="optional identity columns"):
        validate_identity_optional_columns(
            site_metadata=frame,
            field_name="dataset.site_metadata",
            error_type=ValueError,
        )


def test_identity_collision_rejects_same_display_id_with_conflicting_protein_context() -> (
    None
):
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["A" * 31, "A" * 31],
            "protein_id": ["P28482", "P28482"],
            "protein_accession": ["P28482-1", "P28482-2"],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    constructed_site_ids = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=site_metadata.index.copy(),
        name="site_id",
    )

    with pytest.raises(
        ValueError,
        match="conflicting scientific identities for duplicate display site IDs",
    ):
        validate_no_conflicting_identity_collisions(
            site_metadata=site_metadata,
            display_ids=constructed_site_ids,
            field_name="dataset preprocessing.site_matrix",
            error_type=ValueError,
        )


def test_identity_collision_allows_semantically_identical_duplicates() -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["A" * 31, "A" * 31],
            "protein_id": ["P28482-1", "P28482-1"],
            "protein_accession": ["P28482-1", "P28482-1"],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    constructed_site_ids = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=site_metadata.index.copy(),
        name="site_id",
    )

    validate_no_conflicting_identity_collisions(
        site_metadata=site_metadata,
        display_ids=constructed_site_ids,
        field_name="dataset preprocessing.site_matrix",
        error_type=ValueError,
    )
