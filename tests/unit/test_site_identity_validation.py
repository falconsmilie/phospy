from __future__ import annotations

import pandas as pd
import pytest

from phospy.science.sites.identity import (
    validate_identity_optional_columns,
    validate_no_conflicting_identity_collisions,
)
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.validation.datasets.protein_scoped_site_identity import (
    enforce_display_id_column,
    enforce_site_key_column,
    enforce_site_key_index,
    enforce_site_key_matches_metadata,
    enforce_unique_site_key_identity,
)


def _display_id_collision_inputs(
    *,
    protein_ids: list[object] | None = None,
    protein_accessions: list[object] | None = None,
    isoform_ids: list[object] | None = None,
) -> tuple[pd.DataFrame, pd.Series]:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "Y182"]
            ],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    if protein_ids is not None:
        site_metadata.loc[:, "protein_id"] = pd.Series(
            protein_ids,
            index=site_metadata.index,
            dtype="object",
        )
    if protein_accessions is not None:
        site_metadata.loc[:, "protein_accession"] = pd.Series(
            protein_accessions,
            index=site_metadata.index,
            dtype="object",
        )
    if isoform_ids is not None:
        site_metadata.loc[:, "isoform_id"] = pd.Series(
            isoform_ids,
            index=site_metadata.index,
            dtype="object",
        )
    constructed_site_ids = pd.Series(
        ["MAPK14;Y182;", "MAPK14;Y182;"],
        index=site_metadata.index.copy(),
        name="site_id",
    )
    return site_metadata, constructed_site_ids


def test_validate_identity_optional_columns_accepts_missing_and_strings() -> None:
    frame = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "Y182"]
            ],
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
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182"]
            ],
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
    site_metadata, constructed_site_ids = _display_id_collision_inputs(
        protein_ids=["P28482", "P28482"],
        protein_accessions=["P28482-1", "P28482-2"],
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
    site_metadata, constructed_site_ids = _display_id_collision_inputs(
        protein_ids=["P28482-1", "P28482-1"],
        protein_accessions=["P28482-1", "P28482-1"],
    )

    validate_no_conflicting_identity_collisions(
        site_metadata=site_metadata,
        display_ids=constructed_site_ids,
        field_name="dataset preprocessing.site_matrix",
        error_type=ValueError,
    )


def test_identity_collision_allows_same_display_id_with_matching_protein_id() -> None:
    site_metadata, constructed_site_ids = _display_id_collision_inputs(
        protein_ids=["P28482", "P28482"],
    )

    validate_no_conflicting_identity_collisions(
        site_metadata=site_metadata,
        display_ids=constructed_site_ids,
        field_name="dataset preprocessing.site_matrix",
        error_type=ValueError,
    )


def test_identity_collision_rejects_same_display_id_with_conflicting_protein_id() -> (
    None
):
    site_metadata, constructed_site_ids = _display_id_collision_inputs(
        protein_ids=["P28482", "Q12345"],
    )

    with pytest.raises(
        ValueError,
        match="protein_id=\\['P28482', 'Q12345'\\]",
    ):
        validate_no_conflicting_identity_collisions(
            site_metadata=site_metadata,
            display_ids=constructed_site_ids,
            field_name="dataset preprocessing.site_matrix",
            error_type=ValueError,
        )


def test_identity_collision_protein_id_missing_context_follows_strict_policy() -> None:
    site_metadata, constructed_site_ids = _display_id_collision_inputs(
        protein_ids=["P28482", pd.NA],
    )

    with pytest.raises(
        ValueError,
        match="protein_id=\\['P28482', None\\]",
    ):
        validate_no_conflicting_identity_collisions(
            site_metadata=site_metadata,
            display_ids=constructed_site_ids,
            field_name="dataset preprocessing.site_matrix",
            error_type=ValueError,
        )


def test_identity_collision_rejects_same_display_id_with_conflicting_accession() -> (
    None
):
    site_metadata, constructed_site_ids = _display_id_collision_inputs(
        protein_accessions=["P28482-1", "P28482-2"],
    )

    with pytest.raises(
        ValueError,
        match="protein_accession=\\['P28482-1', 'P28482-2'\\]",
    ):
        validate_no_conflicting_identity_collisions(
            site_metadata=site_metadata,
            display_ids=constructed_site_ids,
            field_name="dataset preprocessing.site_matrix",
            error_type=ValueError,
        )


def test_identity_collision_rejects_same_display_id_with_conflicting_isoform() -> None:
    site_metadata, constructed_site_ids = _display_id_collision_inputs(
        isoform_ids=["isoform-1", "isoform-2"],
    )

    with pytest.raises(
        ValueError,
        match="isoform_id=\\['isoform-1', 'isoform-2'\\]",
    ):
        validate_no_conflicting_identity_collisions(
            site_metadata=site_metadata,
            display_ids=constructed_site_ids,
            field_name="dataset preprocessing.site_matrix",
            error_type=ValueError,
        )


def _protein_scoped_metadata_rows() -> pd.DataFrame:
    rows = pd.DataFrame(
        {
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
            "organism": ["rat", "rat"],
            "protein_namespace": ["protein_accession", "protein_accession"],
            "protein_identifier": ["P28482-1", "Q9WVS8-1"],
            "site": ["Y182", "Y182"],
            "residue": ["Y", "Y"],
            "position": [182, 182],
            "isoform_id": [pd.NA, pd.NA],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    rows.loc[:, "site_key"] = [
        encode_site_key(
            build_protein_scoped_site_key(
                organism=str(rows.at[row_id, "organism"]),
                protein_namespace=str(rows.at[row_id, "protein_namespace"]),
                protein_identifier=str(rows.at[row_id, "protein_identifier"]),
                residue=str(rows.at[row_id, "residue"]),
                position=int(rows.at[row_id, "position"]),
                isoform_id=None,
                field_name=f"test[{row_id!r}].site_key",
                error_type=ValueError,
            )
        )
        for row_id in rows.index.tolist()
    ]
    return rows


def test_protein_scoped_identity_rejects_duplicate_site_key() -> None:
    site_metadata = _protein_scoped_metadata_rows()
    site_metadata.loc["row_b", "site_key"] = site_metadata.loc["row_a", "site_key"]

    with pytest.raises(ValueError, match="site_key must be unique"):
        enforce_unique_site_key_identity(
            site_metadata=site_metadata,
            field_name="dataset.site_metadata",
            error_type=ValueError,
        )


def test_protein_scoped_identity_allows_duplicate_display_id_when_site_key_differs() -> (
    None
):
    site_metadata = _protein_scoped_metadata_rows()

    enforce_display_id_column(
        site_metadata=site_metadata,
        field_name="dataset.site_metadata",
        error_type=ValueError,
    )
    enforce_unique_site_key_identity(
        site_metadata=site_metadata,
        field_name="dataset.site_metadata",
        error_type=ValueError,
    )


def test_protein_scoped_identity_rejects_missing_site_key_column() -> None:
    site_metadata = _protein_scoped_metadata_rows().drop(columns=["site_key"])

    with pytest.raises(ValueError, match="missing required columns: site_key"):
        enforce_site_key_column(
            site_metadata=site_metadata,
            field_name="dataset.site_metadata",
            error_type=ValueError,
        )


def test_protein_scoped_identity_rejects_missing_display_id_column() -> None:
    site_metadata = _protein_scoped_metadata_rows().drop(columns=["display_id"])

    with pytest.raises(ValueError, match="missing required columns: display_id"):
        enforce_display_id_column(
            site_metadata=site_metadata,
            field_name="dataset.site_metadata",
            error_type=ValueError,
        )


def test_protein_scoped_identity_rejects_site_key_metadata_mismatch() -> None:
    site_metadata = _protein_scoped_metadata_rows()
    site_metadata.loc["row_b", "site_key"] = encode_site_key(
        build_protein_scoped_site_key(
            organism=str(site_metadata.at["row_b", "organism"]),
            protein_namespace=str(site_metadata.at["row_b", "protein_namespace"]),
            protein_identifier=str(site_metadata.at["row_b", "protein_identifier"]),
            residue=str(site_metadata.at["row_b", "residue"]),
            position=308,
            isoform_id=None,
            field_name="test['row_b'].site_key",
            error_type=ValueError,
        )
    )

    with pytest.raises(ValueError, match="must match metadata-derived"):
        enforce_site_key_matches_metadata(
            site_metadata=site_metadata,
            field_name="dataset.site_metadata",
            error_type=ValueError,
        )


def test_protein_scoped_identity_rejects_mixed_missing_and_specified_isoform_scope() -> (
    None
):
    site_metadata = _protein_scoped_metadata_rows()
    site_metadata.loc[:, "protein_identifier"] = ["P28482-1", "P28482-1"]
    site_metadata.loc[:, "isoform_id"] = [pd.NA, "2"]
    site_metadata.loc[:, "site_key"] = [
        encode_site_key(
            build_protein_scoped_site_key(
                organism=str(site_metadata.at[row_id, "organism"]),
                protein_namespace=str(site_metadata.at[row_id, "protein_namespace"]),
                protein_identifier=str(site_metadata.at[row_id, "protein_identifier"]),
                residue=str(site_metadata.at[row_id, "residue"]),
                position=int(site_metadata.at[row_id, "position"]),
                isoform_id=(
                    None
                    if pd.isna(site_metadata.at[row_id, "isoform_id"])
                    else str(site_metadata.at[row_id, "isoform_id"])
                ),
                field_name=f"test[{row_id!r}].site_key",
                error_type=ValueError,
            )
        )
        for row_id in site_metadata.index.tolist()
    ]

    with pytest.raises(ValueError, match="mixed isoform scope"):
        enforce_unique_site_key_identity(
            site_metadata=site_metadata,
            field_name="dataset.site_metadata",
            error_type=ValueError,
        )


def test_protein_scoped_identity_site_key_index_enforcement() -> None:
    site_metadata = _protein_scoped_metadata_rows().copy(deep=True)
    site_metadata.index = pd.Index(
        site_metadata.loc[:, "site_key"].tolist(), name="row"
    )

    enforce_site_key_index(
        site_metadata=site_metadata,
        field_name="dataset.site_metadata",
        error_type=ValueError,
    )

    invalid = site_metadata.copy(deep=True)
    invalid.index = pd.Index(["row_a", "row_b"], name="row")
    with pytest.raises(ValueError, match="index must match"):
        enforce_site_key_index(
            site_metadata=invalid,
            field_name="dataset.site_metadata",
            error_type=ValueError,
        )
