from __future__ import annotations

import pytest

from phospy.science.sites.identity import (
    PhosphositeIdentity,
    build_phosphosite_identity,
)
from phospy.science.sites.site_keys import (
    ProteinScopedPhosphositeKey,
    build_protein_scoped_site_key,
    decode_site_key,
    encode_site_key,
)


def test_build_phosphosite_identity_parses_canonical_display_id() -> None:
    identity = build_phosphosite_identity(
        display_id="MAPK14;Y182;",
        gene_symbol="MAPK14",
        site="Y182",
        protein_id="P28482",
        protein_accession="P28482-1",
        source_namespace="uniprot",
        source_site_id="P28482-1:Y182",
        field_name="test.site_metadata.row",
        error_type=ValueError,
    )

    assert isinstance(identity, PhosphositeIdentity)
    assert identity.display_id == "MAPK14;Y182;"
    assert identity.gene_symbol == "MAPK14"
    assert identity.residue == "Y"
    assert identity.position == 182
    assert identity.protein_id == "P28482"
    assert identity.protein_accession == "P28482-1"
    assert identity.source_namespace == "uniprot"
    assert identity.source_site_id == "P28482-1:Y182"
    assert identity.has_protein_context()


def test_build_phosphosite_identity_rejects_mismatched_display_and_components() -> None:
    with pytest.raises(ValueError, match="inconsistent display and metadata identity"):
        build_phosphosite_identity(
            display_id="MAPK14;Y182;",
            gene_symbol="MAPK14",
            site="T180",
            field_name="test.site_metadata.row",
            error_type=ValueError,
        )


def test_build_phosphosite_identity_rejects_non_parseable_site_token() -> None:
    with pytest.raises(ValueError, match="must use strict 'S/T/Y<position>'"):
        build_phosphosite_identity(
            display_id="MAPK14;Y182-T180;",
            gene_symbol="MAPK14",
            site="Y182-T180",
            field_name="test.site_metadata.row",
            error_type=ValueError,
        )


@pytest.mark.parametrize("site_token", ["S1", "T45", "Y999"])
def test_build_phosphosite_identity_accepts_strict_sty_site_tokens(
    site_token: str,
) -> None:
    identity = build_phosphosite_identity(
        display_id=f"MAPK14;{site_token};",
        gene_symbol="MAPK14",
        site=site_token,
        field_name="test.site_metadata.row",
        error_type=ValueError,
    )
    assert identity.display_id == f"MAPK14;{site_token};"
    assert identity.residue == site_token[0]
    assert identity.position == int(site_token[1:])


@pytest.mark.parametrize("site_token", ["FOO", "A123", "S0", "S", "123"])
def test_build_phosphosite_identity_rejects_invalid_site_tokens_by_default(
    site_token: str,
) -> None:
    with pytest.raises(ValueError, match="must use strict 'S/T/Y<position>'"):
        build_phosphosite_identity(
            display_id=f"MAPK14;{site_token};",
            gene_symbol="MAPK14",
            site=site_token,
            field_name="test.site_metadata.row",
            error_type=ValueError,
        )


def test_build_phosphosite_identity_allows_opaque_site_token_with_explicit_opt_in() -> (
    None
):
    identity = build_phosphosite_identity(
        display_id="MAPK14;FOO;",
        gene_symbol="MAPK14",
        site="FOO",
        allow_opaque_site_values=True,
        field_name="test.site_metadata.row",
        error_type=ValueError,
    )
    assert identity.display_id == "MAPK14;FOO;"
    assert identity.residue is None
    assert identity.position is None


def test_build_encode_decode_protein_scoped_site_key_round_trips() -> None:
    key = build_protein_scoped_site_key(
        organism="rat",
        protein_namespace="protein_accession",
        protein_identifier="P31750",
        residue="s",
        position=473,
        isoform_id="1",
        field_name="test.site_key",
        error_type=ValueError,
    )

    encoded = encode_site_key(key)
    decoded = decode_site_key(
        encoded,
        field_name="test.site_key",
        error_type=ValueError,
    )

    assert isinstance(key, ProteinScopedPhosphositeKey)
    assert encoded == (
        "phospy:v1|organism=rat|protein_namespace=protein_accession|"
        "protein_identifier=P31750|residue=S|position=473|isoform_id=1"
    )
    assert decoded == key
    assert not hasattr(decoded, "display_id")


def test_encode_site_key_is_deterministic() -> None:
    key = build_protein_scoped_site_key(
        organism="human",
        protein_namespace="uniprot",
        protein_identifier="P04637",
        residue="Y",
        position=220,
        field_name="test.site_key",
        error_type=ValueError,
    )
    assert encode_site_key(key) == encode_site_key(key)


def test_site_key_encoding_escapes_delimiters_and_decodes_losslessly() -> None:
    key = build_protein_scoped_site_key(
        organism="rat|cohort",
        protein_namespace="protein=accession",
        protein_identifier="P31750/alpha",
        residue="T",
        position=308,
        isoform_id="iso|1=alpha",
        field_name="test.site_key",
        error_type=ValueError,
    )
    encoded = encode_site_key(key)

    assert "%7C" in encoded
    assert "%3D" in encoded
    assert (
        decode_site_key(
            encoded,
            field_name="test.site_key",
            error_type=ValueError,
        )
        == key
    )


@pytest.mark.parametrize("invalid_residue", ["A", "X", "ST", ""])
def test_build_protein_scoped_site_key_rejects_invalid_residues(
    invalid_residue: str,
) -> None:
    with pytest.raises(ValueError, match="residue"):
        build_protein_scoped_site_key(
            organism="human",
            protein_namespace="uniprot",
            protein_identifier="P31750",
            residue=invalid_residue,
            position=473,
            field_name="test.site_key",
            error_type=ValueError,
        )


@pytest.mark.parametrize("invalid_position", [0, -1, None, "473", 4.2, True])
def test_build_protein_scoped_site_key_rejects_invalid_positions(
    invalid_position: object,
) -> None:
    with pytest.raises(ValueError, match="position"):
        build_protein_scoped_site_key(
            organism="human",
            protein_namespace="uniprot",
            protein_identifier="P31750",
            residue="S",
            position=invalid_position,
            field_name="test.site_key",
            error_type=ValueError,
        )


def test_build_protein_scoped_site_key_rejects_missing_organism() -> None:
    with pytest.raises(ValueError, match="organism"):
        build_protein_scoped_site_key(
            organism=" ",
            protein_namespace="uniprot",
            protein_identifier="P31750",
            residue="S",
            position=473,
            field_name="test.site_key",
            error_type=ValueError,
        )


def test_build_protein_scoped_site_key_rejects_missing_protein_namespace() -> None:
    with pytest.raises(ValueError, match="protein_namespace"):
        build_protein_scoped_site_key(
            organism="human",
            protein_namespace="",
            protein_identifier="P31750",
            residue="S",
            position=473,
            field_name="test.site_key",
            error_type=ValueError,
        )


def test_build_protein_scoped_site_key_rejects_missing_protein_identifier() -> None:
    with pytest.raises(ValueError, match="protein_identifier"):
        build_protein_scoped_site_key(
            organism="human",
            protein_namespace="uniprot",
            protein_identifier="   ",
            residue="S",
            position=473,
            field_name="test.site_key",
            error_type=ValueError,
        )


def test_decode_site_key_rejects_missing_required_position_field() -> None:
    encoded = (
        "phospy:v1|organism=human|protein_namespace=uniprot|"
        "protein_identifier=P31750|residue=S"
    )
    with pytest.raises(ValueError, match="missing required encoded fields"):
        decode_site_key(
            encoded,
            field_name="test.site_key",
            error_type=ValueError,
        )


def test_decode_site_key_rejects_non_integer_position_field() -> None:
    encoded = (
        "phospy:v1|organism=human|protein_namespace=uniprot|"
        "protein_identifier=P31750|residue=S|position=not-an-integer"
    )
    with pytest.raises(ValueError, match="position must be an integer"):
        decode_site_key(
            encoded,
            field_name="test.site_key",
            error_type=ValueError,
        )
