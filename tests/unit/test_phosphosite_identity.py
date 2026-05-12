from __future__ import annotations

import pytest

from phospy.sites.identity import PhosphositeIdentity, build_phosphosite_identity


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
    with pytest.raises(ValueError, match="must use '<residue><position>'"):
        build_phosphosite_identity(
            display_id="MAPK14;Y182-T180;",
            gene_symbol="MAPK14",
            site="Y182-T180",
            field_name="test.site_metadata.row",
            error_type=ValueError,
        )
