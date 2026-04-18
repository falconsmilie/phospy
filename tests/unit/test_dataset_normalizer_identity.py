from __future__ import annotations

import pandas as pd
import pytest

from phospy.datasets.builders.normalizer import DatasetConventionNormalizer
from phospy.errors.input import UnsupportedInputFormatError


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def test_normalizer_supports_documented_site_metadata_aliases() -> None:
    normalized = DatasetConventionNormalizer().run(
        phospho=_phospho(),
        site_metadata=pd.DataFrame(
            {
                "gene_name": ["MAPK14"],
                "protein_id": ["P28482-2"],
                "residue": ["Y182"],
                "centralized_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        sample_metadata=None,
        total=None,
    )
    site_metadata = normalized.site_metadata
    assert site_metadata.loc["MAPK14;Y182;", "gene_symbol"] == "MAPK14"
    assert site_metadata.loc["MAPK14;Y182;", "site"] == "Y182"
    assert (
        site_metadata.loc["MAPK14;Y182;", "site_sequence"]
        == "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"
    )
    assert site_metadata.loc["MAPK14;Y182;", "protein_id"] == "P28482-2"


def test_normalizer_rejects_unsupported_ambiguous_sequence_column() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="column 'sequence' is ambiguous and unsupported",
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )


def test_normalizer_rejects_unsupported_ambiguous_protein_column() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="column 'protein' is ambiguous and unsupported",
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                    "protein": ["P28482-2"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )


def test_normalizer_rejects_multiple_alias_matches_for_same_target() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="has ambiguous columns for 'gene_symbol': gene_symbol, gene",
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "gene": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )


def test_normalizer_derives_missing_gene_symbol_and_site_from_index() -> None:
    normalized = DatasetConventionNormalizer().run(
        phospho=_phospho(),
        site_metadata=pd.DataFrame(
            {
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        sample_metadata=None,
        total=None,
    )
    site_metadata = normalized.site_metadata
    assert site_metadata.loc["MAPK14;Y182;", "gene_symbol"] == "MAPK14"
    assert site_metadata.loc["MAPK14;Y182;", "site"] == "Y182"


def test_normalizer_fails_when_index_derivation_convention_is_not_met() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="site_metadata is missing required metadata columns",
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )
