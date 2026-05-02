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
                "site": ["Y182"],
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
        match="column 'sequence' is unsupported",
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
        match="column 'protein' is unsupported",
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
        match="has ambiguous columns for 'site_sequence': site_sequence, centralized_sequence",
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                    "centralized_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )


def test_normalizer_rejects_unsupported_historical_site_alias() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="column 'residue' is unsupported",
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "residue": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )


def test_normalizer_rejects_conflicting_site_id_column_and_index() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="conflicting site identifiers between index and 'site_id' column",
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "site_id": ["AKT1;T308;"],
                    "gene_symbol": ["MAPK14"],
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


def test_normalizer_requires_exact_index_derivation_convention() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="site identifiers must use 'GENE;SITE;' format",
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;EXTRA"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )
