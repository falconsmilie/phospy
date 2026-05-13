from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.input import UnsupportedInputFormatError
from phospy.science.datasets.builders.normalizer import DatasetConventionNormalizer


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


def test_normalizer_allows_residue_column_without_treating_it_as_site_alias() -> None:
    normalized = DatasetConventionNormalizer().run(
        phospho=_phospho(),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "residue": ["Y"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        sample_metadata=None,
        total=None,
    )
    assert normalized.site_metadata.loc["MAPK14;Y182;", "residue"] == "Y"


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


def test_normalizer_owns_copies_and_does_not_mutate_inputs() -> None:
    phospho = pd.DataFrame(
        {" sample_a ": [1.0], " sample_b ": [2.0]},
        index=pd.Index([" mapk14 ; y182 "], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_id": [" mapk14 ; y182 "],
            "gene_symbol": ["mapk14"],
            "site": ["y182"],
            "site_sequence": ["A" * 31],
        }
    )
    sample_metadata = pd.DataFrame(
        {"group": ["g1", "g1"]},
        index=pd.Index([" sample_a ", " sample_b "], name="sample_id"),
    )
    total = pd.DataFrame(
        {" sample_a ": [5.0], " sample_b ": [6.0]},
        index=pd.Index(["P28482-2"], name="protein_id"),
    )

    phospho_original = phospho.copy(deep=True)
    site_metadata_original = site_metadata.copy(deep=True)
    sample_metadata_original = sample_metadata.copy(deep=True)
    total_original = total.copy(deep=True)

    normalized = DatasetConventionNormalizer().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
    )

    pd.testing.assert_frame_equal(phospho, phospho_original)
    pd.testing.assert_frame_equal(site_metadata, site_metadata_original)
    pd.testing.assert_frame_equal(sample_metadata, sample_metadata_original)
    pd.testing.assert_frame_equal(total, total_original)

    assert list(normalized.phospho.index) == ["MAPK14;Y182;"]
    assert list(normalized.phospho.columns) == ["sample_a", "sample_b"]
    assert list(normalized.site_metadata.index) == ["MAPK14;Y182;"]
    assert list(normalized.sample_metadata.index) == ["sample_a", "sample_b"]
    assert list(normalized.total.columns) == ["sample_a", "sample_b"]


def test_normalizer_rejects_missing_sample_labels_before_stringification() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match=(
            "dataset build request sample_metadata.index must not contain missing "
            "labels"
        ),
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=pd.DataFrame(
                {"group": ["g1", "g2"]},
                index=pd.Index(["sample_a", pd.NA], name="sample_id"),
            ),
            total=None,
        )


def test_normalizer_rejects_missing_site_labels_before_stringification() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="dataset build request phospho.index must not contain missing labels",
    ):
        DatasetConventionNormalizer().run(
            phospho=pd.DataFrame(
                {"sample_a": [1.0]},
                index=pd.Index([pd.NA], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )


def test_normalizer_rejects_blank_labels_with_field_name() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match=("dataset build request phospho.columns must contain non-blank labels"),
    ):
        DatasetConventionNormalizer().run(
            phospho=pd.DataFrame(
                {"   ": [1.0]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=None,
            total=None,
        )


def test_normalizer_rejects_duplicate_labels_introduced_by_trimming() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match=(
            "dataset build request sample_metadata.index contains duplicate labels "
            "introduced by normalization"
        ),
    ):
        DatasetConventionNormalizer().run(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            sample_metadata=pd.DataFrame(
                {"group": ["g1", "g2"]},
                index=pd.Index(["sample_a", " sample_a "], name="sample_id"),
            ),
            total=None,
        )


def test_normalizer_rejects_duplicate_site_ids_introduced_by_canonicalization() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="contains duplicate site identifiers after canonicalization",
    ):
        DatasetConventionNormalizer().run(
            phospho=pd.DataFrame(
                {"sample_a": [1.0, 2.0]},
                index=pd.Index(
                    ["mapk14;y182", " MAPK14 ; Y182 ;"],
                    name="site_id",
                ),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "MAPK14"],
                    "site": ["Y182", "Y182"],
                    "site_sequence": ["SEQ_A", "SEQ_B"],
                },
                index=pd.Index(
                    ["mapk14;y182", " MAPK14 ; Y182 ;"],
                    name="site_id",
                ),
            ),
            sample_metadata=None,
            total=None,
        )
