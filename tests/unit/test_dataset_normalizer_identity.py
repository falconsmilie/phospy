from __future__ import annotations

import pandas as pd

from phospy.datasets.builders.normalizer import DatasetConventionNormalizer


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def test_normalizer_keeps_gene_symbol_and_protein_id_distinct() -> None:
    normalized = DatasetConventionNormalizer().run(
        phospho=_phospho(),
        site_metadata=pd.DataFrame(
            {
                "gene": ["MAPK14"],
                "protein": ["P28482-2"],
                "residue": ["Y182"],
                "sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        sample_metadata=None,
        total=None,
    )
    site_metadata = normalized.site_metadata
    assert site_metadata.loc["MAPK14;Y182;", "gene_symbol"] == "MAPK14"
    assert site_metadata.loc["MAPK14;Y182;", "protein_id"] == "P28482-2"
    assert (
        site_metadata.loc["MAPK14;Y182;", "gene_symbol"]
        != site_metadata.loc["MAPK14;Y182;", "protein_id"]
    )


def test_normalizer_does_not_backfill_gene_symbol_from_protein_id() -> None:
    normalized = DatasetConventionNormalizer().run(
        phospho=_phospho(),
        site_metadata=pd.DataFrame(
            {
                "protein_id": ["P28482-2"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        sample_metadata=None,
        total=None,
    )
    site_metadata = normalized.site_metadata
    assert site_metadata.loc["MAPK14;Y182;", "protein_id"] == "P28482-2"
    assert site_metadata.loc["MAPK14;Y182;", "gene_symbol"] == "MAPK14"
    assert (
        site_metadata.loc["MAPK14;Y182;", "gene_symbol"]
        != site_metadata.loc["MAPK14;Y182;", "protein_id"]
    )
