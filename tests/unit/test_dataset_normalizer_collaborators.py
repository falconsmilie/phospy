from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors.input import UnsupportedInputFormatError
from phospy.science.datasets.builders.normalization_reporter import (
    DatasetConventionNormalisationReporter,
)
from phospy.science.datasets.builders.sample_metadata_normalizer import (
    SampleMetadataNormalizer,
)
from phospy.science.datasets.builders.site_metadata_normalizer import (
    SiteIdentityFieldNormalizer,
    SiteMetadataColumnAliasResolver,
    SiteMetadataIndexNormalizer,
    SiteMetadataNormalizer,
)
from phospy.science.datasets.builders.total_matrix_normalizer import (
    TotalProteinMatrixNormalizer,
)


def test_site_metadata_column_alias_resolver_maps_supported_aliases() -> None:
    frame = pd.DataFrame(
        {
            "gene_name": ["MAPK14"],
            "site": ["Y182"],
            "centralized_sequence": ["SEQ_A"],
            "localization_probability": [0.9],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    normalized = SiteMetadataColumnAliasResolver().run(frame)

    assert "gene_symbol" in normalized.columns
    assert "site_sequence" in normalized.columns
    assert "localisation_confidence" in normalized.columns


def test_site_metadata_column_alias_resolver_rejects_duplicate_alias_matches() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="ambiguous columns for 'site_sequence'",
    ):
        SiteMetadataColumnAliasResolver().run(
            pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["SEQ_A"],
                    "centralized_sequence": ["SEQ_B"],
                }
            )
        )


def test_site_metadata_index_normalizer_uses_site_id_column_when_index_is_range() -> (
    None
):
    records = []
    normalized = SiteMetadataIndexNormalizer().run(
        pd.DataFrame(
            {
                "site_id": [" mapk14 ; y182 "],
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
            }
        ),
        site_identifier_records=records,
    )

    assert normalized.index.tolist() == ["MAPK14;Y182;"]
    assert normalized.index.name == "site_id"
    assert len(records) == 1


def test_site_metadata_normalizer_derives_gene_and_site_from_index() -> None:
    normalized = SiteMetadataNormalizer().run(
        pd.DataFrame(
            {"site_sequence": ["SEQ_A"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        phospho_index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        site_identifier_records=[],
    )

    assert normalized.loc["MAPK14;Y182;", "gene_symbol"] == "MAPK14"
    assert normalized.loc["MAPK14;Y182;", "site"] == "Y182"


def test_site_metadata_normalizer_rejects_impossible_gene_site_derivation() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="missing required metadata columns",
    ):
        SiteMetadataNormalizer().run(
            pd.DataFrame(
                {"site_sequence": ["SEQ_A"]},
                index=pd.Index(["not_a_site_id"], name="site_id"),
            ),
            phospho_index=pd.Index(["not_a_site_id"], name="site_id"),
            site_identifier_records=[],
        )


def test_sample_metadata_normalizer_trims_and_aligns_to_phospho_columns() -> None:
    normalized = SampleMetadataNormalizer().run(
        pd.DataFrame(
            {"group": ["A", "B"]},
            index=pd.Index([" sample_b ", " sample_a "], name="sample_id"),
        ),
        phospho_columns=pd.Index(["sample_a", "sample_b"], name="sample_id"),
    )

    assert normalized is not None
    assert normalized.index.tolist() == ["sample_a", "sample_b"]


def test_total_protein_matrix_normalizer_aligns_columns_to_phospho_columns() -> None:
    normalized = TotalProteinMatrixNormalizer().run(
        pd.DataFrame(
            {" sample_b ": [2.0], " sample_a ": [1.0]},
            index=pd.Index(["P28482"], name="protein_id"),
        ),
        phospho_columns=pd.Index(["sample_a", "sample_b"], name="sample_id"),
    )

    assert normalized is not None
    assert normalized.columns.tolist() == ["sample_a", "sample_b"]


def test_normalisation_reporter_builds_payload_for_changed_site_identifiers() -> None:
    reporter = DatasetConventionNormalisationReporter()
    records = []
    normalized = reporter.normalize_supported_site_index_if_present(
        pd.Index([" mapk14 ; y182 "], name="site_id"),
        field_name="dataset build request phospho.index",
        site_identifier_records=records,
    )

    assert normalized.tolist() == ["MAPK14;Y182;"]
    report = reporter.build_site_identifier_report(records)
    assert report is not None
    payload = report.to_payload()
    assert payload["changed_identifier_count"] == 1
    assert payload["records"][0]["field_name"] == "dataset build request phospho.index"
    assert payload["records"][0]["normalised_value"] == "MAPK14;Y182;"


def test_site_identity_field_normalizer_canonicalizes_gene_and_site_tokens() -> None:
    normalized = SiteIdentityFieldNormalizer().run(
        pd.DataFrame(
            {"gene_symbol": [" mapk14 "], "site": [" y182 "]},
            index=pd.Index(["row_1"]),
        )
    )

    assert normalized.loc["row_1", "gene_symbol"] == "MAPK14"
    assert normalized.loc["row_1", "site"] == "Y182"
