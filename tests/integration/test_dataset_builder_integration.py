from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetTotalProteinCorrectionConfig,
    Organism,
    ReferencePreset,
)
from phospy.errors import DatasetValidationError, PhosPyInputError
from phospy.references.resolution import ReferenceResolver
from phospy.transformations.models import TransformationState
from tests.support.rewrite_fixture_data import load_rat_l6_phospho, site_metadata_for

pytestmark = pytest.mark.integration


def test_dataset_builder_builds_analysis_ready_dataset_from_fixture() -> None:
    phospho = load_rat_l6_phospho().head(32).copy(deep=True)
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            organism=Organism.RAT,
        )
    )
    pdt.assert_frame_equal(built.phospho, phospho)
    assert list(built.site_metadata.columns) == ["gene_symbol", "site", "site_sequence"]
    assert built.transformation_state == TransformationState.established_raw(
        has_total_matrix=False
    )


def test_dataset_builder_establishes_transformation_state_via_supported_path() -> None:
    phospho = load_rat_l6_phospho().head(8).copy(deep=True)
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            organism=Organism.RAT,
        )
    )
    assert built.transformation_state.label == "linear"
    assert built.transformation_state.is_established
    assert built.transformation_state.established_via is not None
    assert (
        built.transformation_state.phospho.established_by
        == "phospy.transformations.transformers.identity"
    )


def test_dataset_builder_preserves_total_matrix_and_establishes_linear_state() -> None:
    phospho = load_rat_l6_phospho().head(8).copy(deep=True)
    total = pd.DataFrame(
        {
            sample_name: [float(i + 1), float(i + 2)]
            for i, sample_name in enumerate(phospho.columns.astype(str))
        },
        index=pd.Index(["MAPK14", "AKT1"], name="protein_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            total=total,
            organism=Organism.RAT,
        )
    )
    assert built.total is not None
    pdt.assert_frame_equal(built.total, total)
    assert built.transformation_state == TransformationState.established_raw(
        has_total_matrix=True
    )
    assert built.transformation_state.label == "linear"


def test_dataset_builder_applies_total_protein_correction_when_requested() -> None:
    phospho = load_rat_l6_phospho().head(6).copy(deep=True)
    site_metadata = site_metadata_for(phospho)
    gene_symbols = site_metadata.loc[:, "gene_symbol"].astype(str)
    unique_genes = pd.Index(
        dict.fromkeys(gene_symbols.tolist()).keys(), name="protein_id"
    )
    total = pd.DataFrame(
        {
            sample_name: [
                float(gene_position + sample_position + 1)
                for gene_position in range(len(unique_genes))
            ]
            for sample_position, sample_name in enumerate(phospho.columns.astype(str))
        },
        index=unique_genes,
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="ratio_to_total"
                )
            ),
        )
    )

    total_by_site = total.reindex(gene_symbols.tolist())
    total_by_site.index = phospho.index
    expected = phospho - total_by_site
    pdt.assert_frame_equal(built.phospho, expected)
    pdt.assert_frame_equal(built.total, total)
    assert built.transformation_state == TransformationState.established_raw(
        has_total_matrix=True
    )


def test_dataset_builder_requires_total_when_ratio_correction_is_requested() -> None:
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    with pytest.raises(
        PhosPyInputError,
        match="policy='ratio_to_total' requires total input data",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata_for(phospho),
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    total_protein_correction=DatasetTotalProteinCorrectionConfig(
                        policy="ratio_to_total"
                    )
                ),
            )
        )


def test_dataset_builder_supports_documented_alias_and_index_derivation_conventions() -> (
    None
):
    phospho = load_rat_l6_phospho().head(8).copy(deep=True)
    canonical_site_metadata = site_metadata_for(phospho)
    site_metadata = pd.DataFrame(
        {
            "centralized_sequence": canonical_site_metadata.loc[
                :, "site_sequence"
            ].tolist(),
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    assert list(built.phospho.index) == list(phospho.index)
    assert list(built.site_metadata.columns) == ["site_sequence", "gene_symbol", "site"]
    assert (
        built.site_metadata.loc[:, "site_sequence"].tolist()
        == canonical_site_metadata.loc[:, "site_sequence"].tolist()
    )


def test_dataset_builder_preserves_explicit_protein_identity_column() -> None:
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    site_metadata = site_metadata_for(phospho).copy(deep=True)
    site_metadata.loc[:, "protein_id"] = [
        f"PROT_{position:03d}" for position in range(site_metadata.shape[0])
    ]
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    assert "protein_id" in built.site_metadata.columns
    assert (
        built.site_metadata.loc[:, "protein_id"].tolist()
        == site_metadata.loc[:, "protein_id"].tolist()
    )


def test_dataset_builder_supports_row_median_missing_data_preprocessing_policy() -> (
    None
):
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    phospho.iloc[0, 0] = float("nan")
    phospho.iloc[1, :] = float("nan")
    original_index = phospho.index.copy()
    expected_imputed = phospho.loc[original_index[0]].median(skipna=True)

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata_for(phospho),
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                ),
            ),
        )
    )

    assert built.phospho.index.tolist() == [
        site_id for site_id in original_index.tolist() if site_id != original_index[1]
    ]
    assert built.phospho.isna().to_numpy().sum() == 0
    assert built.phospho.loc[original_index[0], phospho.columns[0]] == expected_imputed


def test_dataset_builder_default_forbid_policy_keeps_missingness_strict() -> None:
    phospho = load_rat_l6_phospho().head(4).copy(deep=True)
    phospho.iloc[0, 0] = float("nan")

    with pytest.raises(DatasetValidationError, match="must not contain missing values"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata_for(phospho),
                organism=Organism.RAT,
            )
        )


def test_reference_bundle_rat_tables_are_structurally_coherent() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )
    substrate_sites = set(
        bundle.kinase_substrate_map.loc[:, "substrate_site"].astype(str)
    )
    known_sites = set(bundle.site_sequences.index.astype(str))
    assert substrate_sites.issubset(known_sites)
