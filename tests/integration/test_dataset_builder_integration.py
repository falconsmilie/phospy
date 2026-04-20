from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    Organism,
    ReferencePreset,
)
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
    assert list(built.phospho.index) == list(phospho.index)
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
