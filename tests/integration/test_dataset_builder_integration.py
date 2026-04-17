from __future__ import annotations

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
    assert built.transformation_state == TransformationState.raw(has_total_matrix=False)


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
