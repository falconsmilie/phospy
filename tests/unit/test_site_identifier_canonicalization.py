from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    Organism,
    ReferenceBundle,
)
from phospy.errors import (
    DatasetValidationError,
    ReferenceValidationError,
    UnsupportedInputFormatError,
)
from phospy.references.models import ReferencePreset
from phospy.references.resolution import ReferenceResolver
from tests.support.transformation_states import supported_linear_state


def test_builder_canonicalizes_mixed_site_id_types_and_reorders_site_metadata() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index([101, " 202 "], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_id": ["202", 101],
            "gene_symbol": ["AKT1", "MAPK14"],
            "site": ["T308", "Y182"],
            "site_sequence": ["A" * 31, "B" * 31],
        }
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )

    assert list(built.phospho.index) == ["101", "202"]
    assert list(built.site_metadata.index) == ["101", "202"]
    assert built.site_metadata.loc["101", "gene_symbol"] == "MAPK14"
    assert built.site_metadata.loc["202", "gene_symbol"] == "AKT1"


def test_builder_rejects_ambiguous_site_ids_after_canonicalization() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0]},
        index=pd.Index([101, " 101 "], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["A" * 31, "B" * 31],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        UnsupportedInputFormatError,
        match="duplicate site identifiers after canonicalization",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
            )
        )


def test_reference_bundle_rejects_ambiguous_site_sequence_ids() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="contains colliding site identifiers when stripped",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31, "B" * 31]},
                index=pd.Index(
                    ["MAPK14;Y182;", " MAPK14;Y182;"],
                    name="site_id",
                ),
            ),
        )


def test_dataset_boundary_rejects_non_canonical_site_ids() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.phospho\\.index must contain canonical site identifiers",
    ):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0, 3.0],
                    "sample_b": [2.0, 4.0],
                },
                index=pd.Index([101, 202], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "AKT1"],
                    "site": ["Y182", "T308"],
                    "site_sequence": ["A" * 31, "B" * 31],
                },
                index=pd.Index([101, 202], name="site_id"),
            ),
            organism=Organism.RAT,
            transformation_state=supported_linear_state(has_total_matrix=False),
        )


def test_dataset_boundary_rejects_colliding_dirty_site_ids() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.phospho\\.index contains colliding site identifiers when stripped",
    ):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0, 3.0],
                    "sample_b": [2.0, 4.0],
                },
                index=pd.Index(
                    ["MAPK14;Y182;", " MAPK14;Y182;"],
                    name="site_id",
                ),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14", "MAPK14"],
                    "site": ["Y182", "Y182"],
                    "site_sequence": ["A" * 31, "B" * 31],
                },
                index=pd.Index(
                    ["MAPK14;Y182;", " MAPK14;Y182;"],
                    name="site_id",
                ),
            ),
            organism=Organism.RAT,
            transformation_state=supported_linear_state(has_total_matrix=False),
        )


def test_dataset_boundary_requires_explicit_transformation_state() -> None:
    with pytest.raises(TypeError, match="missing 1 required positional argument"):
        from phospy.datasets.models import AnalysisReadyPhosphoDataset

        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(
                {"sample_a": [1.0]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["A" * 31],
                },
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            organism=Organism.RAT,
        )  # type: ignore[call-arg]


def test_reference_bundle_rejects_duplicate_kinase_substrate_pairs() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="contains duplicate \\(kinase, substrate_site\\) pairs",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["MAP2K6", "MAP2K6"],
                    "substrate_site": ["MAPK14;Y182;", "MAPK14;Y182;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        )


def test_reference_bundle_rejects_colliding_dirty_substrate_site_ids() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match=(
            "references\\.kinase_substrate_map\\.substrate_site contains colliding "
            "site identifiers when stripped"
        ),
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["MAP2K6", "MAP2K7"],
                    "substrate_site": ["MAPK14;Y182;", " MAPK14;Y182;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        )


def test_reference_provider_shapes_bundled_resources_for_strict_boundary() -> None:
    bundle = ReferenceResolver().run(
        ReferencePreset.RAT,
        dataset_organism=Organism.RAT,
    )

    assert bundle.site_sequences.index.is_unique
    assert (
        bundle.kinase_substrate_map.duplicated(
            subset=["kinase", "substrate_site"]
        ).sum()
        == 0
    )
    assert all(
        isinstance(value, str) and value == value.strip() and value != ""
        for value in bundle.site_sequences.index.tolist()
    )
    assert all(
        isinstance(value, str) and value == value.strip() and value != ""
        for value in bundle.kinase_substrate_map.loc[:, "substrate_site"].tolist()
    )
