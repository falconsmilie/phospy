from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.errors import ReferenceValidationError, UnsupportedInputFormatError
from phospy.transformations.models import TransformationState


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
            transformation_state=TransformationState.raw(has_total_matrix=False),
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
                transformation_state=TransformationState.raw(has_total_matrix=False),
            )
        )


def test_reference_bundle_rejects_ambiguous_site_sequence_ids() -> None:
    with pytest.raises(
        ReferenceValidationError,
        match="duplicate site identifiers after canonicalization",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {"kinase": ["MAP2K6"], "substrate_site": ["101"]}
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31, "B" * 31]},
                index=pd.Index([101, " 101 "], name="site_id"),
            ),
        )


def test_kinase_workflow_stably_matches_mixed_identifier_types() -> None:
    dataset = AnalysisReadyPhosphoDataset(
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
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": [101, " 202 "],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "B" * 31]},
            index=pd.Index(["101", 202], name="site_id"),
        ),
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=2),
            activity_config=None,
        )
    )

    assert list(result.scoring_result.profile_scores.index) == ["101", "202"]
    assert list(result.prediction_result.pred_mat.index) == ["101", "202"]
    assert list(result.references.site_sequences.index) == ["101", "202"]
    assert set(result.references.kinase_substrate_map.loc[:, "substrate_site"]) == {
        "101",
        "202",
    }
