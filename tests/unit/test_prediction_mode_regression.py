from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyPhosphoDataset,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.transformations.models import TransformationState


def _dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = pd.Index(
        ["GENEA;S1;", "GENEA;S2;", "GENEB;S3;", "GENEB;S4;"],
        name="site_id",
    )
    phospho = pd.DataFrame(
        {
            "sample_a": [5.0, 4.5, 1.0, 1.2],
            "sample_b": [4.7, 5.1, 1.3, 1.0],
            "sample_c": [1.1, 1.0, 4.8, 4.9],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["GENEA", "GENEA", "GENEB", "GENEB"],
            "site": ["S1", "S2", "S3", "S4"],
            "site_sequence": ["A" * 31, "B" * 31, "C" * 31, "D" * 31],
        },
        index=site_ids,
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        transformation_state=TransformationState.established_raw(
            has_total_matrix=False
        ),
    )


def _references() -> ReferenceBundle:
    site_ids = pd.Index(
        ["GENEA;S1;", "GENEA;S2;", "GENEB;S3;", "GENEB;S4;"],
        name="site_id",
    )
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KINASE_A", "KINASE_A", "KINASE_B", "KINASE_B"],
                "substrate_site": [
                    "GENEA;S1;",
                    "GENEA;S2;",
                    "GENEB;S3;",
                    "GENEB;S4;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "B" * 31, "C" * 31, "D" * 31]},
            index=site_ids,
        ),
    )


def test_prediction_modes_keep_distinct_ensemble_size_semantics() -> None:
    workflow = KinaseWorkflow()
    dataset = _dataset()
    references = _references()
    deterministic = workflow.run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                ensemble_size=1,
                mode="deterministic_ranking",
            ),
            activity_config=None,
        )
    )
    adaptive = workflow.run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                ensemble_size=1,
                mode="adaptive_ensemble",
                n_iterations=2,
                random_state=7,
            ),
            activity_config=None,
        )
    )

    assert deterministic.prediction_result.pred_mat.shape[1] == 1
    assert adaptive.prediction_result.pred_mat.shape[1] >= 2


def test_prediction_config_docs_name_both_supported_lanes() -> None:
    doc = KinasePredictionConfig.__doc__ or ""
    assert "deterministic_ranking" in doc
    assert "adaptive_ensemble" in doc
