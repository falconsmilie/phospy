from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
)
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.api.configs import (
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_PREDICTION_MODES,
)
from phospy.provenance.serialization import from_payload, to_payload
from phospy.scientific_policies import ScientificPolicyId
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


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
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
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


def _run_workflow(
    *,
    workflow: KinaseWorkflow,
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferenceBundle,
    scoring_config: KinaseScoringConfig,
    prediction_mode: str,
):
    return workflow.run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=scoring_config,
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
                mode=prediction_mode,
                n_iterations=2,
                random_state=7,
            ),
            activity_config=None,
        )
    )


@pytest.mark.parametrize(
    "include_diagnostics",
    [False, True],
    ids=["without_diagnostics", "with_diagnostics"],
)
def test_supported_prediction_modes_preserve_scoring_stage_semantics(
    include_diagnostics: bool,
) -> None:
    workflow = KinaseWorkflow()
    dataset = _dataset()
    references = _references()
    scoring_config = KinaseScoringConfig(
        min_substrates=2,
        include_diagnostic_scoring_tables=include_diagnostics,
    )
    baseline = _run_workflow(
        workflow=workflow,
        dataset=dataset,
        references=references,
        scoring_config=scoring_config,
        prediction_mode=KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    )

    for prediction_mode in sorted(KINASE_PREDICTION_MODES):
        result = _run_workflow(
            workflow=workflow,
            dataset=dataset,
            references=references,
            scoring_config=scoring_config,
            prediction_mode=prediction_mode,
        )
        pd.testing.assert_frame_equal(
            result.scoring_result.profile_scores,
            baseline.scoring_result.profile_scores,
        )
        assert result.scoring_result.rank_weighted_fusion_scores is not None
        assert baseline.scoring_result.rank_weighted_fusion_scores is not None
        pd.testing.assert_frame_equal(
            result.scoring_result.rank_weighted_fusion_scores,
            baseline.scoring_result.rank_weighted_fusion_scores,
        )
        if include_diagnostics:
            assert result.scoring_result.motif_scores is not None
            assert baseline.scoring_result.motif_scores is not None
            pd.testing.assert_frame_equal(
                result.scoring_result.motif_scores,
                baseline.scoring_result.motif_scores,
            )
            assert result.scoring_result.score_fusion_weights is not None
            assert baseline.scoring_result.score_fusion_weights is not None
            pd.testing.assert_frame_equal(
                result.scoring_result.score_fusion_weights,
                baseline.scoring_result.score_fusion_weights,
            )
        else:
            assert result.scoring_result.motif_scores is None
            assert result.scoring_result.score_fusion_weights is None


def test_prediction_modes_keep_distinct_mode_specific_size_semantics() -> None:
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
                deterministic_max_selected_kinases=1,
                adaptive_ensemble_runs=99,
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
                deterministic_max_selected_kinases=1,
                adaptive_ensemble_runs=11,
                mode="adaptive_ensemble",
                n_iterations=2,
                random_state=7,
            ),
            activity_config=None,
        )
    )

    assert deterministic.provenance is not None
    assert adaptive.provenance is not None
    deterministic_prediction_config = deterministic.provenance.workflow_parameters[
        "prediction_config"
    ]
    adaptive_prediction_config = adaptive.provenance.workflow_parameters[
        "prediction_config"
    ]
    assert "ensemble_size" not in deterministic_prediction_config
    assert "ensemble_size" not in adaptive_prediction_config
    assert deterministic_prediction_config["deterministic_max_selected_kinases"] == 1
    assert deterministic_prediction_config["adaptive_ensemble_runs"] == 99
    assert adaptive_prediction_config["deterministic_max_selected_kinases"] == 1
    assert adaptive_prediction_config["adaptive_ensemble_runs"] == 11

    assert deterministic.prediction_result.pred_mat.shape[1] == 1
    assert adaptive.prediction_result.pred_mat.shape[1] >= 2


def test_prediction_config_docs_name_both_supported_lanes() -> None:
    doc = KinasePredictionConfig.__doc__ or ""
    assert "deterministic_ranking" in doc
    assert "adaptive_ensemble" in doc


def test_kinase_provenance_records_active_scientific_policies() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
                mode=KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
            ),
            activity_config=None,
        )
    )
    assert result.provenance is not None

    policy_ids = {policy.id for policy in result.provenance.scientific_policies}
    assert ScientificPolicyId.PROFILE_CORRELATION_SHIFTED_UNIT in policy_ids
    assert ScientificPolicyId.KINASE_PROFILE_SCORING in policy_ids
    assert ScientificPolicyId.MOTIF_PROFILE_RANK_FUSION in policy_ids
    assert ScientificPolicyId.CANDIDATE_SUBSTRATE_SELECTION in policy_ids
    assert ScientificPolicyId.SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY not in policy_ids
    candidate_policy = next(
        policy
        for policy in result.provenance.scientific_policies
        if policy.id == ScientificPolicyId.CANDIDATE_SUBSTRATE_SELECTION
    )
    assert candidate_policy.parameters["top_k"] == 2
    assert candidate_policy.parameters["score_threshold"] == pytest.approx(0.0)
    assert candidate_policy.parameters["inclusion"] == 1

    payload = to_payload(result.provenance)
    restored = from_payload(payload)
    restored_ids = {policy.id for policy in restored.scientific_policies}
    assert restored_ids == policy_ids


def test_adaptive_prediction_provenance_records_exact_seed() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
                mode="adaptive_ensemble",
                n_iterations=2,
                random_state=31,
            ),
            activity_config=None,
        )
    )
    assert result.provenance is not None
    prediction_config = result.provenance.workflow_parameters["prediction_config"]
    assert prediction_config["mode"] == "adaptive_ensemble"
    assert prediction_config["random_state"] == 31
    assert result.provenance.random_state == 31
