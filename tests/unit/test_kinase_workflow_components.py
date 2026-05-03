from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import Organism, ReferenceBundle
from phospy.api.configs import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
)
from phospy.api.results import (
    KinaseWorkflowPreprocessingAttritionSummary,
    KinaseWorkflowScoringAttritionSummary,
    KinaseWorkflowSiteAttritionSummary,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.workflows.kinase.activity_runner import KinaseActivityRunner
from phospy.workflows.kinase.component_models import KinaseScoringRunResult
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.prediction_runner import KinasePredictionRunner
from phospy.workflows.kinase.provenance import KinaseProvenanceBuilder
from phospy.workflows.kinase.result_assembly import KinaseResultAssembler
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id")
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0, 0.8],
                "sample_b": [2.0, 1.2],
            },
            index=site_ids,
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "GSK3B"],
                "site": ["Y182", "S9"],
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                ],
            },
            index=site_ids,
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )


def _config(
    *,
    prediction_mode: str = KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    activity: ResolvedKinaseActivityExecutionConfig | None = None,
) -> ResolvedKinaseExecutionConfig:
    return ResolvedKinaseExecutionConfig(
        scoring_min_substrates=2,
        include_diagnostic_scoring_tables=False,
        profile_missing_value_strategy="strict",
        prediction_top_k=2,
        prediction_deterministic_max_selected_kinases=2,
        prediction_adaptive_ensemble_runs=2,
        prediction_mode=prediction_mode,
        prediction_adaptive_policy="stable",
        prediction_n_iterations=3,
        prediction_random_state=7,
        activity=activity,
    )


def _activity_config(
    *,
    method: str,
) -> ResolvedKinaseActivityExecutionConfig:
    return ResolvedKinaseActivityExecutionConfig(
        method=method,  # type: ignore[arg-type]
        threshold=0.6,
        min_substrates=2,
        top_n_substrates=3,
        ksea_min_substrates=2,
        ksea_evidence_threshold=0.6,
        ksea_p_value_method="normal_approximation",
        ksea_adjust_p_values=True,
    )


def _resolved_request(
    *,
    config: ResolvedKinaseExecutionConfig | None = None,
) -> ResolvedKinaseWorkflowRequest:
    dataset = _dataset()
    references = _references()
    scoring_site_index = dataset.phospho.index.copy()
    return ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        kinase_substrate_map=references.kinase_substrate_map,
        site_sequences=references.site_sequences,
        scoring_site_index=scoring_site_index,
        activity_phospho_matrix=dataset.phospho.loc[scoring_site_index, :].copy(
            deep=True
        ),
        execution_config=config or _config(),
    )


def test_scoring_runner_returns_expected_downstream_score_source() -> None:
    request = _resolved_request()
    result = KinaseScoringRunner().run(
        request=request,
        config=request.execution_config,
    )
    assert result.downstream_score_source == "rank_weighted_fusion_scores"
    assert not result.downstream_score_matrix.empty


def test_deterministic_prediction_runner_handles_no_candidates() -> None:
    request = _resolved_request(config=_config())
    scoring_execution = KinaseScoringRunResult(
        scoring_result=KinaseScoringResult._from_owned(
            profile_scores=pd.DataFrame(
                {"MAP2K6": [0.1, 0.2]},
                index=request.scoring_site_index.copy(),
            )
        ),
        downstream_score_matrix=pd.DataFrame(
            {"MAP2K6": [-0.1, -0.2]},
            index=request.scoring_site_index.copy(),
        ),
        downstream_score_source="profile_scores",
        quantified_substrates={"MAP2K6": ["MAPK14;Y182;", "GSK3B;S9;"]},
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinasePredictionRunner().run(
            request=request,
            config=request.execution_config,
            scoring_execution=scoring_execution,
        )

    assert exc_info.value.seam == "kinase.executor.prediction_ensemble"


def test_adaptive_prediction_runner_handles_no_candidates() -> None:
    request = _resolved_request(
        config=_config(prediction_mode=KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE)
    )
    scoring_execution = KinaseScoringRunResult(
        scoring_result=KinaseScoringResult._from_owned(
            profile_scores=pd.DataFrame(
                {"MAP2K6": [0.1, 0.2]},
                index=request.scoring_site_index.copy(),
            )
        ),
        downstream_score_matrix=pd.DataFrame(
            {"MAP2K6": [-0.1, -0.2]},
            index=request.scoring_site_index.copy(),
        ),
        downstream_score_source="profile_scores",
        quantified_substrates={"MAP2K6": ["MAPK14;Y182;", "GSK3B;S9;"]},
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinasePredictionRunner().run(
            request=request,
            config=request.execution_config,
            scoring_execution=scoring_execution,
        )

    assert exc_info.value.seam == "kinase.executor.prediction_adaptive_candidates"


def test_activity_runner_returns_none_when_activity_is_disabled() -> None:
    request = _resolved_request(config=_config(activity=None))
    prediction_result = KinasePredictionResult._from_owned(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.7, 0.6]},
            index=request.scoring_site_index.copy(),
        )
    )
    assert (
        KinaseActivityRunner().run(
            request=request,
            config=request.execution_config,
            prediction_result=prediction_result,
        )
        is None
    )


def test_activity_runner_selects_weighted_method_from_config() -> None:
    request = _resolved_request(
        config=_config(
            activity=_activity_config(
                method=KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
            )
        )
    )
    prediction_result = KinasePredictionResult._from_owned(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.7, 0.6]},
            index=request.scoring_site_index.copy(),
        )
    )

    result = KinaseActivityRunner().run(
        request=request,
        config=request.execution_config,
        prediction_result=prediction_result,
    )
    assert result is not None
    assert result.activity_method.activity_method_id == (
        "simplified_weighted_substrate_activity_v1"
    )


def test_activity_runner_selects_ksea_method_from_config() -> None:
    request = _resolved_request(
        config=_config(
            activity=_activity_config(method=KINASE_ACTIVITY_METHOD_KSEA_ZSCORE)
        )
    )
    prediction_result = KinasePredictionResult._from_owned(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.7, 0.6]},
            index=request.scoring_site_index.copy(),
        )
    )

    result = KinaseActivityRunner().run(
        request=request,
        config=request.execution_config,
        prediction_result=prediction_result,
    )
    assert result is not None
    assert result.activity_method.activity_method_id == "ksea_zscore_v1"
    assert result.statistics_table is not None


def test_provenance_builder_includes_active_scientific_policies() -> None:
    activity = ResolvedKinaseActivityExecutionConfig(
        method="simplified_weighted_substrate_activity",
        threshold=0.6,
        min_substrates=2,
        top_n_substrates=3,
        ksea_min_substrates=5,
        ksea_evidence_threshold=0.6,
        ksea_p_value_method="normal_approximation",
        ksea_adjust_p_values=True,
    )
    request = _resolved_request(config=_config(activity=activity))
    scoring_result = KinaseScoringResult._from_owned(
        profile_scores=pd.DataFrame(
            {"MAP2K6": [0.7, 0.2]},
            index=request.scoring_site_index.copy(),
        )
    )
    prediction_result = KinasePredictionResult._from_owned(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.7, 0.2]},
            index=request.scoring_site_index.copy(),
        )
    )

    provenance = KinaseProvenanceBuilder().run(
        request=request,
        config=request.execution_config,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        activity_result=None,
    )
    policy_ids = {policy.id.value for policy in provenance.scientific_policies}
    assert "candidate_substrate_selection_v1" in policy_ids
    assert "simplified_weighted_substrate_activity_v1" in policy_ids


def test_provenance_builder_includes_ksea_policy_when_selected() -> None:
    activity = _activity_config(method=KINASE_ACTIVITY_METHOD_KSEA_ZSCORE)
    request = _resolved_request(config=_config(activity=activity))
    scoring_result = KinaseScoringResult._from_owned(
        profile_scores=pd.DataFrame(
            {"MAP2K6": [0.7, 0.2]},
            index=request.scoring_site_index.copy(),
        )
    )
    prediction_result = KinasePredictionResult._from_owned(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.7, 0.2]},
            index=request.scoring_site_index.copy(),
        )
    )

    provenance = KinaseProvenanceBuilder().run(
        request=request,
        config=request.execution_config,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        activity_result=None,
    )
    policy_ids = {policy.id.value for policy in provenance.scientific_policies}
    assert "ksea_zscore_activity_v1" in policy_ids


def test_result_assembler_preserves_owned_dataframe_transfer() -> None:
    request = _resolved_request()
    profile_scores = pd.DataFrame(
        {"MAP2K6": [0.8, 0.4]},
        index=request.scoring_site_index.copy(),
    )
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9, 0.3]},
        index=request.scoring_site_index.copy(),
    )
    substrate_list = pd.DataFrame(
        {
            "kinase": ["MAP2K6"],
            "substrate_site": ["MAPK14;Y182;"],
            "score": [0.9],
            "rank": [1],
        }
    )
    scoring_result = KinaseScoringResult._from_owned(profile_scores=profile_scores)
    prediction_result = KinasePredictionResult._from_owned(
        pred_mat=pred_mat,
        substrate_list=substrate_list,
    )
    provenance = KinaseProvenanceBuilder().run(
        request=request,
        config=request.execution_config,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        activity_result=None,
    )

    assembled = KinaseResultAssembler().run(
        request=request,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        site_attrition_summary=KinaseWorkflowSiteAttritionSummary(
            preprocessing=KinaseWorkflowPreprocessingAttritionSummary(
                input_rows=2,
                rows_removed_during_preprocessing=0,
                rows_removed_invalid_or_missing_site_identifiers=0,
                duplicate_sites_merged_or_resolved=0,
                output_rows=2,
            ),
            scoring=KinaseWorkflowScoringAttritionSummary(
                rows_removed_invalid_or_missing_site_identifiers=0,
                final_quantitative_sites_entering_scoring=2,
                sites_with_valid_site_sequence=2,
                sites_without_usable_site_sequence=0,
                sites_eligible_for_motif_scoring=2,
                sites_with_kinase_substrate_reference_profile_evidence=2,
                sites_contributing_to_final_fused_prediction_scoring_output=2,
                sites_contributing_to_activity_scoring=None,
            ),
        ),
        activity_result=None,
        provenance=provenance,
    )

    assert assembled.scoring_result is scoring_result
    assert assembled.prediction_result is prediction_result
    assert assembled.site_attrition_summary is not None
    assert assembled.scoring_result._profile_scores is scoring_result._profile_scores
    assert assembled.prediction_result._pred_mat is prediction_result._pred_mat
    assert (
        assembled.prediction_result._substrate_list is prediction_result._substrate_list
    )
