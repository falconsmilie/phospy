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
from phospy.api.requests import KinaseWorkflowRequest
from phospy.api.results import (
    KinaseEligibilityReport,
    KinaseWorkflowPreprocessingAttritionSummary,
    KinaseWorkflowScoringAttritionSummary,
    KinaseWorkflowSiteAttritionSummary,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from phospy.workflows.kinase.activity_runner import KinaseActivityRunner
from phospy.workflows.kinase.component_models import KinaseScoringRunResult
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.prediction_runner import KinasePredictionRunner
from phospy.workflows.kinase.provenance import KinaseProvenanceBuilder
from phospy.workflows.kinase.result_assembly import KinaseResultAssembler
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    build_prediction_outputs,
)
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;"]
    site_keys = [
        encode_site_key(
            build_protein_scoped_site_key(
                organism="rat",
                protein_namespace="protein_id",
                protein_identifier="MAPK14",
                residue="Y",
                position=182,
                field_name="test.dataset.site_key.mapk14",
                error_type=ValueError,
            )
        ),
        encode_site_key(
            build_protein_scoped_site_key(
                organism="rat",
                protein_namespace="protein_id",
                protein_identifier="GSK3B",
                residue="S",
                position=9,
                field_name="test.dataset.site_key.gsk3b",
                error_type=ValueError,
            )
        ),
    ]
    site_ids = pd.Index(site_keys, name="site_key")
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
                "site_key": site_keys,
                "display_id": display_ids,
                "gene_symbol": ["MAPK14", "GSK3B"],
                "site": ["Y182", "S9"],
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
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
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )


def _site_identity_map(dataset: AnalysisReadyPhosphoDataset) -> pd.DataFrame:
    metadata = dataset._borrow_site_metadata_frame().reindex(dataset.phospho.index)
    site_keys = dataset.phospho.index.astype(str).tolist()
    return pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": metadata.loc[:, "display_id"].astype(str).tolist(),
        },
        index=pd.Index(site_keys, name=dataset.phospho.index.name),
    )


def _project_references_to_site_key(
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferenceBundle,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    identity_map = _site_identity_map(dataset)
    display_lookup = {
        str(display_id): str(site_key)
        for site_key, display_id in identity_map.loc[
            :, ["site_key", "display_id"]
        ].itertuples(index=False)
    }
    site_sequence_rows: list[dict[str, str]] = []
    for site_key, display_id in identity_map.loc[
        :, ["site_key", "display_id"]
    ].itertuples(index=False):
        if display_id not in references.site_sequences.index:
            continue
        site_sequence_rows.append(
            {
                "site_key": str(site_key),
                "display_id": str(display_id),
                "site_sequence": str(
                    references.site_sequences.at[display_id, "site_sequence"]
                ),
            }
        )
    site_sequences = pd.DataFrame.from_records(site_sequence_rows).set_index("site_key")
    site_sequences.index.name = dataset.phospho.index.name

    map_rows: list[dict[str, str]] = []
    for kinase, display_id in references.kinase_substrate_map.loc[
        :, ["kinase", "substrate_site"]
    ].itertuples(index=False):
        site_key = display_lookup.get(str(display_id))
        if site_key is None:
            continue
        map_rows.append(
            {
                "kinase": str(kinase),
                "substrate_site": site_key,
                "display_id": str(display_id),
            }
        )
    kinase_substrate_map = pd.DataFrame.from_records(map_rows).drop_duplicates(
        ignore_index=True
    )
    return kinase_substrate_map, site_sequences


def _mixed_case_references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["map2k6", "Map2K6"],
                "substrate_site": [" mapk14 ; y182 ", "gsk3b;s9"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ]
            },
            index=pd.Index(["mapk14 ; y182", "GSK3B;S9"], name="site_id"),
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
    references: ReferenceBundle | None = None,
) -> ResolvedKinaseWorkflowRequest:
    dataset = _dataset()
    references = references or _references()
    kinase_substrate_map, site_sequences = _project_references_to_site_key(
        dataset,
        references,
    )
    scoring_site_index = site_sequences.index.copy()
    return ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        kinase_substrate_map=kinase_substrate_map,
        site_sequences=site_sequences,
        site_identity_map=_site_identity_map(dataset).loc[scoring_site_index],
        scoring_site_index=scoring_site_index,
        activity_phospho_matrix=dataset.phospho.loc[scoring_site_index].copy(deep=True),
        execution_config=config or _config(),
    )


def test_direct_unnormalised_kinase_ids_are_rejected_by_execution_contract() -> None:
    dataset = _dataset()
    references = _references()
    _, site_sequences = _project_references_to_site_key(dataset, references)

    with pytest.raises(
        WorkflowBoundaryError,
        match="kinase.contracts.kinase_substrate_map_schema",
    ):
        ResolvedKinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": [" map2k6 ", "MAP2K6"],
                    "substrate_site": dataset.phospho.index.astype(str).tolist(),
                }
            ),
            site_sequences=site_sequences,
            site_identity_map=_site_identity_map(dataset),
            scoring_site_index=dataset.phospho.index.copy(),
            activity_phospho_matrix=dataset.phospho.copy(deep=True),
            execution_config=_config(),
        )


def test_direct_unnormalised_site_ids_are_rejected_by_execution_contract() -> None:
    dataset = _dataset()
    references = _references()
    _, mapped_site_sequences = _project_references_to_site_key(dataset, references)
    unnormalised_site_sequences = mapped_site_sequences.copy(deep=True)
    unnormalised_site_sequences.index = pd.Index(
        [" phospy:v1|broken=mapk14 ", "phospy:v1|broken=gsk3b"],
        name=dataset.phospho.index.name,
    )

    with pytest.raises(
        WorkflowBoundaryError, match="kinase.contracts.site_sequence_schema"
    ):
        ResolvedKinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            kinase_substrate_map=_project_references_to_site_key(dataset, references)[
                0
            ],
            site_sequences=unnormalised_site_sequences,
            site_identity_map=_site_identity_map(dataset),
            scoring_site_index=dataset.phospho.index.copy(),
            activity_phospho_matrix=dataset.phospho.copy(deep=True),
            execution_config=_config(),
        )


def test_already_normalised_execution_inputs_are_accepted() -> None:
    request = _resolved_request()

    assert set(request.kinase_substrate_map.loc[:, "display_id"]) == {
        "MAPK14;Y182;",
        "GSK3B;S9;",
    }
    mapk14_site_key = request.site_identity_map.loc[
        request.site_identity_map.loc[:, "display_id"] == "MAPK14;Y182;", "site_key"
    ].iloc[0]
    assert request.site_sequences.loc[mapk14_site_key, "site_sequence"] == (
        "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"
    )


def test_workflow_execution_with_valid_reference_bundle_still_succeeds() -> None:
    interpreted = KinaseWorkflowInterpreter().run(
        KinaseWorkflowRequest(dataset=_dataset(), references=_references())
    )
    scoring = KinaseScoringRunner().run(
        request=interpreted,
        config=interpreted.execution_config,
    )

    assert not scoring.scoring_result.profile_scores.empty
    assert scoring.downstream_score_source == "rank_weighted_fusion_scores"


def test_interpreter_overlap_uses_normalised_reference_tables_after_bundle_construction() -> (
    None
):
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_mixed_case_references(),
    )

    interpreted = KinaseWorkflowInterpreter().run(request)
    assert set(interpreted.kinase_substrate_map.loc[:, "kinase"]) == {"MAP2K6"}
    assert set(interpreted.kinase_substrate_map.loc[:, "display_id"]) == {
        "MAPK14;Y182;",
        "GSK3B;S9;",
    }
    assert set(interpreted.site_sequences.loc[:, "display_id"]) == {
        "MAPK14;Y182;",
        "GSK3B;S9;",
    }
    overlap_sites = interpreted.scoring_site_index.intersection(
        interpreted.kinase_substrate_map.loc[:, "substrate_site"]
    )
    assert len(overlap_sites) == 2


def test_scoring_runner_returns_expected_downstream_score_source() -> None:
    request = _resolved_request()
    result = KinaseScoringRunner().run(
        request=request,
        config=request.execution_config,
    )
    assert result.downstream_score_source == "rank_weighted_fusion_scores"
    assert not result.downstream_score_matrix.empty


def test_scoring_runner_receives_normalised_reference_identifiers() -> None:
    request = _resolved_request(references=_mixed_case_references())
    captured_map: pd.DataFrame | None = None

    def _capture_build_profiles(
        *,
        phospho: pd.DataFrame,
        kinase_substrate_map: pd.DataFrame,
        min_substrates: int,
        allow_single_substrate_profiles: bool,
        profile_missing_value_strategy: str,
    ):
        nonlocal captured_map
        captured_map = kinase_substrate_map.copy(deep=True)
        return build_kinase_profiles(
            phospho=phospho,
            kinase_substrate_map=kinase_substrate_map,
            min_substrates=min_substrates,
            allow_single_substrate_profiles=allow_single_substrate_profiles,
            profile_missing_value_strategy=profile_missing_value_strategy,
        )

    KinaseScoringRunner(build_profiles=_capture_build_profiles).run(
        request=request,
        config=request.execution_config,
    )

    assert captured_map is not None
    assert set(captured_map.loc[:, "kinase"]) == {"MAP2K6"}
    assert set(captured_map.loc[:, "substrate_site"]) == set(
        request.scoring_site_index.astype(str)
    )


def test_prediction_output_construction_receives_normalised_kinase_ids() -> None:
    request = _resolved_request(references=_mixed_case_references())
    scoring_execution = KinaseScoringRunner().run(
        request=request,
        config=request.execution_config,
    )
    captured_selected_kinases: pd.Index | None = None
    captured_candidate_keys: set[str] | None = None

    def _capture_build_outputs(
        *,
        prediction_score_matrix: pd.DataFrame,
        selected_kinases: pd.Index,
        candidate_substrates: dict[str, list[str]],
        top_k: int,
        retain_full_scores: bool = False,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        nonlocal captured_selected_kinases, captured_candidate_keys
        captured_selected_kinases = selected_kinases.copy()
        captured_candidate_keys = set(candidate_substrates)
        return build_prediction_outputs(
            prediction_score_matrix=prediction_score_matrix,
            selected_kinases=selected_kinases,
            candidate_substrates=candidate_substrates,
            top_k=top_k,
            retain_full_scores=retain_full_scores,
        )

    prediction_result = KinasePredictionRunner(
        build_outputs=_capture_build_outputs
    ).run(
        request=request,
        config=request.execution_config,
        scoring_execution=scoring_execution,
    )

    assert captured_selected_kinases is not None
    assert captured_candidate_keys is not None
    assert set(captured_selected_kinases.astype(str)) == {"MAP2K6"}
    assert captured_candidate_keys == {"MAP2K6"}
    assert set(prediction_result.pred_mat.columns.astype(str)) == {"MAP2K6"}
    if prediction_result.substrate_list is not None:
        assert set(prediction_result.substrate_list.loc[:, "kinase"]) <= {"MAP2K6"}


def test_activity_scoring_receives_only_normalised_kinase_ids() -> None:
    request = _resolved_request(
        config=_config(
            activity=_activity_config(
                method=KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
            )
        ),
        references=_mixed_case_references(),
    )
    scoring_execution = KinaseScoringRunner().run(
        request=request,
        config=request.execution_config,
    )
    prediction_result = KinasePredictionRunner().run(
        request=request,
        config=request.execution_config,
        scoring_execution=scoring_execution,
    )

    class _ActivityValidatorSpy(KinaseActivityInputValidator):
        def run(self, **kwargs):
            pred_mat = kwargs["pred_mat"]
            assert set(pred_mat.columns.astype(str)) == {"MAP2K6"}
            return super().run(**kwargs)

    result = KinaseActivityRunner(activity_input_validator=_ActivityValidatorSpy()).run(
        request=request,
        config=request.execution_config,
        prediction_result=prediction_result,
    )
    assert result is not None


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
    # Intentional private-seam coverage: this protects zero-copy ownership transfer
    # across internal execution stages, a performance/copy-budget contract.
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
        eligibility_report=KinaseEligibilityReport(
            total_dataset_sites=2,
            sequence_complete_sites=2,
            localisation_eligible_sites=None,
            reference_overlap_sites=2,
            excluded_no_reference_match=0,
            excluded_low_localisation=None,
            eligible_kinases=1,
            excluded_kinases_below_min_substrates=0,
        ),
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
