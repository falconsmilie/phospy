from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import KinaseWorkflow, Organism, ReferenceBundle, ReferencePreset
from phospy.api.configs import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
    KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS,
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.requests import KinaseWorkflowRequest
from phospy.api.results import (
    KinaseEligibilityReport,
    KinaseWorkflowPreprocessingAttritionSummary,
    KinaseWorkflowScoringAttritionSummary,
    KinaseWorkflowSiteAttritionSummary,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.activities.methods import SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE
from phospy.science.activities.scientific_policies import (
    SSGSEA_PERMUTATION_RNG_SEED_POLICY,
    SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION,
)
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.science.transformations.models import QuantitativeMeaning
from phospy.tables.kinase import KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from phospy.workflows.kinase.activity_runner import KinaseActivityRunner
from phospy.workflows.kinase.component_models import KinaseScoringRunResult
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.executor import KinaseWorkflowExecutor
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.prediction_runner import KinasePredictionRunner
from phospy.workflows.kinase.provenance import KinaseProvenanceBuilder
from phospy.workflows.kinase.reference_projection import KinaseReferenceProjector
from phospy.workflows.kinase.resolved_validator import (
    ResolvedKinaseEligibilityValidator,
)
from phospy.workflows.kinase.result_assembly import KinaseResultAssembler
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    build_prediction_outputs,
)
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state_with_meaning,
    supported_log2_processing_state_with_meaning,
)
from tests.support.site_keys import site_key_context_columns


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
    return trusted_analysis_ready_dataset_from_tables(
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
                **site_key_context_columns(site_keys),
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


def _effect_dataset() -> AnalysisReadyPhosphoDataset:
    base = _dataset()
    return trusted_analysis_ready_dataset_from_tables(
        phospho=base.phospho,
        site_metadata=base.site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state_with_meaning(
            has_total_matrix=False,
            meaning=QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        ),
        processing_state=supported_log2_processing_state_with_meaning(
            has_total_matrix=False,
            meaning=QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        ),
    )


def _dataset_with_duplicate_display_ids() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "MAPK14;Y182;", "GSK3B;S9;"]
    protein_ids = ["MAPK14_PROTEIN_A", "MAPK14_PROTEIN_B", "GSK3B"]
    site_keys = [
        encode_site_key(
            build_protein_scoped_site_key(
                organism="rat",
                protein_namespace="protein_id",
                protein_identifier=protein_id,
                residue=display_id.split(";")[1][0],
                position=int(display_id.split(";")[1][1:]),
                field_name="test.dataset.duplicate_display.site_key",
                error_type=ValueError,
            )
        )
        for display_id, protein_id in zip(display_ids, protein_ids, strict=True)
    ]
    site_index = pd.Index(site_keys, name="site_key")
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0, 1.1, 0.8],
                "sample_b": [2.0, 2.1, 1.2],
            },
            index=site_index,
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_keys,
                "display_id": display_ids,
                **site_key_context_columns(site_keys),
                "gene_symbol": ["MAPK14", "MAPK14", "GSK3B"],
                "site": ["Y182", "Y182", "S9"],
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                ],
                "protein_id": protein_ids,
            },
            index=site_index,
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


def _allow_unknown_reference_context_scoring_config() -> KinaseScoringConfig:
    return KinaseScoringConfig(
        reliability_profile="custom",
        min_substrates=2,
        reference_context_compatibility_policy=(
            ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
        ),
    )


def _config(
    *,
    prediction_mode: str = KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    activity: ResolvedKinaseActivityExecutionConfig | None = None,
    include_substrate_contributions: bool = False,
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
        include_substrate_contributions=include_substrate_contributions,
    )


def _activity_config(
    *,
    method: str,
    ssgsea_min_substrates: int = 5,
    ssgsea_permutations: int = 0,
    ssgsea_random_seed: int | None = 0,
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
        ssgsea_min_substrates=ssgsea_min_substrates,
        ssgsea_ranking_direction="descending",
        ssgsea_permutations=ssgsea_permutations,
        ssgsea_random_seed=ssgsea_random_seed,
        ssgsea_adjust_p_values=True,
    )


def _resolved_request(
    *,
    config: ResolvedKinaseExecutionConfig | None = None,
    references: ReferenceBundle | None = None,
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> ResolvedKinaseWorkflowRequest:
    dataset = dataset or _dataset()
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


def test_kinase_workflow_calls_validator_interpreter_executor_in_order() -> None:
    events: list[str] = []
    validated = object()
    interpreted = object()
    expected_result = object()

    class _Validator:
        def run(self, request: object) -> object:
            events.append("validator")
            return validated

    class _Interpreter:
        def run(self, request: object) -> object:
            events.append("interpreter")
            assert request is validated
            return interpreted

    class _Executor:
        def run(self, request: object) -> object:
            events.append("executor")
            assert request is interpreted
            return expected_result

    result = KinaseWorkflow._with_components(
        validator=_Validator(),  # type: ignore[arg-type]
        interpreter=_Interpreter(),  # type: ignore[arg-type]
        executor=_Executor(),  # type: ignore[arg-type]
    ).run(object())  # type: ignore[arg-type]

    assert events == ["validator", "interpreter", "executor"]
    assert result is expected_result


def test_kinase_stage_components_expose_run() -> None:
    for stage_type in (
        KinaseWorkflow,
        KinaseWorkflowValidator,
        KinaseWorkflowInterpreter,
        ResolvedKinaseEligibilityValidator,
        KinaseWorkflowExecutor,
    ):
        assert callable(getattr(stage_type, "run", None))


def test_kinase_validator_preserves_reference_preset_until_interpretation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from phospy.science.references.resolution import ReferenceResolver

    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("reference resolution must not run during validation")

    monkeypatch.setattr(ReferenceResolver, "run", fail_if_called)

    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.RAT,
        scoring_config=KinaseScoringConfig.exploratory(),
    )
    validated = KinaseWorkflowValidator().run(request)

    assert validated is request
    assert validated.references is ReferencePreset.RAT


def test_kinase_interpreter_resolves_and_projects_references_for_execution() -> None:
    dataset = _dataset()
    resolved_references = _mixed_case_references()
    calls: list[tuple[object, object]] = []

    class _ReferenceResolverSpy:
        def run(self, reference_input: object, *, dataset_organism: object):
            calls.append((reference_input, dataset_organism))
            return resolved_references

    interpreted = KinaseWorkflowInterpreter(
        reference_resolver=_ReferenceResolverSpy(),  # type: ignore[arg-type]
    ).run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.RAT,
            scoring_config=_allow_unknown_reference_context_scoring_config(),
        )
    )

    assert calls == [(ReferencePreset.RAT, Organism.RAT)]
    assert interpreted.references is resolved_references
    assert set(interpreted.kinase_substrate_map.loc[:, "kinase"]) == {"MAP2K6"}
    assert set(interpreted.kinase_substrate_map.loc[:, "substrate_site"]) == set(
        dataset.phospho.index.astype(str)
    )
    assert set(interpreted.site_sequences.index.astype(str)) == set(
        interpreted.scoring_site_index.astype(str)
    )


def test_kinase_executor_consumes_resolved_request_without_reference_discovery() -> (
    None
):
    resolved = _resolved_request()
    events: list[str] = []
    expected_eligibility_report = object()
    expected_scoring_result = object()
    expected_prediction_result = object()
    expected_activity_result = object()
    expected_site_attrition_summary = object()
    expected_provenance = object()
    expected_result = object()

    class _ScoringExecution:
        def __init__(self, scoring_result: object) -> None:
            self.scoring_result = scoring_result

    expected_scoring_execution = _ScoringExecution(expected_scoring_result)

    class _EligibilityReportComposer:
        def run(self, *, request: object, config: object) -> object:
            events.append("eligibility")
            assert request is resolved
            assert config is resolved.execution_config
            return expected_eligibility_report

    class _ScoringRunner:
        def run(
            self,
            *,
            request: object,
            config: object,
            collect_substrate_contributions: bool,
        ) -> object:
            events.append("scoring")
            assert request is resolved
            assert config is resolved.execution_config
            assert collect_substrate_contributions is False
            return expected_scoring_execution

    class _PredictionRunner:
        def run(
            self,
            *,
            request: object,
            config: object,
            scoring_execution: object,
        ) -> object:
            events.append("prediction")
            assert request is resolved
            assert config is resolved.execution_config
            assert scoring_execution is expected_scoring_execution
            return expected_prediction_result

    class _ActivityRunner:
        def run(
            self,
            *,
            request: object,
            config: object,
            prediction_result: object,
        ) -> object:
            events.append("activity")
            assert request is resolved
            assert config is resolved.execution_config
            assert prediction_result is expected_prediction_result
            return expected_activity_result

    class _SiteAttritionSummaryComposer:
        def run(
            self,
            *,
            request: object,
            scoring_execution: object,
            prediction_result: object,
            activity_enabled: bool,
        ) -> object:
            events.append("site_attrition")
            assert request is resolved
            assert scoring_execution is expected_scoring_execution
            assert prediction_result is expected_prediction_result
            assert activity_enabled is True
            return expected_site_attrition_summary

    class _ProvenanceBuilder:
        def run(
            self,
            *,
            request: object,
            config: object,
            scoring_result: object,
            prediction_result: object,
            activity_result: object,
            substrate_contributions: object,
        ) -> object:
            events.append("provenance")
            assert request is resolved
            assert config is resolved.execution_config
            assert scoring_result is expected_scoring_result
            assert prediction_result is expected_prediction_result
            assert activity_result is expected_activity_result
            assert substrate_contributions is None
            return expected_provenance

    class _ResultAssembler:
        def run(
            self,
            *,
            request: object,
            scoring_result: object,
            prediction_result: object,
            eligibility_report: object,
            site_attrition_summary: object,
            activity_result: object,
            provenance: object,
            substrate_contributions: object,
        ) -> object:
            events.append("assembly")
            assert request is resolved
            assert scoring_result is expected_scoring_result
            assert prediction_result is expected_prediction_result
            assert eligibility_report is expected_eligibility_report
            assert site_attrition_summary is expected_site_attrition_summary
            assert activity_result is expected_activity_result
            assert provenance is expected_provenance
            assert substrate_contributions is None
            return expected_result

    result = KinaseWorkflowExecutor(
        eligibility_report_composer=_EligibilityReportComposer(),  # type: ignore[arg-type]
        scoring_runner=_ScoringRunner(),  # type: ignore[arg-type]
        prediction_runner=_PredictionRunner(),  # type: ignore[arg-type]
        activity_runner=_ActivityRunner(),  # type: ignore[arg-type]
        site_attrition_summary_composer=_SiteAttritionSummaryComposer(),  # type: ignore[arg-type]
        provenance_builder=_ProvenanceBuilder(),  # type: ignore[arg-type]
        result_assembler=_ResultAssembler(),  # type: ignore[arg-type]
    ).run(resolved)

    assert events == [
        "eligibility",
        "scoring",
        "prediction",
        "activity",
        "site_attrition",
        "provenance",
        "assembly",
    ]
    assert result is expected_result


def test_kinase_executor_scoring_stage_threads_contribution_collection_flag() -> None:
    resolved = _resolved_request(config=_config(include_substrate_contributions=True))
    expected_scoring_execution = object()
    observed_collection_flag: bool | None = None

    class _ScoringRunner:
        def run(
            self,
            *,
            request: object,
            config: object,
            collect_substrate_contributions: bool,
        ) -> object:
            nonlocal observed_collection_flag
            assert request is resolved
            assert config is resolved.execution_config
            observed_collection_flag = collect_substrate_contributions
            return expected_scoring_execution

    result = KinaseWorkflowExecutor(
        scoring_runner=_ScoringRunner(),  # type: ignore[arg-type]
    )._run_scoring_stage(
        request=resolved,
        config=resolved.execution_config,
    )

    assert result is expected_scoring_execution
    assert observed_collection_flag is True


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


def test_execution_contract_rejects_display_indexed_scoring_site_index() -> None:
    dataset = _dataset()
    references = _references()
    kinase_substrate_map, site_sequences = _project_references_to_site_key(
        dataset,
        references,
    )
    display_index = pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id")

    with pytest.raises(WorkflowBoundaryError, match="scoring_site_index"):
        ResolvedKinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            site_identity_map=_site_identity_map(dataset),
            scoring_site_index=display_index,
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
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=_allow_unknown_reference_context_scoring_config(),
        )
    )
    scoring = KinaseScoringRunner().run(
        request=interpreted,
        config=interpreted.execution_config,
    )

    assert not scoring.scoring_result.profile_scores.empty
    assert scoring.downstream_score_source == "rank_weighted_fusion_scores"


def test_kinase_workflow_request_default_activity_config_does_not_emit_activity() -> (
    None
):
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
        scoring_config=_allow_unknown_reference_context_scoring_config(),
    )

    result = KinaseWorkflow().run(request)

    assert request.activity_config is None
    assert result.activity_result is None
    assert result.provenance is not None
    assert result.provenance.workflow_parameters["activity_config"] is None


def test_interpreter_preserves_reference_bundle_validation_report() -> None:
    interpreted = KinaseWorkflowInterpreter().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=_allow_unknown_reference_context_scoring_config(),
        )
    )

    report = interpreted.references.validation_report

    assert report.kinase_substrate_record_count == 2
    assert {item.table_name for item in report.required_tables} == {
        "kinase_substrate_map",
        "site_sequences",
    }


def test_interpreter_overlap_uses_normalised_reference_tables_after_bundle_construction() -> (
    None
):
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_mixed_case_references(),
        scoring_config=_allow_unknown_reference_context_scoring_config(),
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


def test_interpreter_maps_one_display_reference_to_multiple_site_keys() -> None:
    dataset = _dataset_with_duplicate_display_ids()

    interpreted = KinaseWorkflowInterpreter().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=_references(),
            scoring_config=_allow_unknown_reference_context_scoring_config(),
            reference_display_ambiguity_policy=(
                KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS
            ),
        )
    )

    duplicate_rows = interpreted.kinase_substrate_map.loc[
        interpreted.kinase_substrate_map.loc[:, "display_id"] == "MAPK14;Y182;",
        :,
    ]
    assert duplicate_rows.shape[0] == 2
    assert duplicate_rows.loc[:, "substrate_site"].nunique() == 2
    assert set(duplicate_rows.loc[:, "substrate_site"].astype(str)) == set(
        dataset.phospho.index.astype(str)[:2]
    )
    diagnostics = interpreted.site_sequence_merge_diagnostics[
        "display_reference_matching"
    ]
    assert diagnostics["reference_key"] == "display_id"
    assert diagnostics["dataset_row_identity"] == "site_key"
    assert diagnostics["one_to_many_display_reference_match_count"] == 1
    assert diagnostics["one_to_many_display_reference_site_key_rows"] == 2
    matches = diagnostics["one_to_many_display_reference_matches"]
    assert isinstance(matches, list)
    assert matches[0]["display_id"] == "MAPK14;Y182;"
    assert set(matches[0]["site_keys"]) == set(dataset.phospho.index.astype(str)[:2])
    assert matches[0]["matched_row_count"] == 2
    assert matches[0]["reference_row_count"] == 1


def test_duplicate_display_ids_are_accepted_before_reference_projection() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset_with_duplicate_display_ids(),
        references=_references(),
        scoring_config=_allow_unknown_reference_context_scoring_config(),
    )

    validated = KinaseWorkflowValidator().run(request)

    assert validated is request
    assert validated.reference_display_ambiguity_policy == (
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR
    )


def test_ambiguous_display_reference_projection_fails_by_default() -> None:
    dataset = _dataset_with_duplicate_display_ids()
    expected_site_keys = dataset.phospho.index.astype(str).tolist()[:2]

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflowInterpreter().run(
            KinaseWorkflowRequest(
                dataset=dataset,
                references=_references(),
                scoring_config=_allow_unknown_reference_context_scoring_config(),
            )
        )

    error = exc_info.value
    message = str(error)
    assert error.seam == "kinase.interpreter.reference_display_ambiguity"
    assert "MAPK14;Y182;" in message
    for site_key in expected_site_keys:
        assert site_key in message
    assert error.details["ambiguous_display_ids"] == ["MAPK14;Y182;"]
    diagnostics = error.details["ambiguity_diagnostics"]
    assert isinstance(diagnostics, list)
    assert diagnostics[0]["display_id"] == "MAPK14;Y182;"
    assert set(diagnostics[0]["site_keys"]) == set(expected_site_keys)
    assert diagnostics[0]["matched_row_count"] == 2
    assert diagnostics[0]["reference_row_count"] == 1
    assert diagnostics[0]["reference_rows"] == [
        {
            "row_position": 0,
            "row_index": "0",
            "kinase": "MAP2K6",
            "substrate_site": "MAPK14;Y182;",
        }
    ]


def test_projected_kinase_substrate_map_preserves_site_key_identity() -> None:
    dataset = _dataset_with_duplicate_display_ids()

    interpreted = KinaseWorkflowInterpreter().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=_references(),
            scoring_config=_allow_unknown_reference_context_scoring_config(),
            reference_display_ambiguity_policy=(
                KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS
            ),
        )
    )

    substrate_sites = interpreted.kinase_substrate_map.loc[:, "substrate_site"].astype(
        str
    )
    assert set(substrate_sites) <= set(dataset.phospho.index.astype(str))
    assert "MAPK14;Y182;" not in set(substrate_sites)
    assert {"substrate_site", "display_id"} <= set(interpreted.kinase_substrate_map)
    assert set(
        interpreted.kinase_substrate_map.loc[
            interpreted.kinase_substrate_map.loc[:, "display_id"] == "MAPK14;Y182;",
            "substrate_site",
        ].astype(str)
    ) == set(dataset.phospho.index.astype(str)[:2])


def test_one_to_many_reference_matching_diagnostic_includes_display_and_site_keys() -> (
    None
):
    dataset = _dataset_with_duplicate_display_ids()

    interpreted = KinaseWorkflowInterpreter().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=_references(),
            scoring_config=_allow_unknown_reference_context_scoring_config(),
            reference_display_ambiguity_policy=(
                KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS
            ),
        )
    )

    diagnostics = interpreted.site_sequence_merge_diagnostics[
        "display_reference_matching"
    ]
    matches = diagnostics["one_to_many_display_reference_matches"]
    assert isinstance(matches, list)
    assert matches == [
        {
            "display_id": "MAPK14;Y182;",
            "site_keys": dataset.phospho.index.astype(str).tolist()[:2],
            "matched_row_count": 2,
            "reference_row_count": 1,
            "reference_rows": [
                {
                    "row_position": 0,
                    "row_index": "0",
                    "kinase": "MAP2K6",
                    "substrate_site": "MAPK14;Y182;",
                }
            ],
            "projected_rows": 2,
            "interpreter_version": "phospy.workflows.kinase.reference_projector.v1",
        }
    ]
    assert diagnostics["ambiguity_policy"] == (
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS
    )


def test_site_key_based_reference_projection_is_not_ambiguous() -> None:
    dataset = _dataset_with_duplicate_display_ids()
    site_identity_map = _site_identity_map(dataset)
    direct_site_key = dataset.phospho.index.astype(str).tolist()[0]

    result = KinaseReferenceProjector().run(
        reference_kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": [direct_site_key]}
        ),
        site_identity_map=site_identity_map,
        ambiguity_policy=KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    )

    assert result.ambiguity_diagnostics == ()
    assert result.kinase_substrate_map.to_dict(orient="records") == [
        {
            "kinase": "MAP2K6",
            "substrate_site": direct_site_key,
            "display_id": "MAPK14;Y182;",
        }
    ]


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
    assert prediction_result.pred_mat.index.equals(request.scoring_site_index)
    substrate_list = prediction_result.substrate_list
    assert substrate_list is not None
    assert set(substrate_list.loc[:, "kinase"]) <= {"MAP2K6"}
    assert {"site_key", "display_id"} <= set(substrate_list)
    assert set(substrate_list.loc[:, "site_key"].astype(str)) <= set(
        request.scoring_site_index.astype(str)
    )


def test_duplicate_display_id_does_not_overwrite_prediction_rows() -> None:
    dataset = _dataset_with_duplicate_display_ids()
    request = KinaseWorkflowInterpreter().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=_references(),
            scoring_config=_allow_unknown_reference_context_scoring_config(),
            reference_display_ambiguity_policy=(
                KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS
            ),
        )
    )
    scores = pd.DataFrame(
        {"MAP2K6": [0.95, 0.85, 0.75]},
        index=request.scoring_site_index.copy(),
    )
    scoring_execution = KinaseScoringRunResult(
        scoring_result=KinaseScoringResult._from_owned(profile_scores=scores),
        downstream_score_matrix=scores,
        downstream_score_source="profile_scores",
        quantified_substrates={
            "MAP2K6": request.scoring_site_index.astype(str).tolist()
        },
    )

    prediction_result = KinasePredictionRunner().run(
        request=request,
        config=request.execution_config,
        scoring_execution=scoring_execution,
    )

    duplicate_site_keys = request.site_identity_map.loc[
        request.site_identity_map.loc[:, "display_id"] == "MAPK14;Y182;",
        "site_key",
    ].astype(str)
    observed_scores = prediction_result.pred_mat.loc[
        duplicate_site_keys.tolist(), "MAP2K6"
    ].astype(float)
    assert observed_scores.to_dict() == {
        duplicate_site_keys.iloc[0]: pytest.approx(0.95),
        duplicate_site_keys.iloc[1]: pytest.approx(0.85),
    }
    substrate_list = prediction_result.substrate_list
    assert substrate_list is not None
    duplicate_rows = substrate_list.loc[
        substrate_list.loc[:, "display_id"] == "MAPK14;Y182;", :
    ]
    assert set(duplicate_rows.loc[:, "site_key"].astype(str)) == set(
        duplicate_site_keys
    )


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
            {"MAP2K6": [0.9, 0.8]},
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
    assert {"site_key", "display_id"} <= set(result.target_table.columns)
    assert set(result.target_table.loc[:, "site_key"].astype(str)) == set(
        request.scoring_site_index.astype(str)
    )
    assert set(result.target_table.loc[:, "display_id"].astype(str)) == {
        "MAPK14;Y182;",
        "GSK3B;S9;",
    }


def test_activity_runner_selects_ksea_method_from_config() -> None:
    request = _resolved_request(
        config=_config(
            activity=_activity_config(method=KINASE_ACTIVITY_METHOD_KSEA_ZSCORE)
        ),
        dataset=_effect_dataset(),
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


def test_activity_runner_selects_ssgsea_method_from_config() -> None:
    request = _resolved_request(
        config=_config(
            activity=_activity_config(
                method=KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
                ssgsea_min_substrates=2,
            )
        ),
        dataset=_effect_dataset(),
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
        "ssgsea_substrate_enrichment_activity_v1"
    )
    assert result.statistics_table is not None
    assert {"site_key", "display_id"} <= set(result.target_table.columns)
    assert set(result.target_table.loc[:, "site_key"].astype(str)) == set(
        request.scoring_site_index.astype(str)
    )


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


def test_provenance_builder_includes_ssgsea_policy_when_selected() -> None:
    activity = _activity_config(
        method=KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
        ssgsea_min_substrates=2,
        ssgsea_permutations=5,
        ssgsea_random_seed=13,
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
    assert "ssgsea_substrate_enrichment_activity_v1" in policy_ids
    activity_payload = provenance.workflow_parameters["activity_config"]
    assert isinstance(activity_payload, Mapping)
    assert activity_payload["ssgsea_random_seed"] == 13
    assert activity_payload["ssgsea_permutation_rng_seed_policy"] == (
        SSGSEA_PERMUTATION_RNG_SEED_POLICY
    )
    assert activity_payload["ssgsea_permutation_rng_seed_policy_version"] == (
        SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION
    )
    assert activity_payload["ssgsea_significance_status"] == (
        SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE
    )
    assert activity_payload["ssgsea_significance_status_counts"] is None


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
    substrate_contributions = pd.DataFrame.from_records(
        [
            {
                "kinase": "MAP2K6",
                "substrate_site": str(request.scoring_site_index[0]),
                "substrate_identifier": "MAPK14;Y182;",
                "value_used_in_scoring": 0.9,
                "score_component": "rank_weighted_fusion_scores",
                "score_source": "profile_only_motif_missing_or_constant",
                "reference_source_name": "fixture",
                "reference_source_version": "v1",
                "reference_bundle_id": "fixture_bundle",
                "reference_identifier_namespace": "display_id",
                "status": "included",
                "exclusion_reason": None,
                "ambiguous": False,
            }
        ],
        columns=pd.Index(KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS),
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
        substrate_contributions=substrate_contributions,
    )

    assert assembled.scoring_result is scoring_result
    assert assembled.prediction_result is prediction_result
    assert assembled.site_attrition_summary is not None
    assert assembled.scoring_result._profile_scores is scoring_result._profile_scores
    assert assembled.prediction_result._pred_mat is prediction_result._pred_mat
    assert (
        assembled.prediction_result._substrate_list is prediction_result._substrate_list
    )
    assert assembled._substrate_contributions is substrate_contributions
