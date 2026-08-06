from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    Organism,
    ReferenceBundle,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.activities.membership import (
    ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE,
    ActivityMembershipSelection,
)
from phospy.science.activities.method_contracts import (
    ACTIVITY_SITE_UNIVERSE_KSEA_BACKGROUND,
    ACTIVITY_SITE_UNIVERSE_PREDICTED_MEMBERSHIP,
    ACTIVITY_SITE_UNIVERSE_REFERENCE_SUPPORTED_MEMBERSHIP,
    ACTIVITY_SITE_UNIVERSE_SSGSEA_EFFECT_RANKING,
    kinase_activity_method_universe_contract,
)
from phospy.science.activities.methods import SsgseaSubstrateEnrichmentActivityMethod
from phospy.science.activities.methods.ksea_zscore import (
    KSEA_STATUS_COMPUTED,
    KseaZScoreActivityMethod,
)
from phospy.science.activities.models import (
    KinaseActivityInputs,
    PredMatOverlapSummary,
)
from phospy.science.activities.semantics import ActivityInputMatrix
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.models import KinasePredictionResult
from phospy.workflows._row_attrition import make_row_attrition_record
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseSiteUniverses,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.executor import KinaseWorkflowExecutor
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def test_resolved_request_separates_measured_activity_and_scoring_universes() -> None:
    dataset = _dataset_with_three_measured_sites()
    scoring_site_index = dataset.phospho.index[:2].copy()
    site_sequences = _site_sequences_for(dataset, scoring_site_index)
    request = ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=_references_with_three_members(),
        kinase_substrate_map=_projected_membership(dataset),
        site_sequences=site_sequences,
        site_identity_map=_site_identity_map(dataset),
        scoring_site_index=scoring_site_index,
        activity_phospho_matrix=dataset.phospho.copy(deep=True),
        execution_config=_execution_config(),
    )

    assert request.activity_phospho_matrix.index.equals(dataset.phospho.index)
    assert request.scoring_phospho_matrix.index.equals(scoring_site_index)
    assert not request.activity_phospho_matrix.index.equals(request.scoring_site_index)
    assert request.site_universes is not None
    assert request.site_universes.ksea_background_sites.equals(dataset.phospho.index)
    assert set(
        request.scoring_kinase_substrate_map.loc[:, "substrate_site"].astype(str)
    ) == set(scoring_site_index.astype(str))
    assert set(
        request.reference_membership_map.loc[:, "substrate_site"].astype(str)
    ) == set(dataset.phospho.index.astype(str))


def test_resolved_request_rejects_accidental_measured_scoring_universe_substitution() -> (
    None
):
    dataset = _dataset_with_three_measured_sites()
    scoring_site_index = dataset.phospho.index[:2].copy()
    site_sequences = _site_sequences_for(dataset, scoring_site_index)
    substituted_universes = ResolvedKinaseSiteUniverses(
        measured_quantitative_sites=scoring_site_index.copy(),
        sequence_supported_scoring_sites=scoring_site_index.copy(),
        reference_supported_membership_sites=scoring_site_index.copy(),
        predicted_membership_sites=scoring_site_index.copy(),
        ksea_background_sites=scoring_site_index.copy(),
        ssgsea_effect_ranking_sites=scoring_site_index.copy(),
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="measured_universe_alignment",
    ):
        ResolvedKinaseWorkflowRequest(
            dataset=dataset,
            references=_references_with_three_members(),
            kinase_substrate_map=_projected_membership(dataset),
            site_sequences=site_sequences,
            site_identity_map=_site_identity_map(dataset),
            scoring_site_index=scoring_site_index,
            activity_phospho_matrix=dataset.phospho.copy(deep=True),
            execution_config=_execution_config(),
            site_universes=substituted_universes,
        )


def test_ksea_background_uses_full_quantitative_universe_not_prediction_rows() -> None:
    sites = _site_index()
    phospho = pd.DataFrame({"profile_a": [1.0, 2.0, 100.0]}, index=sites)
    pred_mat = pd.DataFrame({"K1": [0.9, 0.8]}, index=sites[:2])
    inputs = KinaseActivityInputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho,
        threshold=0.5,
        min_substrates=2,
        top_n_substrates=1,
        overlap_summary=PredMatOverlapSummary(
            overlap_count=2,
            pred_mat_rows=2,
            phospho_rows=3,
        ),
        activity_input=ActivityInputMatrix.sample_level_abundance(
            phospho,
            _assume_owned=True,
        ),
        membership_selection=ActivityMembershipSelection(
            source_category=ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE,
            selection_method="test_fixed_membership",
            selection_method_version="1",
            score_source="test_reference",
            consumed_tested_matrix=False,
            selected_kinase_universe=("K1",),
            selected_substrate_universe=tuple(sites[:2].astype(str).tolist()),
        ),
    )

    result = KseaZScoreActivityMethod(
        evidence_threshold=0.5,
        min_substrates=2,
    ).run(inputs)

    statistics = result.statistics_table
    assert statistics is not None
    row = statistics.iloc[0]
    assert row["computability_status"] == KSEA_STATUS_COMPUTED
    assert int(row["n_substrates"]) == 2
    assert int(row["n_background_sites"]) == 3


def test_ssgsea_intersects_reference_membership_with_effect_universe_explicitly() -> (
    None
):
    sites = _site_index()
    effect_matrix = pd.DataFrame({"contrast_a": [3.0, 2.0, -1.0]}, index=sites)
    membership = pd.DataFrame(
        {
            "kinase": ["K1", "K1", "K1", "K2"],
            "substrate_site": [
                str(sites[0]),
                str(sites[1]),
                "OUTSIDE;S4;",
                "OUTSIDE;S4;",
            ],
        }
    )

    result = SsgseaSubstrateEnrichmentActivityMethod(
        min_substrates=2,
        permutation_count=0,
    ).run(
        activity_input=ActivityInputMatrix.standardised_effect(
            effect_matrix,
            _assume_owned=True,
        ),
        kinase_substrate_membership=membership,
    )

    assert int(result.target_counts.loc["K1"]) == 2
    assert int(result.target_counts.loc["K2"]) == 0
    assert "OUTSIDE;S4;" not in set(result.target_table.loc[:, "site_id"].astype(str))
    statistics = result.statistics_table
    assert statistics is not None
    assert (
        int(statistics.loc[statistics["kinase"] == "K1", "n_background_sites"].iloc[0])
        == 3
    )


def test_activity_method_universe_contracts_declare_activity_owned_universes() -> None:
    ksea = kinase_activity_method_universe_contract("ksea_zscore")
    ssgsea = kinase_activity_method_universe_contract("ssgsea_substrate_enrichment")

    assert ksea.quantitative_universe == ACTIVITY_SITE_UNIVERSE_KSEA_BACKGROUND
    assert ksea.membership_universe == ACTIVITY_SITE_UNIVERSE_PREDICTED_MEMBERSHIP
    assert ksea.requires_sequence_context is False
    assert ssgsea.quantitative_universe == ACTIVITY_SITE_UNIVERSE_SSGSEA_EFFECT_RANKING
    assert (
        ssgsea.membership_universe
        == ACTIVITY_SITE_UNIVERSE_REFERENCE_SUPPORTED_MEMBERSHIP
    )
    assert ssgsea.requires_sequence_context is False


def test_ksea_fixed_membership_is_stable_when_scoring_sequence_support_changes() -> (
    None
):
    dataset = _dataset_with_three_measured_sites()
    references = _references_with_three_members()
    activity = ResolvedKinaseActivityExecutionConfig(
        method="ksea_zscore",
        threshold=0.5,
        min_substrates=2,
        top_n_substrates=2,
        ksea_min_substrates=2,
        ksea_evidence_threshold=0.5,
        ksea_p_value_method="normal_approximation",
        ksea_adjust_p_values=True,
    )
    full_sequence_request = _resolved_request_for(
        dataset=dataset,
        references=references,
        scoring_site_index=dataset.phospho.index.copy(),
        predicted_membership_sites=dataset.phospho.index.copy(),
        activity=activity,
    )
    missing_sequence_request = _resolved_request_for(
        dataset=dataset,
        references=references,
        scoring_site_index=dataset.phospho.index[:2].copy(),
        predicted_membership_sites=dataset.phospho.index.copy(),
        activity=activity,
    )
    fixed_prediction = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"K_REF": [0.9, 0.8, 0.7]},
            index=dataset.phospho.index.copy(),
        )
    )

    full_sequence_result = KinaseWorkflowExecutor()._run_activity_stage(
        request=full_sequence_request,
        config=full_sequence_request.execution_config,
        prediction_result=fixed_prediction,
    )
    missing_sequence_result = KinaseWorkflowExecutor()._run_activity_stage(
        request=missing_sequence_request,
        config=missing_sequence_request.execution_config,
        prediction_result=fixed_prediction,
    )

    assert full_sequence_result is not None
    assert missing_sequence_result is not None
    pd.testing.assert_frame_equal(
        full_sequence_result.activity_matrix,
        missing_sequence_result.activity_matrix,
    )
    full_stats = full_sequence_result.statistics_table
    missing_stats = missing_sequence_result.statistics_table
    assert full_stats is not None
    assert missing_stats is not None
    assert int(full_stats["n_background_sites"].iloc[0]) == 3
    assert int(missing_stats["n_background_sites"].iloc[0]) == 3


def test_workflow_serializes_separate_universe_attrition_records() -> None:
    dataset = _dataset_with_three_measured_sites()
    scoring_site_index = dataset.phospho.index[:2].copy()
    resolved = ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=_references_with_three_members(),
        kinase_substrate_map=_projected_membership(dataset),
        site_sequences=_site_sequences_for(dataset, scoring_site_index),
        site_identity_map=_site_identity_map(dataset),
        scoring_site_index=scoring_site_index,
        activity_phospho_matrix=dataset.phospho.copy(deep=True),
        row_attrition_records=tuple(
            record
            for record in (
                make_row_attrition_record(
                    workflow="kinase",
                    stage="kinase_sequence_context",
                    reason="sites_missing_valid_centered_sequence",
                    input_site_ids=dataset.phospho.index,
                    output_site_ids=scoring_site_index,
                ),
            )
            if record is not None
        ),
        execution_config=_execution_config(
            activity=ResolvedKinaseActivityExecutionConfig(
                method="ksea_zscore",
                threshold=0.5,
                min_substrates=2,
                top_n_substrates=2,
                ksea_min_substrates=2,
                ksea_evidence_threshold=0.0,
                ksea_p_value_method="normal_approximation",
                ksea_adjust_p_values=True,
            )
        ),
    )
    result = KinaseWorkflowExecutor().run(resolved)

    workflow_parameters = result.provenance.workflow_parameters
    universe_attrition = workflow_parameters["universe_attrition"]
    assert universe_attrition["sequence_attrition"][0]["removed_sites"] == 1
    assert universe_attrition["membership_attrition"]
    assert universe_attrition["finite_value_attrition"]
    background_records = universe_attrition["activity_background_attrition"]
    assert background_records[0]["output_universe"] == "ksea_background_sites"
    assert background_records[0]["output_sites"] == 3
    activity_config = workflow_parameters["activity_config"]
    assert activity_config["method_universe_contract"]["background_universe"] == (
        "ksea_background_sites"
    )


def _resolved_request_for(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferenceBundle,
    scoring_site_index: pd.Index,
    predicted_membership_sites: pd.Index,
    activity: ResolvedKinaseActivityExecutionConfig,
) -> ResolvedKinaseWorkflowRequest:
    return ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        kinase_substrate_map=_projected_membership(dataset),
        site_sequences=_site_sequences_for(dataset, scoring_site_index),
        site_identity_map=_site_identity_map(dataset),
        scoring_site_index=scoring_site_index,
        activity_phospho_matrix=dataset.phospho.copy(deep=True),
        execution_config=_execution_config(activity=activity),
        site_universes=ResolvedKinaseSiteUniverses(
            measured_quantitative_sites=dataset.phospho.index.copy(),
            sequence_supported_scoring_sites=scoring_site_index.copy(),
            reference_supported_membership_sites=dataset.phospho.index.copy(),
            predicted_membership_sites=predicted_membership_sites.copy(),
            ksea_background_sites=dataset.phospho.index.copy(),
            ssgsea_effect_ranking_sites=dataset.phospho.index.copy(),
        ),
    )


def _dataset_with_three_measured_sites() -> AnalysisReadyPhosphoDataset:
    sites = _site_index()
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {
                "profile_a": [1.0, 2.0, 4.0],
                "profile_b": [4.0, 2.0, 1.0],
            },
            index=sites.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": sites.astype(str).tolist(),
                "display_id": _display_ids(),
                **site_key_context_columns(sites),
                "gene_symbol": ["GENE1", "GENE2", "GENE3"],
                "protein_id": ["GENE1", "GENE2", "GENE3"],
                "site": ["S10", "T20", "S30"],
                "site_sequence": [_window("S"), _window("T"), _window("S")],
            },
            index=sites.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _references_with_three_members() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K_REF", "K_REF", "K_REF"],
                "substrate_site": _display_ids(),
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [_window("S"), _window("T"), _window("S")]},
            index=pd.Index(_display_ids(), name="site_id"),
        ),
    )


def _projected_membership(dataset: AnalysisReadyPhosphoDataset) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kinase": ["K_REF", "K_REF", "K_REF"],
            "substrate_site": dataset.phospho.index.astype(str).tolist(),
            "display_id": _display_ids(),
        }
    )


def _site_sequences_for(
    dataset: AnalysisReadyPhosphoDataset,
    scoring_site_index: pd.Index,
) -> pd.DataFrame:
    identity = _site_identity_map(dataset)
    display_ids = identity.loc[scoring_site_index, "display_id"].astype(str).tolist()
    residues = [
        display_id.split(";")[1].strip().upper()[0] for display_id in display_ids
    ]
    return pd.DataFrame(
        {
            "site_sequence": [_window(residue) for residue in residues],
            "display_id": display_ids,
        },
        index=scoring_site_index.copy(),
    )


def _site_identity_map(dataset: AnalysisReadyPhosphoDataset) -> pd.DataFrame:
    sites = dataset.phospho.index.astype(str).tolist()
    return pd.DataFrame(
        {
            "site_key": sites,
            "display_id": dataset.site_metadata.loc[:, "display_id"]
            .astype(str)
            .tolist(),
        },
        index=dataset.phospho.index.copy(),
    )


def _execution_config(
    *,
    activity: ResolvedKinaseActivityExecutionConfig | None = None,
) -> ResolvedKinaseExecutionConfig:
    return ResolvedKinaseExecutionConfig(
        scoring_min_substrates=2,
        include_diagnostic_scoring_tables=False,
        profile_missing_value_strategy="strict",
        prediction_top_k=2,
        prediction_deterministic_max_selected_kinases=1,
        prediction_adaptive_ensemble_runs=2,
        prediction_mode="deterministic_ranking",
        prediction_adaptive_policy="stable",
        prediction_n_iterations=5,
        prediction_random_state=None,
        activity=activity,
    )


def _site_index() -> pd.Index:
    return site_key_index_from_display_ids(_display_ids())


def _display_ids() -> tuple[str, str, str]:
    return ("GENE1;S10;", "GENE2;T20;", "GENE3;S30;")


def _window(residue: str) -> str:
    return ("A" * 15) + residue + ("A" * 15)
