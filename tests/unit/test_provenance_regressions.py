from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.api import (
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinasePredictionResult,
    KinaseScoringConfig,
    KinaseScoringResult,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
)
from phospy.provenance.hashing import hash_table
from phospy.scientific_policies import ScientificPolicyId
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset
from tests.support.signalome_config import build_signalome_config


def _base_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=pd.Index(["A;S1;", "B;S2;"], name="site_id"),
    )


def _dataset_for_workflows() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 0.8],
            "sample_b": [2.0, 1.2],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
            "protein_id": ["P28482", "Q9Y243"],
        },
        index=phospho.index.copy(),
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


def _reference_bundle() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1", "K2", "K2"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                ],
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


def test_table_hash_changes_when_values_change() -> None:
    table = _base_table()
    changed = table.copy(deep=True)
    changed.at["A;S1;", "sample_a"] = 99.0
    assert hash_table(table, name="dataset.phospho") != hash_table(
        changed, name="dataset.phospho"
    )


def test_table_hash_changes_when_dtypes_change() -> None:
    table = _base_table()
    changed = table.astype({"sample_b": "int64"})
    assert hash_table(table, name="dataset.phospho") != hash_table(
        changed, name="dataset.phospho"
    )


def test_table_hash_changes_when_column_label_type_changes() -> None:
    numeric_label = pd.DataFrame(
        {1: [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    string_label = pd.DataFrame(
        {"1": [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    assert hash_table(numeric_label, name="table") != hash_table(
        string_label, name="table"
    )


def test_table_hash_changes_when_axis_names_change() -> None:
    first = pd.DataFrame(
        {"sample_a": [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    second = first.copy(deep=True)
    second.index = second.index.rename("site_id")
    assert hash_table(first, name="table") != hash_table(second, name="table")

    third = first.copy(deep=True)
    third.columns = third.columns.rename("sample_id")
    assert hash_table(first, name="table") != hash_table(third, name="table")


def test_preprocessing_policy_changes_are_visible_in_provenance() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 2.5, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )

    default_result = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    imputed_result = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                )
            ),
        )
    )

    default_stage = next(
        stage
        for stage in default_result.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )
    imputed_stage = next(
        stage
        for stage in imputed_result.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )

    assert default_stage.operation == "forbid"
    assert imputed_stage.operation == "impute_row_median"
    assert default_stage.parameters != imputed_stage.parameters


def test_active_scientific_policy_ids_are_present_in_kinase_and_signalome_provenance() -> (
    None
):
    dataset = _dataset_for_workflows()
    references = _reference_bundle()

    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            activity_config=KinaseActivityConfig(enabled=False),
        )
    )
    kinase_policy_ids = {
        policy.id for policy in kinase_result.provenance.scientific_policies
    }
    assert {
        ScientificPolicyId.PROFILE_CORRELATION_SHIFTED_UNIT,
        ScientificPolicyId.KINASE_PROFILE_SCORING,
        ScientificPolicyId.MOTIF_PROFILE_RANK_FUSION,
        ScientificPolicyId.CANDIDATE_SUBSTRATE_SELECTION,
    }.issubset(kinase_policy_ids)

    signalome_kinase_result = KinaseWorkflowResult(
        dataset=dataset,
        references=references,
        scoring_result=KinaseScoringResult(
            profile_scores=pd.DataFrame(
                {"K1": [0.1, 0.9], "K2": [0.9, 0.1]},
                index=dataset.phospho.index.copy(),
            ),
            rank_weighted_fusion_scores=pd.DataFrame(
                {"K1": [0.1, 0.9], "K2": [0.9, 0.1]},
                index=dataset.phospho.index.copy(),
            ),
        ),
        prediction_result=KinasePredictionResult(
            pred_mat=pd.DataFrame(
                {"K1": [0.9, 0.8], "K2": [0.8, 0.9]},
                index=dataset.phospho.index.copy(),
            )
        ),
        activity_result=None,
    )

    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=signalome_kinase_result,
            config=build_signalome_config(module_count=1),
        )
    )
    signalome_policy_ids = {
        policy.id for policy in signalome_result.provenance.scientific_policies
    }
    assert {
        ScientificPolicyId.SIGNALOME_MODULE_CANDIDATE_SCORE,
        ScientificPolicyId.SIGNALOME_MISSING_VALUE_CLUSTERING,
        ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
        ScientificPolicyId.PROTEIN_MODULE_FROM_SITE_MEMBERSHIP,
    }.issubset(signalome_policy_ids)


def test_adaptive_kinase_prediction_records_random_seed_provenance() -> None:
    dataset = build_rat_l6_dataset(n_sites=120)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                mode="adaptive_ensemble",
                top_k=4,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
                n_iterations=2,
                adaptive_policy="stable",
                random_state=17,
            ),
            activity_config=None,
        )
    )
    provenance = result.provenance
    assert provenance is not None
    assert provenance.random_state == 17
    assert provenance.random_seed_policy == "stable_by_kinase"
    assert provenance.workflow_parameters["prediction_config"]["random_state"] == 17


def test_workflow_provenance_fingerprints_and_policy_versions_are_stable() -> None:
    dataset = build_rat_l6_dataset(n_sites=180)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=4,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
            ),
            activity_config=None,
        )
    )
    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(module_count=1),
        )
    )

    for provenance in (kinase_result.provenance, signalome_result.provenance):
        assert provenance is not None
        assert provenance.workflow_name in {"kinase_workflow", "signalome_workflow"}
        assert provenance.input_tables
        assert provenance.output_tables
        for fingerprint in (*provenance.input_tables, *provenance.output_tables):
            assert fingerprint.hash_algorithm == "sha256"
            assert len(fingerprint.hash_value) == 64
        for policy in provenance.scientific_policies:
            assert policy.name
            assert policy.version
