from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
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
from phospy.io.bundles._shared.processing_state import (
    processing_state_from_payload,
    processing_state_to_payload,
)
from phospy.provenance.hashing import hash_table, hash_table_exact, hash_table_tolerance
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset
from tests.support.signalome_config import build_signalome_config

pytestmark = [pytest.mark.reproducibility, pytest.mark.release_gate]


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
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
            "protein_id": ["P28482", "Q9Y243"],
            "localisation_confidence": [0.95, 0.9],
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
                    "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
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


def test_dataset_output_fingerprints_match_observed_numeric_outputs() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )

    output_fingerprint = next(
        item
        for item in built.provenance.output_tables
        if item.name == "dataset.phospho"
    )

    assert output_fingerprint.hash_value == hash_table(
        built.phospho,
        name="dataset.phospho",
    )
    assert output_fingerprint.exact_hash_value == hash_table_exact(
        built.phospho,
        name="dataset.phospho",
    )
    assert output_fingerprint.tolerance_hash_value == hash_table_tolerance(
        built.phospho,
        name="dataset.phospho",
    )


def test_dataset_provenance_exact_hash_changes_for_tiny_shift_but_tolerance_hash_can_stay_stable() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.123456781, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    changed = phospho.copy(deep=True)
    changed.at["MAPK14;Y182;", "sample_a"] = 1.123456784
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=phospho.index.copy(),
    )

    first = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    second = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=changed,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    first_fingerprint = next(
        item
        for item in first.provenance.output_tables
        if item.name == "dataset.phospho"
    )
    second_fingerprint = next(
        item
        for item in second.provenance.output_tables
        if item.name == "dataset.phospho"
    )

    assert first_fingerprint.exact_hash_value != second_fingerprint.exact_hash_value
    assert (
        first_fingerprint.tolerance_hash_value
        == second_fingerprint.tolerance_hash_value
    )


def test_processing_state_bundle_round_trip_from_real_preprocessing_output() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan")],
            "sample_b": [2.0, 3.0],
            "sample_c": [3.0, 5.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=phospho.index.copy(),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                ),
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
            ),
        )
    )

    payload = processing_state_to_payload(built.processing_state)
    restored = processing_state_from_payload(payload)

    assert processing_state_to_payload(restored) == payload


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
            "localisation_confidence": [0.95, 0.9, 0.92],
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
        ScientificPolicyId.SIGNALOME_DOWNSTREAM_SCORE_SELECTION,
        ScientificPolicyId.SIGNALOME_CANDIDATE_SCORING,
        ScientificPolicyId.SIGNALOME_MISSING_VALUE_CLUSTERING,
        ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
        ScientificPolicyId.SIGNALOME_ASSIGNMENT_POLICY,
        ScientificPolicyId.SIGNALOME_NETWORK_POLICY,
        ScientificPolicyId.PROTEIN_MODULE_FROM_SITE_MEMBERSHIP,
    }.issubset(signalome_policy_ids)


def test_provenance_policy_metadata_includes_stable_name_and_version() -> None:
    dataset = _dataset_for_workflows()
    references = _reference_bundle()

    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            activity_config=KinaseActivityConfig(enabled=False),
        )
    )
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

    observed_names = set()
    for provenance in (kinase_result.provenance, signalome_result.provenance):
        assert provenance is not None
        for policy in provenance.scientific_policies:
            assert policy.name
            assert policy.version == "1"
            observed_names.add(policy.name)
    assert {
        "profile_correlation_v1",
        "rank_weighted_motif_profile_fusion_v1",
        "signalome_downstream_score_rank_weighted_preferred_v1",
        "signalome_candidate_scoring_full_v1",
        "score_preconditioning_error_on_drop_v1",
        "protein_module_from_site_membership_v1",
    }.issubset(observed_names)


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


def test_kinase_provenance_includes_duplicate_site_resolution_policy_when_preprocessed() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 1.2, 2.0],
            "sample_b": [1.5, 1.7, 2.5],
        },
        index=pd.Index(["r1", "r2", "r3"], name="row_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "GSK3B"],
            "site": ["Y182", "Y182", "S9"],
            "site_sequence": ["SEQ_A", "SEQ_A", "SEQ_B"],
            "protein_id": ["P28482", "P28482", "Q9Y243"],
            "localisation_confidence": [0.95, 0.94, 0.9],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="forbid",
                    min_observed_values=None,
                ),
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_policy="aggregate_mean",
                    missing_data_policy="drop_any_missing",
                ),
            ),
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K_DUP", "K_DUP"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["SEQ_A", "SEQ_B"]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            activity_config=KinaseActivityConfig(enabled=False),
        )
    )

    duplicate_policies = [
        policy
        for policy in kinase_result.provenance.scientific_policies
        if policy.id == ScientificPolicyId.DUPLICATE_SITE_RESOLUTION
    ]
    assert duplicate_policies
    assert duplicate_policies[0].name == "duplicate_site_resolution_aggregate_mean_v1"


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
            assert fingerprint.exact_hash_algorithm == "sha256-stable-json-v1"
            assert len(str(fingerprint.exact_hash_value)) == 64
            assert fingerprint.tolerance_hash_algorithm == "sha256-float-round-8dp-v1"
            assert len(str(fingerprint.tolerance_hash_value)) == 64
        for policy in provenance.scientific_policies:
            assert policy.name
            assert policy.version
