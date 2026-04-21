from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.api.requests import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.results import KinaseWorkflowResult
from phospy.datasets.builders.validator import DatasetBuildRequestValidator
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors import (
    PhosPyInputError,
    ReferenceCompatibilityError,
    UnsupportedInputFormatError,
    WorkflowValidationError,
)
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.transformations.models import TransformationState
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator


def test_kinase_scoring_default_sets_two_substrate_support_floor() -> None:
    assert KinaseScoringConfig().min_substrates == 2
    assert KinaseScoringConfig().include_diagnostic_scoring_tables is False
    assert KinaseScoringConfig().profile_missing_value_strategy == "strict"


def _dataset() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            "protein_id": ["MAPK14"],
        },
        index=phospho.index,
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
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )


def _kinase_result() -> KinaseWorkflowResult:
    dataset = _dataset()
    prediction_matrix = pd.DataFrame(
        {"MAP2K6": [0.9]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    score_matrix = pd.DataFrame(
        {"MAP2K6": [1.5]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_references(),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            combined_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def test_dataset_build_request_rejects_invalid_source_types_at_validator_boundary() -> (
    None
):
    request = DatasetBuildRequest(
        phospho=object(),
        site_metadata=object(),
    )
    with pytest.raises(
        UnsupportedInputFormatError,
        match="dataset build request phospho must be a pandas DataFrame or a file path",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_checks_organism_type_at_validator_boundary() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        organism="rat",
    )
    with pytest.raises(PhosPyInputError, match="organism must be an Organism"):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_has_default_preprocessing_config() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
    )
    validated = DatasetBuildRequestValidator().run(request)
    assert validated.preprocessing_config.missing_data.policy == "forbid"
    assert validated.preprocessing_config.missing_data.min_observed_values is None


def test_dataset_build_request_rejects_unknown_missing_data_policy() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(
                policy="unsupported"  # type: ignore[arg-type]
            )
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match="preprocessing_config.missing_data.policy must be one of",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_rejects_impute_policy_without_min_observed_values() -> (
    None
):
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(policy="impute_row_median")
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match=(
            "missing_data.min_observed_values must be an int when "
            "missing_data.policy='impute_row_median'"
        ),
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_rejects_min_observed_values_for_forbid_policy() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(
                policy="forbid",
                min_observed_values=1,
            ),
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match=(
            "missing_data.min_observed_values must be None when "
            "missing_data.policy='forbid'"
        ),
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_allows_ratio_total_protein_correction_policy() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        total=pd.DataFrame({"sample_a": [2.0]}, index=["MAPK14"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="ratio_to_total"
            )
        ),
    )
    validated = DatasetBuildRequestValidator().run(request)
    assert validated is request


def test_dataset_build_request_requires_total_for_ratio_total_protein_correction() -> (
    None
):
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="ratio_to_total"
            )
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match="policy='ratio_to_total' requires total input data",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_allows_site_matrix_build_from_metadata_policy() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
        ),
    )
    validated = DatasetBuildRequestValidator().run(request)
    assert validated is request


def test_dataset_build_request_allows_site_matrix_policy_overrides() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0], "sample_b": [2.0]}, index=["ROW_1"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["ROW_1"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            site_matrix=DatasetSiteMatrixConfig(
                policy="build_from_metadata",
                duplicate_site_strategy="aggregate_mean",
                missing_data_policy="require_min_observed_values",
                minimum_observed_values=1,
            )
        ),
    )

    validated = DatasetBuildRequestValidator().run(request)
    assert validated is request


def test_dataset_build_request_rejects_site_matrix_minimum_observed_without_policy() -> (
    None
):
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["ROW_1"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["ROW_1"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            site_matrix=DatasetSiteMatrixConfig(
                policy="build_from_metadata",
                missing_data_policy="drop_any_missing",
                minimum_observed_values=1,
            )
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="minimum_observed_values must be None unless",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_rejects_site_matrix_require_min_without_threshold() -> (
    None
):
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["ROW_1"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["ROW_1"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            site_matrix=DatasetSiteMatrixConfig(
                policy="build_from_metadata",
                missing_data_policy="require_min_observed_values",
            )
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="minimum_observed_values must be an int when",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_rejects_site_matrix_missing_policy_overrides_when_as_input() -> (
    None
):
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["ROW_1"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["ROW_1"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            site_matrix=DatasetSiteMatrixConfig(
                policy="as_input",
                missing_data_policy="retain_missing",
            )
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="missing_data_policy is only valid when site_matrix.policy='build_from_metadata'",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_rejects_site_matrix_duplicate_overrides_when_as_input() -> (
    None
):
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["ROW_1"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["ROW_1"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            site_matrix=DatasetSiteMatrixConfig(
                policy="as_input",
                duplicate_site_strategy="first",
            )
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="duplicate_site_strategy is only valid when site_matrix.policy='build_from_metadata'",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_allows_sample_metadata_pairs_comparison_policy() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0], "sample_b": [2.0]},
            index=["MAPK14;Y182;"],
        ),
        sample_metadata=pd.DataFrame(
            {"comparison_group": ["group_a", "group_b"]},
            index=["sample_a", "sample_b"],
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs")
        ),
    )
    validated = DatasetBuildRequestValidator().run(request)
    assert validated is request


def test_dataset_build_request_requires_sample_metadata_for_comparison_building() -> (
    None
):
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs")
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match="policy='sample_metadata_pairs' requires sample_metadata input data",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_rejects_duplicate_comparison_pairs_in_config() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0], "sample_b": [2.0]},
            index=["MAPK14;Y182;"],
        ),
        sample_metadata=pd.DataFrame(
            {"comparison_group": ["group_a", "group_b"]},
            index=["sample_a", "sample_b"],
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            comparisons=DatasetComparisonBuildingConfig(
                policy="sample_metadata_pairs",
                pairs=(("group_a", "group_b"), ("group_b", "group_a")),
            )
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match="contains duplicate pairs regardless of direction",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_rejects_pairs_when_comparison_policy_is_none() -> None:
    request = DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        preprocessing_config=DatasetPreprocessingConfig(
            comparisons=DatasetComparisonBuildingConfig(
                policy="none",
                pairs=(("group1", "group4"),),
            )
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match="comparisons.pairs must be None when comparisons.policy='none'",
    ):
        DatasetBuildRequestValidator().run(request)


def test_kinase_request_config_policy_fails_at_validator_boundary() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=1),
        prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=5),
        activity_config=None,
    )
    with pytest.raises(
        WorkflowValidationError,
        match="scoring_config.min_substrates must be greater than or equal to 2",
    ):
        KinaseWorkflowValidator().run(request)


def test_kinase_request_rejects_non_bool_diagnostic_scoring_policy() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables="yes",  # type: ignore[arg-type]
        ),
        prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=5),
        activity_config=None,
    )
    with pytest.raises(
        WorkflowValidationError,
        match="scoring_config.include_diagnostic_scoring_tables must be a bool",
    ):
        KinaseWorkflowValidator().run(request)


def test_kinase_request_rejects_unknown_profile_missing_value_strategy() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            profile_missing_value_strategy="unsupported",  # type: ignore[arg-type]
        ),
        prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=5),
        activity_config=None,
    )
    with pytest.raises(
        WorkflowValidationError,
        match="scoring_config.profile_missing_value_strategy must be one of",
    ):
        KinaseWorkflowValidator().run(request)


def test_kinase_request_rejects_unknown_prediction_mode() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            ensemble_size=5,
            mode="unsupported",  # type: ignore[arg-type]
        ),
        activity_config=None,
    )
    with pytest.raises(
        WorkflowValidationError,
        match="prediction_config.mode must be one of",
    ):
        KinaseWorkflowValidator().run(request)


def test_kinase_request_rejects_non_positive_adaptive_iterations() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            ensemble_size=5,
            mode="adaptive_ensemble",
            n_iterations=0,
        ),
        activity_config=None,
    )
    with pytest.raises(
        WorkflowValidationError,
        match="prediction_config.n_iterations must be greater than or equal to 1",
    ):
        KinaseWorkflowValidator().run(request)


def test_kinase_request_reference_compatibility_is_enforced_in_resolver() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.HUMAN,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=5),
        activity_config=None,
    )
    validated = KinaseWorkflowValidator().run(request)
    assert validated is request
    with pytest.raises(ReferenceCompatibilityError):
        KinaseWorkflow().run(request)


def test_kinase_activity_config_policy_fails_at_validator_boundary() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=5),
        activity_config=KinaseActivityConfig(
            enabled=True,
            threshold=0.6,
            min_substrates=0,
            top_n_substrates=20,
        ),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="activity_config.min_substrates must be greater than or equal to 1",
    ):
        KinaseWorkflowValidator().run(request)


def test_kinase_activity_top_n_config_policy_fails_at_validator_boundary() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(top_k=5, ensemble_size=5),
        activity_config=KinaseActivityConfig(
            enabled=True,
            threshold=0.6,
            min_substrates=1,
            top_n_substrates=0,
        ),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="activity_config.top_n_substrates must be greater than or equal to 1",
    ):
        KinaseWorkflowValidator().run(request)


def test_signalome_request_support_cutoff_policy_fails_at_validator_boundary() -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=SignalomeConfig(substrate_support_cutoff=1.5),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="signalome workflow request config.substrate_support_cutoff",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_request_network_threshold_policy_fails_at_validator_boundary() -> (
    None
):
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=SignalomeConfig(network_correlation_threshold=-0.1),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="signalome workflow request config.network_correlation_threshold",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_request_assignment_policy_fails_at_validator_boundary() -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=SignalomeConfig(assignment_policy="invalid"),  # type: ignore[arg-type]
    )
    with pytest.raises(
        WorkflowValidationError,
        match="signalome workflow request config.assignment_policy",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_request_network_policy_fails_at_validator_boundary() -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=SignalomeConfig(network_policy="invalid"),  # type: ignore[arg-type]
    )
    with pytest.raises(
        WorkflowValidationError,
        match="signalome workflow request config.network_policy",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_request_module_count_policy_fails_at_validator_boundary() -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=SignalomeConfig(module_count=0),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="signalome workflow request config.module_count",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_request_module_selection_threshold_policy_fails_at_boundary() -> (
    None
):
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=SignalomeConfig(module_selection_primary_correlation_threshold=1.2),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="module_selection_primary_correlation_threshold",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_request_module_selection_max_clusters_policy_fails_at_boundary() -> (
    None
):
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=SignalomeConfig(module_selection_max_clusters=0),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="module_selection_max_clusters",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_validator_requires_explicit_site_metadata_protein_id_column() -> (
    None
):
    kinase_result = _kinase_result()
    dataset_without_protein_id = AnalysisReadyPhosphoDataset(
        phospho=kinase_result.dataset.phospho,
        site_metadata=kinase_result.dataset.site_metadata.drop(columns=["protein_id"]),
        sample_metadata=kinase_result.dataset.sample_metadata,
        total=kinase_result.dataset.total,
        organism=kinase_result.dataset.organism,
        transformation_state=kinase_result.dataset.transformation_state,
    )
    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset_without_protein_id,
            references=kinase_result.references,
            scoring_result=kinase_result.scoring_result,
            prediction_result=kinase_result.prediction_result,
            activity_result=kinase_result.activity_result,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="site_metadata is missing required columns: protein_id",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_validator_rejects_empty_site_metadata_protein_id_values() -> None:
    kinase_result = _kinase_result()
    site_metadata = kinase_result.dataset.site_metadata.copy(deep=True)
    site_metadata.loc[:, "protein_id"] = [""]
    dataset_with_empty_protein_id = AnalysisReadyPhosphoDataset(
        phospho=kinase_result.dataset.phospho,
        site_metadata=site_metadata,
        sample_metadata=kinase_result.dataset.sample_metadata,
        total=kinase_result.dataset.total,
        organism=kinase_result.dataset.organism,
        transformation_state=kinase_result.dataset.transformation_state,
    )
    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset_with_empty_protein_id,
            references=kinase_result.references,
            scoring_result=kinase_result.scoring_result,
            prediction_result=kinase_result.prediction_result,
            activity_result=kinase_result.activity_result,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="site_metadata.protein_id must contain non-empty string values",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_validator_rejects_non_string_site_metadata_protein_id_values() -> (
    None
):
    kinase_result = _kinase_result()
    site_metadata = kinase_result.dataset.site_metadata.copy(deep=True)
    site_metadata = site_metadata.astype({"protein_id": object})
    site_metadata.loc[:, "protein_id"] = [123]  # type: ignore[list-item]
    dataset_with_numeric_protein_id = AnalysisReadyPhosphoDataset(
        phospho=kinase_result.dataset.phospho,
        site_metadata=site_metadata,
        sample_metadata=kinase_result.dataset.sample_metadata,
        total=kinase_result.dataset.total,
        organism=kinase_result.dataset.organism,
        transformation_state=kinase_result.dataset.transformation_state,
    )
    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset_with_numeric_protein_id,
            references=kinase_result.references,
            scoring_result=kinase_result.scoring_result,
            prediction_result=kinase_result.prediction_result,
            activity_result=kinase_result.activity_result,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="site_metadata.protein_id must contain non-empty string values",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_validator_does_not_cast_numeric_matrices(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    def _fail_astype(*args, **kwargs):
        raise AssertionError("validator must not coerce numeric matrices")

    monkeypatch.setattr(pd.DataFrame, "astype", _fail_astype)
    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request


def test_signalome_validator_allows_prediction_matrix_missingness() -> None:
    kinase_result = _kinase_result()
    prediction_with_missing = kinase_result.prediction_result.pred_mat.copy(deep=True)
    prediction_with_missing.iloc[0, 0] = float("nan")
    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=kinase_result.dataset,
            references=kinase_result.references,
            scoring_result=kinase_result.scoring_result,
            prediction_result=KinasePredictionResult(pred_mat=prediction_with_missing),
            activity_result=kinase_result.activity_result,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request


def test_signalome_validator_allows_downstream_score_matrix_missingness() -> None:
    kinase_result = _kinase_result()
    assert kinase_result.scoring_result.combined_scores is not None
    combined_with_missing = kinase_result.scoring_result.combined_scores.copy(deep=True)
    combined_with_missing.iloc[0, 0] = float("nan")

    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=kinase_result.dataset,
            references=kinase_result.references,
            scoring_result=KinaseScoringResult(
                profile_scores=kinase_result.scoring_result.profile_scores,
                combined_scores=combined_with_missing,
            ),
            prediction_result=kinase_result.prediction_result,
            activity_result=kinase_result.activity_result,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request


def test_signalome_validator_prefers_combined_scores_when_available() -> None:
    kinase_result = _kinase_result()
    assert kinase_result.scoring_result.combined_scores is not None
    invalid_combined = kinase_result.scoring_result.combined_scores.copy(deep=True)
    invalid_combined.iloc[0, 0] = float("inf")

    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=kinase_result.dataset,
            references=kinase_result.references,
            scoring_result=KinaseScoringResult(
                profile_scores=kinase_result.scoring_result.profile_scores,
                combined_scores=invalid_combined,
            ),
            prediction_result=kinase_result.prediction_result,
            activity_result=kinase_result.activity_result,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="kinase_result.scoring_result.combined_scores",
    ):
        SignalomeWorkflowValidator().run(request)
