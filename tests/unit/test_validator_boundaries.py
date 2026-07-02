from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
    DatasetTotalProteinCorrectionIdentityConfig,
    KinaseActivityConfig,
    KinaseAttritionPolicy,
    KinasePredictionConfig,
    KinaseScoringConfig,
    LocalisationRequirement,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
)
from phospy.api.requests import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.results import KinaseWorkflowResult
from phospy.errors import (
    PhosPyInputError,
    PhosPyValidationError,
    ReferenceCompatibilityError,
    UnsupportedInputFormatError,
    WorkflowValidationError,
)
from phospy.science.datasets.builders.validator import DatasetBuildRequestValidator
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.science.transformations.models import QuantitativeMeaning
from phospy.validation.datasets.preprocessing import (
    DatasetPreprocessingConfigValidator,
)
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)
from tests.support.unsafe_dataset_states import (
    unsafe_set_dataset_site_metadata_columns,
    unsafe_set_dataset_site_metadata_index,
)


def _site_keys(display_ids: list[str]) -> pd.Index:
    return site_key_index_from_display_ids(display_ids)


def _site_key(display_id: str) -> str:
    return str(_site_keys([display_id])[0])


def _dataset_state_kwargs(*, has_total_matrix: bool) -> dict[str, object]:
    return {
        "intensity_scale_state": supported_linear_intensity_scale_state(
            has_total_matrix=has_total_matrix
        ),
        "processing_state": supported_linear_processing_state(
            has_total_matrix=has_total_matrix
        ),
    }


def test_kinase_scoring_default_sets_two_substrate_support_floor() -> None:
    assert KinaseScoringConfig().min_substrates == 2
    assert KinaseScoringConfig().include_diagnostic_scoring_tables is False
    assert KinaseScoringConfig().include_substrate_contributions is False
    assert KinaseScoringConfig().profile_missing_value_strategy == "strict"
    assert KinaseScoringConfig().attrition_policy == KinaseAttritionPolicy()


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;"]
    site_index = _site_keys(display_ids)
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            "protein_id": ["MAPK14"],
        },
        index=site_index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        **_dataset_state_kwargs(has_total_matrix=False),
    )


def _rebuild_dataset_with_site_metadata(
    dataset: AnalysisReadyPhosphoDataset,
    site_metadata: pd.DataFrame,
) -> AnalysisReadyPhosphoDataset:
    return AnalysisReadyPhosphoDataset(
        phospho=dataset.phospho,
        site_metadata=site_metadata,
        sample_metadata=dataset.sample_metadata,
        total=dataset.total,
        organism=dataset.organism,
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=dataset.processing_state,
        preprocessing_report=dataset.preprocessing_report,
        protein_aware_preparation=dataset.protein_aware_preparation,
        provenance=dataset.provenance,
        allow_opaque_site_values=dataset.allow_opaque_site_values,
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


def _mixed_total_correction_dataset() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {"sample_a": [15.0, 7.0], "sample_b": [31.0, 15.0]},
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAATAAAAAAA"],
            "protein_id": ["MAPK14", "AKT1"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=phospho.index.copy(),
    )
    total = pd.DataFrame(
        {"sample_a": [3.0], "sample_b": [7.0]},
        index=pd.Index(["MAPK14"], name="protein_id"),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total",
                    identity=DatasetTotalProteinCorrectionIdentityConfig(
                        unmatched_policy=DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED
                    ),
                ),
            ),
        )
    )


def _kinase_result(
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> KinaseWorkflowResult:
    resolved_dataset = _dataset() if dataset is None else dataset
    site_index = resolved_dataset.phospho.index.copy()
    prediction_matrix = pd.DataFrame(
        {"MAP2K6": [0.9]},
        index=site_index,
    )
    score_matrix = pd.DataFrame(
        {"MAP2K6": [1.5]},
        index=site_index.copy(),
    )
    return KinaseWorkflowResult(
        dataset=resolved_dataset,
        references=_references(),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def _two_site_kinase_result() -> KinaseWorkflowResult:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    site_index = _site_keys(display_ids)
    dataset = AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAATAAAAAAA"],
                "protein_id": ["MAPK14", "AKT1"],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        **_dataset_state_kwargs(has_total_matrix=False),
    )
    prediction_matrix = pd.DataFrame(
        {"MAP2K6": [0.9, 0.1], "AKT1": [0.2, 0.8]},
        index=site_index.copy(),
    )
    score_matrix = pd.DataFrame(
        {"MAP2K6": [1.5, 0.5], "AKT1": [0.25, 1.25]},
        index=site_index.copy(),
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_references(),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def _mixed_total_kinase_request(
    *, allow_mixed_total_protein_quantitative_meaning: bool = False
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_mixed_total_correction_dataset(),
        references=_references(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            allow_mixed_total_protein_quantitative_meaning=allow_mixed_total_protein_quantitative_meaning,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )


def _mixed_total_signalome_request(
    *, allow_mixed_total_protein_quantitative_meaning: bool = False
) -> SignalomeWorkflowRequest:
    return SignalomeWorkflowRequest(
        kinase_result=_kinase_result(_mixed_total_correction_dataset()),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            allow_mixed_total_protein_quantitative_meaning=allow_mixed_total_protein_quantitative_meaning,
        ),
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
    assert validated.preprocessing_config.intensity_transform.policy == "identity"
    assert (
        validated.preprocessing_config.intensity_transform.pseudocount
        == pytest.approx(1.0)
    )
    assert validated.preprocessing_config.normalisation.policy == "none"
    assert validated.preprocessing_config.missing_data.policy == "forbid"
    assert validated.preprocessing_config.missing_data.min_observed_values is None


def test_dataset_build_request_requires_boolean_allow_opaque_site_values() -> None:
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
        allow_opaque_site_values="yes",  # type: ignore[arg-type]
    )
    with pytest.raises(
        PhosPyInputError,
        match="allow_opaque_site_values must be a bool",
    ):
        DatasetBuildRequestValidator().run(request)


@pytest.mark.parametrize(
    "value",
    [
        pytest.param("true", id="truthy-string"),
        pytest.param("false", id="falsey-string"),
        pytest.param(1, id="integer-one"),
        pytest.param(0, id="integer-zero"),
    ],
)
def test_dataset_build_request_requires_boolean_suspicious_declared_scale_override(
    value: object,
) -> None:
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
        allow_suspicious_declared_input_intensity_scale=value,  # type: ignore[arg-type]
    )
    with pytest.raises(
        PhosPyInputError,
        match=("allow_suspicious_declared_input_intensity_scale must be a bool"),
    ):
        DatasetBuildRequestValidator().run(request)


@pytest.mark.parametrize(
    "quantitative_meaning",
    [
        QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE.value,
        QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE.value,
    ],
)
def test_dataset_build_request_allows_supported_quantitative_meaning_literal(
    quantitative_meaning: str,
) -> None:
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
        quantitative_meaning=quantitative_meaning,
    )
    validated = DatasetBuildRequestValidator().run(request)
    assert validated is request


def test_dataset_build_request_rejects_unknown_quantitative_meaning_literal() -> None:
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
        quantitative_meaning="not_a_supported_meaning",
    )
    with pytest.raises(
        PhosPyInputError,
        match="dataset build request quantitative_meaning must be one of:",
    ):
        DatasetBuildRequestValidator().run(request)


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        pytest.param(
            lambda: DatasetIntensityTransformConfig(
                policy="unsupported",  # type: ignore[arg-type]
                pseudocount=1.0,
            ),
            "preprocessing_config.intensity_transform.policy must be one of",
            id="intensity-transform-policy-unsupported",
        ),
        pytest.param(
            lambda: DatasetNormalisationConfig(
                policy="unsupported",  # type: ignore[arg-type]
            ),
            "preprocessing_config.normalisation.policy must be one of",
            id="normalisation-policy-unsupported",
        ),
        pytest.param(
            lambda: DatasetMissingDataConfig(
                policy="unsupported"  # type: ignore[arg-type]
            ),
            "preprocessing_config.missing_data.policy must be one of",
            id="missing-data-policy-unsupported",
        ),
    ],
)
def test_dataset_build_request_rejects_unsupported_policy_literals(
    factory: object, pattern: str
) -> None:
    assert callable(factory)
    with pytest.raises(PhosPyInputError, match=pattern):
        factory()


@pytest.mark.parametrize(
    ("pseudocount", "pattern"),
    [
        pytest.param(
            -0.1,
            "intensity_transform.pseudocount must be greater than or equal to 0",
            id="below-min",
        ),
        pytest.param(0.0, None, id="at-min-valid"),
        pytest.param(1.0, None, id="inside-valid"),
        pytest.param(
            True,
            "intensity_transform.pseudocount must be a float or int",
            id="wrong-type",
        ),
        pytest.param(
            float("nan"), "intensity_transform.pseudocount must be finite", id="nan"
        ),
        pytest.param(
            float("inf"),
            "intensity_transform.pseudocount must be finite",
            id="infinite",
        ),
    ],
)
def test_dataset_build_request_log2_pseudocount_range_boundary(
    pseudocount: object, pattern: str | None
) -> None:
    if pattern is None:
        config = DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=pseudocount,  # type: ignore[arg-type]
        )
        assert config.pseudocount == pseudocount
        return

    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=pseudocount,  # type: ignore[arg-type]
        )


def test_dataset_build_request_rejects_impute_policy_without_min_observed_values() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match=(
            "missing_data.min_observed_values must be an int when "
            "missing_data.policy='impute_row_median'"
        ),
    ):
        DatasetMissingDataConfig(policy="impute_row_median")


def test_dataset_build_request_rejects_min_observed_values_for_forbid_policy() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "missing_data.min_observed_values must be None when "
            "missing_data.policy='forbid'"
        ),
    ):
        DatasetMissingDataConfig(
            policy="forbid",
            min_observed_values=1,
        )


def test_dataset_preprocessing_config_validator_allows_log2_subtract_log_total() -> (
    None
):
    config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        ),
        total_protein_correction=DatasetTotalProteinCorrectionConfig(
            policy="subtract_log_total"
        ),
    )
    validated = DatasetPreprocessingConfigValidator().run(config)
    assert validated is config


def test_dataset_preprocessing_config_validator_rejects_subtract_log_total_without_log2() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match=(
            "total_protein_correction.policy='subtract_log_total' requires "
            "log2-scale phospho and total values"
        ),
    ):
        DatasetPreprocessingConfigValidator().run(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="identity",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            )
        )


def test_dataset_build_request_allows_subtract_log_total_policy() -> None:
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
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total"
            ),
        ),
    )
    validated = DatasetBuildRequestValidator().run(request)
    assert validated is request


def test_dataset_build_request_requires_total_for_subtract_log_total() -> None:
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
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total"
            ),
        ),
    )
    with pytest.raises(
        PhosPyInputError,
        match="policy='subtract_log_total' requires total input data",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_build_request_rejects_removed_ratio_to_total_alias() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="preprocessing_config.total_protein_correction.policy must be one of",
    ):
        DatasetTotalProteinCorrectionConfig(policy="ratio_to_total")  # type: ignore[arg-type]


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
    with pytest.raises(
        PhosPyInputError,
        match=(
            "missing_data_policy='require_min_observed_values' is not supported "
            "for strict AnalysisReadyPhosphoDataset construction"
        ),
    ):
        DatasetSiteMatrixConfig(
            policy="build_from_metadata",
            duplicate_site_policy="aggregate_mean",
            missing_data_policy="require_min_observed_values",
            minimum_observed_values=1,
        )


def test_dataset_build_request_rejects_site_matrix_minimum_observed_without_policy() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match=(
            "minimum_observed_values is not supported for strict "
            "AnalysisReadyPhosphoDataset construction and must be None"
        ),
    ):
        DatasetSiteMatrixConfig(
            policy="build_from_metadata",
            missing_data_policy="drop_any_missing",
            minimum_observed_values=1,
        )


def test_dataset_build_request_rejects_site_matrix_require_min_without_threshold() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match=(
            "missing_data_policy='require_min_observed_values' is not supported "
            "for strict AnalysisReadyPhosphoDataset construction"
        ),
    ):
        DatasetSiteMatrixConfig(
            policy="build_from_metadata",
            missing_data_policy="require_min_observed_values",
        )


def test_dataset_build_request_rejects_site_matrix_missing_policy_overrides_when_as_input() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match=(
            "missing_data_policy='retain_missing' is not supported for strict "
            "AnalysisReadyPhosphoDataset construction"
        ),
    ):
        DatasetSiteMatrixConfig(
            policy="as_input",
            missing_data_policy="retain_missing",
        )


def test_dataset_build_request_rejects_site_matrix_duplicate_overrides_when_as_input() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match="duplicate_site_policy is only valid when site_matrix.policy='build_from_metadata'",
    ):
        DatasetSiteMatrixConfig(
            policy="as_input",
            duplicate_site_policy="first",
        )


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
    with pytest.raises(
        PhosPyInputError,
        match="contains duplicate pairs regardless of direction",
    ):
        DatasetComparisonBuildingConfig(
            policy="sample_metadata_pairs",
            pairs=(("group_a", "group_b"), ("group_b", "group_a")),
        )


def test_dataset_build_request_rejects_pairs_when_comparison_policy_is_none() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="comparisons.pairs must be None when comparisons.policy='none'",
    ):
        DatasetComparisonBuildingConfig(
            policy="none",
            pairs=(("group1", "group4"),),
        )


def test_kinase_request_config_policy_fails_at_validator_boundary() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="scoring_config.min_substrates must be greater than or equal to 2",
    ):
        KinaseScoringConfig(min_substrates=1)


def test_kinase_request_rejects_non_bool_diagnostic_scoring_policy() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="scoring_config.include_diagnostic_scoring_tables must be a bool",
    ):
        KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables="yes",  # type: ignore[arg-type]
        )


def test_kinase_request_rejects_non_bool_substrate_contribution_policy() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="scoring_config.include_substrate_contributions must be a bool",
    ):
        KinaseScoringConfig(
            min_substrates=2,
            include_substrate_contributions="yes",  # type: ignore[arg-type]
        )


def test_kinase_validator_rejects_unknown_reference_display_ambiguity_policy() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
        reference_display_ambiguity_policy="warn",  # type: ignore[arg-type]
    )

    with pytest.raises(
        WorkflowValidationError,
        match="reference_display_ambiguity_policy must be one of",
    ):
        KinaseWorkflowValidator().run(request)


def test_kinase_request_default_reference_display_ambiguity_policy_is_error() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
    )

    assert request.reference_display_ambiguity_policy == (
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR
    )


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        pytest.param(
            lambda: KinaseScoringConfig(
                min_substrates=2,
                profile_missing_value_strategy="unsupported",  # type: ignore[arg-type]
            ),
            "scoring_config.profile_missing_value_strategy must be one of",
            id="profile-missing-value-strategy-unsupported",
        ),
        pytest.param(
            lambda: KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=5,
                adaptive_ensemble_runs=5,
                mode="unsupported",  # type: ignore[arg-type]
            ),
            "prediction_config.mode must be one of",
            id="prediction-mode-unsupported",
        ),
    ],
)
def test_kinase_request_rejects_unsupported_literal_policies(
    factory: object, pattern: str
) -> None:
    assert callable(factory)
    with pytest.raises(WorkflowValidationError, match=pattern):
        factory()


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        pytest.param(
            lambda: KinasePredictionConfig(
                top_k=5,
                deterministic_max_selected_kinases=5,
                adaptive_ensemble_runs=5,
                mode="adaptive_ensemble",
                n_iterations=0,
                random_state=1,
            ),
            "prediction_config.n_iterations must be greater than or equal to 1",
            id="prediction-n-iterations-zero",
        ),
        pytest.param(
            lambda: KinaseActivityConfig(
                enabled=True,
                threshold=0.6,
                min_substrates=0,
                top_n_substrates=20,
            ),
            "activity_config.min_substrates must be greater than or equal to 1",
            id="activity-min-substrates-zero",
        ),
        pytest.param(
            lambda: KinaseActivityConfig(
                enabled=True,
                threshold=0.6,
                min_substrates=1,
                top_n_substrates=0,
            ),
            "activity_config.top_n_substrates must be greater than or equal to 1",
            id="activity-top-n-substrates-zero",
        ),
    ],
)
def test_positive_integer_policy_fails_at_validator_boundary(
    factory: object,
    pattern: str,
) -> None:
    # Consolidated matrix preserves field-specific boundary messages.
    assert callable(factory)
    with pytest.raises(WorkflowValidationError, match=pattern):
        factory()


@pytest.mark.parametrize(
    ("reference_source", "expected_exception"),
    [
        pytest.param(ReferencePreset.RAT, None, id="preset-rat-compatible"),
        pytest.param(
            ReferencePreset.HUMAN,
            ReferenceCompatibilityError,
            id="preset-human-incompatible",
        ),
        pytest.param(
            _references(),
            None,
            id="explicit-bundle-compatible",
        ),
        pytest.param(
            ReferenceBundle(
                organism=Organism.MOUSE,
                kinase_substrate_map=pd.DataFrame(
                    {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
                ),
                site_sequences=pd.DataFrame(
                    {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
                    index=pd.Index(["MAPK14;Y182;"], name="site_id"),
                ),
            ),
            ReferenceCompatibilityError,
            id="explicit-bundle-organism-mismatch",
        ),
    ],
)
def test_kinase_reference_compatibility_boundary_matrix(
    reference_source: object,
    expected_exception: type[Exception] | None,
) -> None:
    # Compatibility remains a domain boundary contract, not executor behavior.
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=reference_source,  # type: ignore[arg-type]
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    validated = KinaseWorkflowValidator().run(request)
    assert validated is request

    if expected_exception is None:
        try:
            KinaseWorkflow().run(request)
        except ReferenceCompatibilityError as exc:  # pragma: no cover - defensive
            raise AssertionError(
                "compatible references must not fail compatibility"
            ) from exc
        except Exception:
            # Other boundary failures (for example eligible kinase floor) are
            # outside the reference-compatibility contract in this matrix.
            pass
        return

    with pytest.raises(expected_exception):
        KinaseWorkflow().run(request)


def test_kinase_validator_rejects_mixed_total_protein_quantitative_meaning_by_default() -> (
    None
):
    request = _mixed_total_kinase_request()
    with pytest.raises(
        WorkflowValidationError,
        match=(
            "kinase workflow request dataset received a dataset with mixed "
            "total-protein quantitative meaning"
        ),
    ) as exc_info:
        KinaseWorkflowValidator().run(request)
    assert "uncorrected_rows=1" in str(exc_info.value)
    assert "unmatched_policy='allow_uncorrected'" in str(exc_info.value)


def test_kinase_validator_allows_mixed_total_protein_quantitative_meaning_with_opt_in() -> (
    None
):
    request = _mixed_total_kinase_request(
        allow_mixed_total_protein_quantitative_meaning=True
    )
    validated = KinaseWorkflowValidator().run(request)
    assert validated is request


def test_kinase_validator_allows_unknown_localisation_by_default() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    validated = KinaseWorkflowValidator().run(request)
    assert validated is request


def test_kinase_validator_can_require_localisation_probability_presence() -> None:
    request = KinaseWorkflowRequest(
        dataset=_dataset(),
        references=_references(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            localisation_requirement=LocalisationRequirement(require_present=True),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(request)
    message = str(exc_info.value)
    assert "kinase workflow request requires localisation metadata policy" in message
    assert (
        "missing required column=kinase workflow request "
        "dataset.site_metadata.localisation_confidence"
    ) in message


def test_kinase_validator_can_require_localisation_probability_threshold() -> None:
    source = _dataset()
    site_metadata = source.site_metadata
    site_metadata.loc[:, "localisation_probability"] = [0.6]
    dataset = _rebuild_dataset_with_site_metadata(source, site_metadata)
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=_references(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            localisation_requirement=LocalisationRequirement(minimum_probability=0.75),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )
    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(request)
    message = str(exc_info.value)
    assert "policy=require_threshold" in message
    assert "must be >= 0.750" in message


def test_kinase_validator_rejects_malformed_site_tokens() -> None:
    dataset = _dataset()
    unsafe_set_dataset_site_metadata_columns(dataset, {"site": ["FOO"]})
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=_references(),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="cannot resolve expected phosphosite residue",
    ):
        KinaseWorkflowValidator().run(request)


@pytest.mark.parametrize(
    ("sequence_value", "pattern"),
    [
        pytest.param(pd.NA, "missing or blank", id="missing"),
        pytest.param("AAAAAAAYAAAAAAAA", "must be odd length", id="even-length"),
        pytest.param(
            "AAAAAAASAAAAAAA",
            "centre residue must match the site token residue",
            id="centre-mismatch",
        ),
        pytest.param(
            "AAAAAA_SAAAAAAA",
            "centre residue must match the site token residue",
            id="underscore-gapped-mismatch",
        ),
    ],
)
def test_kinase_validator_requires_centred_sequence_context(
    sequence_value: object,
    pattern: str,
) -> None:
    dataset = _dataset()
    unsafe_set_dataset_site_metadata_columns(
        dataset,
        {"site_sequence": [sequence_value]},
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=_references(),
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )

    with pytest.raises(
        WorkflowValidationError,
        match=f"requires centred sequence context.*{pattern}",
    ) as exc_info:
        KinaseWorkflowValidator().run(request)
    assert _site_key("MAPK14;Y182;") in str(exc_info.value)


def test_kinase_validator_allows_gapped_flanks_when_centre_is_valid() -> None:
    source = _dataset()
    site_metadata = source.site_metadata
    site_metadata.loc[:, "site_sequence"] = ["________Y________"]
    dataset = _rebuild_dataset_with_site_metadata(source, site_metadata)
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=_references(),
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )

    validated = KinaseWorkflowValidator().run(request)
    assert validated is request


@pytest.mark.parametrize(
    (
        "config_path",
        "invalid_value",
        "factory",
        "expected_exception",
        "expected_message_fragment",
    ),
    [
        pytest.param(
            "config.scientific.substrate_support_cutoff",
            1.5,
            lambda value: SignalomeScientificConfig(substrate_support_cutoff=value),
            WorkflowValidationError,
            "signalome workflow request config.scientific.substrate_support_cutoff",
            id="substrate-support-cutoff-above-max",
        ),
        pytest.param(
            "config.output.network_correlation_threshold",
            -0.1,
            lambda value: SignalomeOutputConfig(network_correlation_threshold=value),
            WorkflowValidationError,
            "signalome workflow request config.output.network_correlation_threshold",
            id="network-correlation-threshold-below-min",
        ),
        pytest.param(
            "config.scientific.assignment_policy",
            "invalid",
            lambda value: SignalomeScientificConfig(assignment_policy=value),  # type: ignore[arg-type]
            WorkflowValidationError,
            "signalome workflow request config.scientific.assignment_policy",
            id="assignment-policy-invalid",
        ),
        pytest.param(
            "config.output.network_policy",
            "invalid",
            lambda value: SignalomeOutputConfig(network_policy=value),  # type: ignore[arg-type]
            WorkflowValidationError,
            "signalome workflow request config.output.network_policy",
            id="network-policy-invalid",
        ),
        pytest.param(
            "config.validation.score_preconditioning_policy",
            "invalid",
            lambda value: SignalomeValidationConfig(
                score_preconditioning_policy=value  # type: ignore[arg-type]
            ),
            WorkflowValidationError,
            "signalome workflow request config.validation.score_preconditioning_policy",
            id="score-preconditioning-policy-invalid",
        ),
        pytest.param(
            "config.clustering.module_count",
            True,
            lambda value: SignalomeClusteringConfig(module_count=value),  # type: ignore[arg-type]
            WorkflowValidationError,
            "signalome workflow request config.clustering.module_count must be an int",
            id="module-count-wrong-type-bool",
        ),
        pytest.param(
            "config.clustering.module_count",
            0,
            lambda value: SignalomeClusteringConfig(module_count=value),
            WorkflowValidationError,
            "signalome workflow request config.clustering.module_count",
            id="module-count-zero",
        ),
        pytest.param(
            "config.clustering.module_count",
            -1,
            lambda value: SignalomeClusteringConfig(module_count=value),
            WorkflowValidationError,
            "signalome workflow request config.clustering.module_count",
            id="module-count-negative",
        ),
        pytest.param(
            "config.clustering.module_selection_primary_correlation_threshold",
            1.2,
            lambda value: SignalomeClusteringConfig(
                module_selection_primary_correlation_threshold=value
            ),
            WorkflowValidationError,
            "module_selection_primary_correlation_threshold",
            id="module-selection-primary-threshold-above-max",
        ),
        pytest.param(
            "config.clustering.module_selection_fallback_correlation_threshold",
            -0.1,
            lambda value: SignalomeClusteringConfig(
                module_selection_fallback_correlation_threshold=value
            ),
            WorkflowValidationError,
            "module_selection_fallback_correlation_threshold",
            id="module-selection-fallback-threshold-below-min",
        ),
        pytest.param(
            "config.clustering.module_selection_max_clusters",
            0,
            lambda value: SignalomeClusteringConfig(
                module_selection_max_clusters=value
            ),
            WorkflowValidationError,
            "module_selection_max_clusters",
            id="module-selection-max-clusters-zero",
        ),
        pytest.param(
            "config.clustering.clustering_engine",
            "unsupported",
            lambda value: SignalomeClusteringConfig(
                clustering_engine=value  # type: ignore[arg-type]
            ),
            WorkflowValidationError,
            "signalome workflow request config.clustering.clustering_engine",
            id="clustering-engine-unsupported",
        ),
    ],
)
def test_signalome_config_boundary_invalid_case_matrix(
    config_path: str,
    invalid_value: object,
    factory: object,
    expected_exception: type[Exception],
    expected_message_fragment: str,
) -> None:
    # Domain boundary matrix keeps each public signalome field explicit.
    assert config_path.startswith("config.")
    assert callable(factory)
    with pytest.raises(expected_exception, match=expected_message_fragment):
        factory(invalid_value)  # type: ignore[misc]


def test_signalome_validator_rejects_mixed_total_protein_quantitative_meaning_by_default() -> (
    None
):
    request = _mixed_total_signalome_request()
    with pytest.raises(
        WorkflowValidationError,
        match=(
            "signalome workflow request kinase_result.dataset received a dataset "
            "with mixed total-protein quantitative meaning"
        ),
    ) as exc_info:
        SignalomeWorkflowValidator().run(request)
    assert "uncorrected_rows=1" in str(exc_info.value)
    assert "unmatched_policy='allow_uncorrected'" in str(exc_info.value)


def test_signalome_validator_allows_mixed_total_protein_quantitative_meaning_with_opt_in() -> (
    None
):
    request = _mixed_total_signalome_request(
        allow_mixed_total_protein_quantitative_meaning=True
    )
    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request


@pytest.mark.parametrize(
    ("module_count", "pattern"),
    [
        pytest.param(None, None, id="none-allowed"),
        pytest.param(1, None, id="valid-positive"),
    ],
)
def test_signalome_request_module_count_optional_positive_integer_boundary(
    module_count: object,
    pattern: str | None,
) -> None:
    if pattern is None:
        config = SignalomeClusteringConfig(module_count=module_count)  # type: ignore[arg-type]
        assert config.module_count == module_count
        return

    with pytest.raises(WorkflowValidationError, match=pattern):
        SignalomeClusteringConfig(module_count=module_count)  # type: ignore[arg-type]


def test_signalome_validator_accepts_network_minimum_paired_observation_setting() -> (
    None
):
    request = SignalomeWorkflowRequest(
        kinase_result=_two_site_kinase_result(),
        config=build_signalome_config(
            network_policy="absolute_threshold",
            network_min_paired_finite_observations=2,
            substrate_support_cutoff=0.5,
        ),
    )

    validated = SignalomeWorkflowValidator().run(request)

    assert validated is request


def test_signalome_validator_rejects_network_minimum_above_score_observations() -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_two_site_kinase_result(),
        config=build_signalome_config(
            network_min_paired_finite_observations=3,
            substrate_support_cutoff=0.5,
        ),
    )

    with pytest.raises(
        WorkflowValidationError,
        match=(
            "config.output.network_min_paired_finite_observations "
            "\\(3\\) cannot exceed available downstream score observations \\(2\\)"
        ),
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_request_max_exact_clustering_sites_policy_fails_at_boundary() -> (
    None
):
    with pytest.raises(
        TypeError,
        match="unexpected keyword argument 'max_exact_clustering_sites'",
    ):
        SignalomeConfig(max_exact_clustering_sites=0)  # type: ignore[call-arg]


def test_signalome_validator_accepts_complete_signalome_grouping_metadata() -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    validated = SignalomeWorkflowValidator().run(request)

    assert validated is request
    assert validated.kinase_result.dataset.site_metadata.loc[
        :, "protein_id"
    ].tolist() == ["MAPK14"]


def test_signalome_validator_requires_signalome_protein_grouping_metadata_column() -> (
    None
):
    kinase_result = _kinase_result()
    dataset_without_protein_id = AnalysisReadyPhosphoDataset(
        phospho=kinase_result.dataset.phospho,
        site_metadata=kinase_result.dataset.site_metadata.drop(columns=["protein_id"]),
        sample_metadata=kinase_result.dataset.sample_metadata,
        total=kinase_result.dataset.total,
        organism=kinase_result.dataset.organism,
        intensity_scale_state=kinase_result.dataset.intensity_scale_state,
        processing_state=kinase_result.dataset.processing_state,
    )
    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset_without_protein_id,
            references=kinase_result.references,
            scoring_result=kinase_result.scoring_result,
            prediction_result=kinase_result.prediction_result,
            activity_result=kinase_result.activity_result,
        ),
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert "Missing signalome protein grouping metadata: protein_id" in message
    assert "signalome protein grouping metadata requirement failed" in message
    assert "dataset.site_metadata.protein_id" in message
    assert "dataset-level protein-scoped row identity" in message
    assert "gene_symbol" in message
    assert "display_id" in message
    assert "identity requirement failed" not in message


def test_signalome_validator_reports_site_metadata_index_alignment_details() -> None:
    kinase_result = _kinase_result()
    mismatched_site_key = _site_key("AKT1;T308;")
    unsafe_set_dataset_site_metadata_index(
        kinase_result.dataset,
        pd.Index([mismatched_site_key], name="site_key"),
        update_site_key_column=True,
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert (
        "signalome workflow request kinase_result.dataset.site_metadata.index "
        "must exactly match signalome workflow request "
        "kinase_result.dataset.phospho.index"
    ) in message
    assert (
        "Only in signalome workflow request kinase_result.dataset.site_metadata.index: "
        f"'{mismatched_site_key}'"
    ) in message
    assert (
        "Only in signalome workflow request kinase_result.dataset.phospho.index: "
        f"'{_site_key('MAPK14;Y182;')}'"
    ) in message
    assert "First positional mismatch: position 0" in message


def _signalome_request_with_site_metadata(
    site_metadata: pd.DataFrame,
) -> SignalomeWorkflowRequest:
    kinase_result = _kinase_result()
    dataset = AnalysisReadyPhosphoDataset(
        phospho=kinase_result.dataset.phospho,
        site_metadata=site_metadata,
        sample_metadata=kinase_result.dataset.sample_metadata,
        total=kinase_result.dataset.total,
        organism=kinase_result.dataset.organism,
        intensity_scale_state=kinase_result.dataset.intensity_scale_state,
        processing_state=kinase_result.dataset.processing_state,
    )
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=kinase_result.references,
            scoring_result=kinase_result.scoring_result,
            prediction_result=kinase_result.prediction_result,
            activity_result=kinase_result.activity_result,
        ),
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )


def test_signalome_validator_rejects_blank_site_metadata_protein_id_values() -> None:
    kinase_result = _kinase_result()
    site_metadata = kinase_result.dataset.site_metadata.copy(deep=True)
    site_metadata.loc[:, "protein_id"] = [""]
    request = _signalome_request_with_site_metadata(site_metadata)

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert "Missing signalome protein grouping metadata: protein_id" in message
    assert "dataset.site_metadata.protein_id" in message
    assert "non-empty string values" in message
    assert "identity requirement failed" not in message


def test_signalome_validator_rejects_non_string_site_metadata_protein_id_values() -> (
    None
):
    kinase_result = _kinase_result()
    site_metadata = kinase_result.dataset.site_metadata.copy(deep=True)
    site_metadata = site_metadata.astype({"protein_id": object})
    site_metadata.loc[:, "protein_id"] = [123]  # type: ignore[list-item]
    request = _signalome_request_with_site_metadata(site_metadata)

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert "Missing signalome protein grouping metadata: protein_id" in message
    assert "dataset.site_metadata.protein_id" in message
    assert "non-empty string values" in message
    assert "identity requirement failed" not in message


def test_signalome_validator_can_require_localisation_probability_presence() -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            localisation_requirement=LocalisationRequirement(require_present=True),
        ),
    )
    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)
    message = str(exc_info.value)
    assert "signalome workflow request requires localisation metadata policy" in message
    assert (
        "missing required column=signalome workflow request "
        "kinase_result.dataset.site_metadata.localisation_confidence"
    ) in message
    assert "affected_rows=1" in message
    assert f"example_site_ids=['{_site_key('MAPK14;Y182;')}']" in message


def test_signalome_validator_can_require_localisation_probability_threshold() -> None:
    source = _dataset()
    site_metadata = source.site_metadata
    site_metadata.loc[:, "localisation_probability"] = [0.5]
    kinase_result = _kinase_result(
        _rebuild_dataset_with_site_metadata(source, site_metadata)
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            localisation_requirement=LocalisationRequirement(minimum_probability=0.75),
        ),
    )
    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)
    message = str(exc_info.value)
    assert "policy=require_threshold" in message
    assert "must be >= 0.750" in message
    assert f"example_site_ids=['{_site_key('MAPK14;Y182;')}']" in message


def test_signalome_validator_reports_invalid_localisation_probability_values() -> None:
    kinase_result = _kinase_result()
    unsafe_set_dataset_site_metadata_columns(
        kinase_result.dataset,
        {"localisation_probability": ["high"]},
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            localisation_requirement=LocalisationRequirement(require_present=True),
        ),
    )
    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)
    message = str(exc_info.value)
    assert (
        "invalid values in signalome workflow request "
        "kinase_result.dataset.site_metadata.localisation_probability"
    ) in message
    assert "affected_rows=1" in message
    assert f"example_site_ids=['{_site_key('MAPK14;Y182;')}']" in message


@pytest.mark.parametrize(
    ("sequence_value", "pattern"),
    [
        pytest.param(pd.NA, "missing or blank", id="missing"),
        pytest.param("AAAAAAAYAAAAAAAA", "must be odd length", id="even-length"),
        pytest.param(
            "AAAAAAASAAAAAAA",
            "centre residue must match the site token residue",
            id="centre-mismatch",
        ),
        pytest.param(
            "AAAAAA-SAAAAAAA",
            "centre residue must match the site token residue",
            id="gapped-mismatch",
        ),
    ],
)
def test_signalome_validator_requires_centred_sequence_context(
    sequence_value: object,
    pattern: str,
) -> None:
    kinase_result = _kinase_result()
    unsafe_set_dataset_site_metadata_columns(
        kinase_result.dataset,
        {"site_sequence": [sequence_value]},
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    with pytest.raises(
        WorkflowValidationError,
        match=f"requires centred sequence context.*{pattern}",
    ) as exc_info:
        SignalomeWorkflowValidator().run(request)
    message = str(exc_info.value)
    assert _site_key("MAPK14;Y182;") in message
    assert "protein grouping metadata requirement failed" not in message


def test_signalome_validator_allows_gapped_flanks_when_centre_is_valid() -> None:
    source = _dataset()
    site_metadata = source.site_metadata
    site_metadata.loc[:, "site_sequence"] = ["________Y________"]
    kinase_result = _kinase_result(
        _rebuild_dataset_with_site_metadata(source, site_metadata)
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request


def test_kinase_validator_does_not_filter_rows_for_localisation_policy() -> None:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    site_ids = _site_keys(display_ids)
    dataset = AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0]},
            index=site_ids.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_ids),
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAATAAAAAAA"],
                "protein_id": ["MAPK14", "AKT1"],
                "localisation_confidence": [0.2, pd.NA],
            },
            index=site_ids.copy(),
        ),
        organism=Organism.RAT,
        **_dataset_state_kwargs(has_total_matrix=False),
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["SEQ_A", "SEQ_R"]},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(
            localisation_requirement=LocalisationRequirement()
        ),
    )

    validated = KinaseWorkflowValidator().run(request)
    assert validated is request
    assert validated.dataset._borrow_phospho_frame().index.tolist() == site_ids.tolist()
    assert (
        validated.dataset._borrow_site_metadata_frame().index.tolist()
        == site_ids.tolist()
    )


def test_signalome_validator_does_not_filter_rows_for_localisation_policy() -> None:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    site_ids = _site_keys(display_ids)
    dataset = AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0]},
            index=site_ids.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_ids),
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAATAAAAAAA"],
                "protein_id": ["MAPK14", "AKT1"],
                "localisation_confidence": [0.2, pd.NA],
            },
            index=site_ids.copy(),
        ),
        organism=Organism.RAT,
        **_dataset_state_kwargs(has_total_matrix=False),
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["SEQ_A", "SEQ_R"]},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )
    signalome_scores = pd.DataFrame(
        {"MAP2K6": [1.0, 0.5]},
        index=site_ids.copy(),
    )
    signalome_predictions = pd.DataFrame(
        {"MAP2K6": [0.8, 0.2]},
        index=site_ids.copy(),
    )
    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=references,
            scoring_result=KinaseScoringResult(
                profile_scores=signalome_scores,
                rank_weighted_fusion_scores=signalome_scores,
            ),
            prediction_result=KinasePredictionResult(pred_mat=signalome_predictions),
            activity_result=None,
        ),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            localisation_requirement=LocalisationRequirement(),
        ),
    )

    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request
    assert (
        validated.kinase_result.dataset._borrow_phospho_frame().index.tolist()
        == site_ids.tolist()
    )
    assert (
        validated.kinase_result.dataset._borrow_site_metadata_frame().index.tolist()
        == site_ids.tolist()
    )


def test_signalome_validator_does_not_cast_numeric_matrices(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    def _fail_astype(*args, **kwargs):
        raise AssertionError("validator must not coerce numeric matrices")

    monkeypatch.setattr(pd.DataFrame, "astype", _fail_astype)
    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request


def test_signalome_validator_uses_internal_borrowed_dataframe_access(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(),
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    def _blocked_property(name: str) -> property:
        def _raiser(_self: object) -> object:
            raise AssertionError(f"validator must not use public property: {name}")

        return property(_raiser)

    monkeypatch.setattr(
        AnalysisReadyPhosphoDataset,
        "phospho",
        _blocked_property("dataset.phospho"),
    )
    monkeypatch.setattr(
        AnalysisReadyPhosphoDataset,
        "site_metadata",
        _blocked_property("dataset.site_metadata"),
    )
    monkeypatch.setattr(
        KinasePredictionResult,
        "pred_mat",
        _blocked_property("prediction_result.pred_mat"),
    )
    monkeypatch.setattr(
        KinaseScoringResult,
        "profile_scores",
        _blocked_property("scoring_result.profile_scores"),
    )
    monkeypatch.setattr(
        KinaseScoringResult,
        "rank_weighted_fusion_scores",
        _blocked_property("scoring_result.rank_weighted_fusion_scores"),
    )

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
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request


def test_signalome_validator_allows_downstream_score_matrix_missingness() -> None:
    kinase_result = _kinase_result()
    assert kinase_result.scoring_result.rank_weighted_fusion_scores is not None
    combined_with_missing = (
        kinase_result.scoring_result.rank_weighted_fusion_scores.copy(deep=True)
    )
    combined_with_missing.iloc[0, 0] = float("nan")

    request = SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=kinase_result.dataset,
            references=kinase_result.references,
            scoring_result=KinaseScoringResult(
                profile_scores=kinase_result.scoring_result.profile_scores,
                rank_weighted_fusion_scores=combined_with_missing,
            ),
            prediction_result=kinase_result.prediction_result,
            activity_result=kinase_result.activity_result,
        ),
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request


def test_signalome_validator_rejects_missing_site_metadata_protein_values() -> None:
    kinase_result = _kinase_result()
    site_metadata = kinase_result.dataset.site_metadata
    site_metadata.loc[:, "protein_id"] = np.nan
    request = _signalome_request_with_site_metadata(site_metadata)
    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert "Missing signalome protein grouping metadata: protein_id" in message
    assert "dataset.site_metadata.protein_id" in message
    assert "non-empty string values" in message
    assert "identity requirement failed" not in message


def test_signalome_validator_prefers_rank_weighted_fusion_scores_when_available() -> (
    None
):
    kinase_result = _kinase_result()
    assert kinase_result.scoring_result.rank_weighted_fusion_scores is not None
    invalid_combined = kinase_result.scoring_result.rank_weighted_fusion_scores.copy(
        deep=True
    )
    invalid_combined.iloc[0, 0] = float("inf")
    with pytest.raises(
        PhosPyValidationError,
        match="scoring_result.rank_weighted_fusion_scores must contain finite numeric values",
    ):
        SignalomeWorkflowRequest(
            kinase_result=KinaseWorkflowResult(
                dataset=kinase_result.dataset,
                references=kinase_result.references,
                scoring_result=KinaseScoringResult(
                    profile_scores=kinase_result.scoring_result.profile_scores,
                    rank_weighted_fusion_scores=invalid_combined,
                ),
                prediction_result=kinase_result.prediction_result,
                activity_result=kinase_result.activity_result,
            ),
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
