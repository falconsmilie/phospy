from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY,
    DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    DATASET_MISSING_DATA_POLICY_FORBID,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
    DATASET_NORMALISATION_POLICY_MEDIAN_CENTER,
    DATASET_NORMALISATION_POLICY_NONE,
    DATASET_SITE_MATRIX_POLICY_AS_INPUT,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
    KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
    KINASE_RELIABILITY_PROFILE_CUSTOM,
    KINASE_RELIABILITY_PROFILE_EXPLORATORY,
    KINASE_RELIABILITY_PROFILE_PRODUCTION,
    LOCALISATION_POLICY_REQUIRE_THRESHOLD,
    LOCALISATION_PRODUCTION_MINIMUM_PROBABILITY,
    REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT,
    SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT,
    SIGNALOME_MODE_EXPLORATORY_COMPATIBILITY,
    SIGNALOME_MODE_PRODUCTION,
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteSequenceResolutionConfig,
    DatasetTotalProteinCorrectionConfig,
    KinaseActivityConfig,
    KinaseAttritionPolicy,
    KinasePredictionConfig,
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    LocalisationRequirement,
    ProfileSelfInclusionPolicy,
    ReferenceContextCompatibilityPolicy,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
)
from phospy.errors import ContractValidationError, PhosPyInputError


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"policy": "unsupported"},
            "preprocessing_config.missing_data.policy must be one of",
        ),
        (
            {"policy": "forbid", "min_observed_values": 2},
            "missing_data.min_observed_values must be None when missing_data.policy='forbid'",
        ),
        (
            {"policy": "impute_row_median", "min_observed_values": None},
            "missing_data.min_observed_values must be an int",
        ),
        (
            {"policy": "impute_row_median", "min_observed_values": True},
            "missing_data.min_observed_values must be an int",
        ),
        (
            {"policy": "impute_row_median", "min_observed_values": 0},
            "missing_data.min_observed_values must be greater than or equal to 1",
        ),
    ],
)
def test_dataset_missing_data_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetMissingDataConfig(**kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"policy": "unsupported", "pseudocount": 1.0},
            "preprocessing_config.intensity_transform.policy must be one of",
        ),
    ],
)
def test_dataset_intensity_transform_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetIntensityTransformConfig(**kwargs)  # type: ignore[arg-type]


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
            float("inf"), "intensity_transform.pseudocount must be finite", id="inf"
        ),
    ],
)
def test_dataset_intensity_transform_pseudocount_range_matrix(
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


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        pytest.param(
            lambda: DatasetNormalisationConfig(
                policy="unsupported"  # type: ignore[arg-type]
            ),
            "preprocessing_config.normalisation.policy must be one of",
            id="normalisation-policy-unsupported",
        ),
        pytest.param(
            lambda: DatasetTotalProteinCorrectionConfig(
                policy="unsupported"  # type: ignore[arg-type]
            ),
            "preprocessing_config.total_protein_correction.policy must be one of",
            id="total-protein-correction-policy-unsupported",
        ),
    ],
)
def test_dataset_preprocessing_literal_policies_reject_unsupported_values(
    factory: object, pattern: str
) -> None:
    assert callable(factory)
    with pytest.raises(PhosPyInputError, match=pattern):
        factory()


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"policy": "invalid"},
            "preprocessing_config.site_matrix.policy must be one of",
        ),
        (
            {"policy": "build_from_metadata", "duplicate_site_policy": "invalid"},
            "site_matrix.duplicate_site_policy must be one of",
        ),
        (
            {"policy": "build_from_metadata", "missing_data_policy": "retain_missing"},
            "missing_data_policy='retain_missing' is not supported",
        ),
        (
            {
                "policy": "build_from_metadata",
                "missing_data_policy": "require_min_observed_values",
            },
            "missing_data_policy='require_min_observed_values' is not supported",
        ),
        (
            {"policy": "build_from_metadata", "minimum_observed_values": 1},
            "minimum_observed_values is not supported",
        ),
        (
            {"policy": "as_input", "duplicate_site_policy": "first"},
            "duplicate_site_policy is only valid when site_matrix.policy='build_from_metadata'",
        ),
    ],
)
def test_dataset_site_matrix_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetSiteMatrixConfig(**kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"policy": "invalid"},
            "preprocessing_config.comparisons.policy must be one of",
        ),
        (
            {"policy": "sample_metadata_pairs", "sample_group_column": ""},
            "comparisons.sample_group_column must be a non-empty string",
        ),
        (
            {"policy": "none", "pairs": (("a", "b"),)},
            "comparisons.pairs must be None when comparisons.policy='none'",
        ),
        (
            {"policy": "sample_metadata_pairs", "pairs": ()},
            "comparisons.pairs must contain at least one pair when provided",
        ),
        (
            {"policy": "sample_metadata_pairs", "pairs": ("a",)},
            "comparisons.pairs must contain only",
        ),
        (
            {"policy": "sample_metadata_pairs", "pairs": (("a", ""),)},
            "comparisons.pairs must contain non-empty right_group strings",
        ),
        (
            {"policy": "sample_metadata_pairs", "pairs": (("a", "a"),)},
            "comparisons.pairs cannot contain self-comparison pairs",
        ),
        (
            {
                "policy": "sample_metadata_pairs",
                "pairs": (("a", "b"), ("b", "a")),
            },
            "contains duplicate pairs regardless of direction",
        ),
    ],
)
def test_dataset_comparison_building_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetComparisonBuildingConfig(**kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"mode": "invalid"},
            "site_sequence_resolution.mode must be one of",
        ),
        (
            {"conflict_policy": "invalid"},
            "site_sequence_resolution.conflict_policy must be one of",
        ),
        (
            {"fasta_path": "https://example.org/test.fasta"},
            "fasta_path must be a local filesystem path",
        ),
        (
            {"fasta_path": "local.fasta", "flank_size": 0},
            "flank_size must be greater than or equal to 1",
        ),
        (
            {"accession_column": ""},
            "site_sequence_resolution.accession_column must be a non-empty string",
        ),
        (
            {"site_column": ""},
            "site_sequence_resolution.site_column must be a non-empty string",
        ),
    ],
)
def test_dataset_site_sequence_resolution_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetSiteSequenceResolutionConfig(**kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"missing_data": object()},
            "preprocessing_config.missing_data must be a DatasetMissingDataConfig",
        ),
        (
            {"intensity_transform": object()},
            "preprocessing_config.intensity_transform must be a DatasetIntensityTransformConfig",
        ),
        (
            {"comparisons": object()},
            "preprocessing_config.comparisons must be a DatasetComparisonBuildingConfig",
        ),
        (
            {"site_sequence_resolution": object()},
            "preprocessing_config.site_sequence_resolution must be a DatasetSiteSequenceResolutionConfig",
        ),
    ],
)
def test_dataset_preprocessing_config_rejects_wrong_nested_types(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetPreprocessingConfig(**kwargs)  # type: ignore[arg-type]


def test_dataset_preprocessing_config_presets_return_expected_values() -> None:
    default = DatasetPreprocessingConfig.default()
    strict = DatasetPreprocessingConfig.strict()
    raw = DatasetPreprocessingConfig.from_raw_phosphosite_table()

    assert (
        default.intensity_transform.policy
        == DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY
    )
    assert default.normalisation.policy == DATASET_NORMALISATION_POLICY_NONE
    assert default.missing_data.policy == DATASET_MISSING_DATA_POLICY_FORBID
    assert (
        default.total_protein_correction.policy
        == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )
    assert default.site_matrix.policy == DATASET_SITE_MATRIX_POLICY_AS_INPUT
    assert default.comparisons.policy == DATASET_COMPARISON_BUILDING_POLICY_NONE

    assert (
        strict.intensity_transform.policy == DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY
    )
    assert strict.normalisation.policy == DATASET_NORMALISATION_POLICY_NONE
    assert strict.missing_data.policy == DATASET_MISSING_DATA_POLICY_FORBID
    assert (
        strict.total_protein_correction.policy
        == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )
    assert strict.site_matrix.policy == DATASET_SITE_MATRIX_POLICY_AS_INPUT
    assert strict.comparisons.policy == DATASET_COMPARISON_BUILDING_POLICY_NONE

    assert raw.intensity_transform.policy == DATASET_INTENSITY_TRANSFORM_POLICY_LOG2
    assert raw.normalisation.policy == DATASET_NORMALISATION_POLICY_MEDIAN_CENTER
    assert raw.missing_data.policy == DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN
    assert raw.missing_data.min_observed_values == 1
    assert (
        raw.total_protein_correction.policy
        == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )
    assert raw.site_matrix.policy == DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA
    assert raw.comparisons.policy == DATASET_COMPARISON_BUILDING_POLICY_NONE


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"min_substrates": 1},
            "scoring_config.min_substrates must be greater than or equal to 2",
        ),
        (
            {"include_diagnostic_scoring_tables": "yes"},
            "scoring_config.include_diagnostic_scoring_tables must be a bool",
        ),
        (
            {"include_substrate_contributions": "yes"},
            "scoring_config.include_substrate_contributions must be a bool",
        ),
        (
            {"profile_missing_value_strategy": "invalid"},
            "scoring_config.profile_missing_value_strategy must be one of",
        ),
        (
            {"profile_self_inclusion_policy": "invalid"},
            "scoring_config.profile_self_inclusion_policy must be one of",
        ),
        (
            {"reference_context_compatibility_policy": "invalid"},
            "scoring_config.reference_context_compatibility_policy must be one of",
        ),
        (
            {"allow_mixed_total_protein_quantitative_meaning": "yes"},
            "scoring_config.allow_mixed_total_protein_quantitative_meaning must be a bool",
        ),
        (
            {"attrition_policy": object()},
            "scoring_config.attrition_policy must be KinaseAttritionPolicy",
        ),
    ],
)
def test_kinase_scoring_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(ContractValidationError, match=pattern):
        KinaseScoringConfig(reliability_profile="custom", **kwargs)  # type: ignore[arg-type]


def test_kinase_attrition_policy_accepts_valid_thresholds() -> None:
    policy = KinaseAttritionPolicy(
        minimum_reference_overlap_fraction=0.25,
        minimum_sequence_supported_fraction=0.5,
        minimum_scored_fraction=1.0,
        on_violation="error",
    )

    assert policy.minimum_reference_overlap_fraction == pytest.approx(0.25)
    assert policy.minimum_sequence_supported_fraction == pytest.approx(0.5)
    assert policy.minimum_scored_fraction == pytest.approx(1.0)
    assert policy.on_violation == "error"


def test_kinase_attrition_policy_rejects_negative_fraction() -> None:
    with pytest.raises(
        ContractValidationError,
        match="attrition_policy.minimum_reference_overlap_fraction",
    ):
        KinaseAttritionPolicy(minimum_reference_overlap_fraction=-0.1)


def test_kinase_attrition_policy_rejects_fraction_above_one() -> None:
    with pytest.raises(
        ContractValidationError,
        match="attrition_policy.minimum_sequence_supported_fraction",
    ):
        KinaseAttritionPolicy(minimum_sequence_supported_fraction=1.1)


@pytest.mark.parametrize("value", [float("nan"), float("inf")])
def test_kinase_attrition_policy_rejects_nonfinite_fraction(value: float) -> None:
    with pytest.raises(
        ContractValidationError,
        match="attrition_policy.minimum_scored_fraction must be finite",
    ):
        KinaseAttritionPolicy(minimum_scored_fraction=value)


def test_kinase_attrition_policy_rejects_invalid_violation_mode() -> None:
    with pytest.raises(
        ContractValidationError,
        match="attrition_policy.on_violation must be one of",
    ):
        KinaseAttritionPolicy(on_violation="ignore")  # type: ignore[arg-type]


def test_kinase_config_uses_default_attrition_policy() -> None:
    config = KinaseScoringConfig.exploratory()

    assert config.attrition_policy == KinaseAttritionPolicy(
        minimum_reference_overlap_fraction=0.0,
        minimum_sequence_supported_fraction=0.0,
        minimum_scored_fraction=0.0,
        on_violation="warn",
    )
    assert config.profile_self_inclusion_policy is ProfileSelfInclusionPolicy.ALLOW
    assert (
        config.reference_context_compatibility_policy
        is ReferenceContextCompatibilityPolicy.REQUIRE_KNOWN_MATCH
    )
    assert config.reliability_profile is KINASE_RELIABILITY_PROFILE_EXPLORATORY
    assert (
        config.effective_reliability_profile is KINASE_RELIABILITY_PROFILE_EXPLORATORY
    )
    assert (
        config.requested_reliability_profile is KINASE_RELIABILITY_PROFILE_EXPLORATORY
    )


def test_kinase_scoring_config_requires_explicit_reliability_profile() -> None:
    with pytest.raises(
        ContractValidationError, match="reliability_profile is required"
    ):
        KinaseScoringConfig()  # type: ignore[call-arg]


def test_kinase_scoring_config_accepts_attrition_policy() -> None:
    policy = KinaseAttritionPolicy(
        minimum_reference_overlap_fraction=0.2,
        minimum_sequence_supported_fraction=0.3,
        minimum_scored_fraction=0.4,
        on_violation="error",
    )

    config = KinaseScoringConfig(reliability_profile="custom", attrition_policy=policy)

    assert config.attrition_policy is policy


def test_localisation_requirement_production_requires_probability_threshold() -> None:
    requirement = LocalisationRequirement.production_site_level()

    assert requirement.policy == LOCALISATION_POLICY_REQUIRE_THRESHOLD
    assert requirement.require_present is True
    assert requirement.requires_probability_column is True
    assert requirement.minimum_probability == pytest.approx(
        LOCALISATION_PRODUCTION_MINIMUM_PROBABILITY
    )
    assert LocalisationRequirement().requires_probability_column is False
    assert LocalisationRequirement().minimum_probability is None


def test_kinase_scoring_exploratory_matches_historical_default() -> None:
    direct = KinaseScoringConfig.exploratory()
    exploratory = KinaseScoringConfig.exploratory()

    assert exploratory == direct
    assert exploratory.reliability_profile is KINASE_RELIABILITY_PROFILE_EXPLORATORY
    assert (
        exploratory.requested_reliability_profile
        is KINASE_RELIABILITY_PROFILE_EXPLORATORY
    )


def test_kinase_scoring_default_is_deprecated_exploratory_alias() -> None:
    with pytest.warns(DeprecationWarning, match="exploratory"):
        default = KinaseScoringConfig.default()

    assert default == KinaseScoringConfig.exploratory()
    assert default.reliability_profile is KINASE_RELIABILITY_PROFILE_EXPLORATORY


def test_kinase_production_requires_explicit_attrition_thresholds() -> None:
    with pytest.raises(TypeError):
        KinaseScoringConfig.production()  # type: ignore[call-arg]


def test_kinase_production_config_uses_strict_reliability_invariants() -> None:
    production = KinaseScoringConfig.production(
        minimum_reference_overlap_fraction=0.25,
        minimum_sequence_supported_fraction=0.5,
        minimum_scored_fraction=0.75,
    )

    assert (
        production.localisation_requirement
        == LocalisationRequirement.production_site_level()
    )
    assert production.min_substrates == 5
    assert (
        production.profile_self_inclusion_policy
        is ProfileSelfInclusionPolicy.LEAVE_ONE_OUT
    )
    assert production.attrition_policy == KinaseAttritionPolicy(
        minimum_reference_overlap_fraction=0.25,
        minimum_sequence_supported_fraction=0.5,
        minimum_scored_fraction=0.75,
        on_violation="error",
    )
    assert production.reliability_profile is KINASE_RELIABILITY_PROFILE_PRODUCTION
    assert (
        production.requested_reliability_profile
        is KINASE_RELIABILITY_PROFILE_PRODUCTION
    )
    assert (
        production.profile_missing_value_strategy
        == KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    )


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"min_substrates": 4},
            "scoring_config.min_substrates must be at least 5",
        ),
        (
            {"profile_self_inclusion_policy": "allow"},
            "profile_self_inclusion_policy must be leave_one_out",
        ),
        (
            {"localisation_requirement": LocalisationRequirement()},
            "must reject unknown localisation",
        ),
        (
            {
                "localisation_requirement": LocalisationRequirement(
                    require_present=True,
                    minimum_probability=0.5,
                )
            },
            "minimum_probability must be at least",
        ),
        (
            {"attrition_policy": KinaseAttritionPolicy(on_violation="warn")},
            "attrition_policy.on_violation must be error",
        ),
        (
            {
                "attrition_policy": KinaseAttritionPolicy(
                    minimum_reference_overlap_fraction=0.25,
                    minimum_sequence_supported_fraction=0.0,
                    minimum_scored_fraction=0.75,
                    on_violation="error",
                )
            },
            "minimum_sequence_supported_fraction must be set above 0.0",
        ),
    ],
)
def test_kinase_production_label_rejects_weakened_invariants(
    kwargs: dict[str, object],
    pattern: str,
) -> None:
    valid_values: dict[str, object] = {
        "min_substrates": 5,
        "profile_self_inclusion_policy": "leave_one_out",
        "localisation_requirement": LocalisationRequirement.production_site_level(),
        "attrition_policy": KinaseAttritionPolicy(
            minimum_reference_overlap_fraction=0.25,
            minimum_sequence_supported_fraction=0.5,
            minimum_scored_fraction=0.75,
            on_violation="error",
        ),
    }
    valid_values.update(kwargs)

    with pytest.raises(ContractValidationError, match=pattern):
        KinaseScoringConfig(
            **valid_values,
            reliability_profile=KINASE_RELIABILITY_PROFILE_PRODUCTION,
        )


def test_kinase_custom_profile_accepts_modified_exploratory_values() -> None:
    modified = KinaseScoringConfig(
        min_substrates=3,
        reliability_profile=KinaseReliabilityProfile.CUSTOM,
    )

    assert modified.reliability_profile is KINASE_RELIABILITY_PROFILE_CUSTOM
    assert modified.requested_reliability_profile is KINASE_RELIABILITY_PROFILE_CUSTOM


def test_kinase_explicit_exploratory_rejects_modified_values() -> None:
    with pytest.raises(ContractValidationError, match="requires the exploratory"):
        KinaseScoringConfig(
            min_substrates=3,
            reliability_profile=KinaseReliabilityProfile.EXPLORATORY,
        )


def test_kinase_explicit_custom_accepts_modified_values() -> None:
    config = KinaseScoringConfig(
        min_substrates=3,
        reliability_profile=KinaseReliabilityProfile.CUSTOM,
    )

    assert config.reliability_profile is KINASE_RELIABILITY_PROFILE_CUSTOM
    assert config.requested_reliability_profile is KINASE_RELIABILITY_PROFILE_CUSTOM


def test_kinase_scoring_config_presets_return_expected_values() -> None:
    default = KinaseScoringConfig.exploratory()
    strict_missing = KinaseScoringConfig.strict_missing_values()

    assert (
        default.profile_missing_value_strategy
        == KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    )
    assert default.include_substrate_contributions is False
    assert default.attrition_policy == KinaseAttritionPolicy()
    assert default.profile_self_inclusion_policy is ProfileSelfInclusionPolicy.ALLOW
    assert (
        strict_missing.profile_missing_value_strategy
        == KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    )
    assert strict_missing.include_substrate_contributions is False
    assert strict_missing.attrition_policy == KinaseAttritionPolicy()
    assert (
        strict_missing.profile_self_inclusion_policy is ProfileSelfInclusionPolicy.ALLOW
    )
    assert (
        default.reference_context_compatibility_policy
        == REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH
    )


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"mode": "invalid"},
            "prediction_config.mode must be one of",
        ),
        (
            {"adaptive_policy": "invalid"},
            "prediction_config.adaptive_policy must be one of",
        ),
        (
            {"top_k": 0},
            "prediction_config.top_k must be greater than or equal to 1",
        ),
        (
            {"deterministic_max_selected_kinases": 0},
            "prediction_config.deterministic_max_selected_kinases must be greater than or equal to 1",
        ),
        (
            {"adaptive_ensemble_runs": 0},
            "prediction_config.adaptive_ensemble_runs must be greater than or equal to 1",
        ),
        (
            {"n_iterations": 0},
            "prediction_config.n_iterations must be greater than or equal to 1",
        ),
    ],
)
def test_kinase_prediction_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(ContractValidationError, match=pattern):
        KinasePredictionConfig(**kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("factory", "value", "attribute_name", "pattern"),
    [
        pytest.param(
            KinasePredictionConfig,
            None,
            "random_state",
            None,
            id="kinase-random-state-none",
        ),
        pytest.param(
            KinasePredictionConfig,
            True,
            "random_state",
            "prediction_config.random_state must be an int",
            id="kinase-random-state-wrong-type",
        ),
        pytest.param(
            KinasePredictionConfig,
            -1,
            "random_state",
            "prediction_config.random_state must be greater than or equal to 0",
            id="kinase-random-state-negative",
        ),
        pytest.param(
            KinasePredictionConfig,
            1,
            "random_state",
            None,
            id="kinase-random-state-valid-positive",
        ),
        pytest.param(
            SignalomeClusteringConfig,
            None,
            "module_count",
            None,
            id="signalome-module-count-none",
        ),
        pytest.param(
            SignalomeClusteringConfig,
            True,
            "module_count",
            "signalome workflow request config.clustering.module_count must be an int",
            id="signalome-module-count-wrong-type",
        ),
        pytest.param(
            SignalomeClusteringConfig,
            0,
            "module_count",
            "signalome workflow request config.clustering.module_count",
            id="signalome-module-count-zero",
        ),
        pytest.param(
            SignalomeClusteringConfig,
            -1,
            "module_count",
            "signalome workflow request config.clustering.module_count",
            id="signalome-module-count-negative",
        ),
        pytest.param(
            SignalomeClusteringConfig,
            1,
            "module_count",
            None,
            id="signalome-module-count-valid-positive",
        ),
    ],
)
def test_optional_positive_integer_config_fields_self_validate(
    factory: type[KinasePredictionConfig] | type[SignalomeClusteringConfig],
    value: object,
    attribute_name: str,
    pattern: str | None,
) -> None:
    # Consolidated matrix for optional integer boundaries keeps field-level messages explicit.
    kwargs = {attribute_name: value}
    if pattern is None:
        config = factory(**kwargs)  # type: ignore[arg-type]
        assert getattr(config, attribute_name) == value
        return

    with pytest.raises(ContractValidationError, match=pattern):
        factory(**kwargs)  # type: ignore[arg-type]


def test_kinase_prediction_config_has_explicit_mode_specific_sizes() -> None:
    config = KinasePredictionConfig()

    assert config.deterministic_max_selected_kinases >= 1
    assert config.adaptive_ensemble_runs >= 1


def test_kinase_prediction_config_presets_return_expected_values() -> None:
    deterministic = KinasePredictionConfig.deterministic()
    adaptive = KinasePredictionConfig.adaptive_reproducible(random_state=1)

    assert deterministic.mode == KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING
    assert deterministic.random_state is None
    assert adaptive.mode == KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE
    assert adaptive.random_state == 1


def test_kinase_prediction_adaptive_reproducible_rejects_invalid_seed() -> None:
    with pytest.raises(
        ContractValidationError,
        match="prediction_config.random_state must be greater than or equal to 0",
    ):
        KinasePredictionConfig.adaptive_reproducible(random_state=-1)


def test_config_presets_return_frozen_dataclass_objects() -> None:
    config = KinasePredictionConfig.deterministic()
    with pytest.raises(FrozenInstanceError):
        config.mode = KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE  # type: ignore[misc]


def test_kinase_prediction_config_rejects_removed_ensemble_size_alias() -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument 'ensemble_size'"):
        KinasePredictionConfig(ensemble_size=25)  # type: ignore[call-arg]


@pytest.mark.parametrize(
    ("kwargs", "pattern"),
    [
        (
            {"enabled": "yes"},
            "activity_config.enabled must be a bool",
        ),
        (
            {"method": "invalid"},
            "activity_config.method must be one of",
        ),
        (
            {"min_substrates": 0},
            "activity_config.min_substrates must be greater than or equal to 1",
        ),
        (
            {"top_n_substrates": 0},
            "activity_config.top_n_substrates must be greater than or equal to 1",
        ),
        (
            {"ksea_min_substrates": 0},
            "activity_config.ksea_min_substrates must be greater than or equal to 1",
        ),
        (
            {"ksea_p_value_method": "invalid"},
            "activity_config.ksea_p_value_method must be one of",
        ),
        (
            {"ksea_adjust_p_values": "yes"},
            "activity_config.ksea_adjust_p_values must be a bool",
        ),
        (
            {"ssgsea_min_substrates": 0},
            "activity_config.ssgsea_min_substrates must be greater than or equal to 1",
        ),
        (
            {"ssgsea_ranking_direction": "invalid"},
            "activity_config.ssgsea_ranking_direction must be one of",
        ),
        (
            {"ssgsea_permutations": -1},
            "activity_config.ssgsea_permutations must be greater than or equal to 0",
        ),
        (
            {"ssgsea_random_seed": -1},
            "activity_config.ssgsea_random_seed must be greater than or equal to 0",
        ),
        (
            {"ssgsea_permutations": 1, "ssgsea_random_seed": None},
            "activity_config.ssgsea_random_seed must be set",
        ),
        (
            {"ssgsea_adjust_p_values": "yes"},
            "activity_config.ssgsea_adjust_p_values must be a bool",
        ),
    ],
)
def test_kinase_activity_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(ContractValidationError, match=pattern):
        KinaseActivityConfig(**kwargs)  # type: ignore[arg-type]


def test_kinase_activity_ssgsea_permutation_helper_requires_seed_and_is_reproducible() -> (
    None
):
    with pytest.raises(TypeError):
        KinaseActivityConfig.ssgsea_with_permutation_significance(  # type: ignore[call-arg]
            permutations=10
        )
    with pytest.raises(ContractValidationError, match="must be greater than 0"):
        KinaseActivityConfig.ssgsea_with_permutation_significance(
            permutations=0,
            random_seed=3,
        )

    first = KinaseActivityConfig.ssgsea_with_permutation_significance(
        permutations=10,
        random_seed=3,
    )
    second = KinaseActivityConfig.ssgsea_with_permutation_significance(
        permutations=10,
        random_seed=3,
    )

    assert first == second
    assert first.method == KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT
    assert first.ssgsea_permutations == 10
    assert first.ssgsea_random_seed == 3


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        (
            lambda: SignalomeOutputConfig(network_policy="invalid"),  # type: ignore[arg-type]
            "signalome workflow request config.output.network_policy",
        ),
        (
            lambda: SignalomeOutputConfig(network_min_paired_finite_observations=2),
            (
                "signalome workflow request config.output."
                "network_min_paired_finite_observations"
            ),
        ),
        (
            lambda: SignalomeOutputConfig(network_min_paired_finite_observations=1),
            (
                "signalome workflow request config.output."
                "network_min_paired_finite_observations"
            ),
        ),
        (
            lambda: SignalomeOutputConfig(
                network_min_paired_finite_observations=True  # type: ignore[arg-type]
            ),
            (
                "signalome workflow request config.output."
                "network_min_paired_finite_observations must be an int"
            ),
        ),
        (
            lambda: SignalomeOutputConfig(network_min_paired_finite_observations=0),
            (
                "signalome workflow request config.output."
                "network_min_paired_finite_observations"
            ),
        ),
        (
            lambda: SignalomeOutputConfig(network_min_paired_finite_observations=-1),
            (
                "signalome workflow request config.output."
                "network_min_paired_finite_observations"
            ),
        ),
        (
            lambda: SignalomeScientificConfig(assignment_policy="invalid"),  # type: ignore[arg-type]
            "signalome workflow request config.scientific.assignment_policy",
        ),
        (
            lambda: SignalomeValidationConfig(
                score_preconditioning_policy="invalid"  # type: ignore[arg-type]
            ),
            "signalome workflow request config.validation.score_preconditioning_policy",
        ),
        (
            lambda: SignalomeValidationConfig(
                allow_mixed_total_protein_quantitative_meaning="yes"  # type: ignore[arg-type]
            ),
            "signalome workflow request config.validation.allow_mixed_total_protein_quantitative_meaning",
        ),
        (
            lambda: SignalomeValidationConfig(
                reference_context_compatibility_policy="invalid"  # type: ignore[arg-type]
            ),
            (
                "signalome workflow request config.validation."
                "reference_context_compatibility_policy"
            ),
        ),
        (
            lambda: SignalomeClusteringConfig(module_selection_max_clusters=0),
            "module_selection_max_clusters",
        ),
        (
            lambda: SignalomeClusteringConfig(clustering_engine="invalid"),  # type: ignore[arg-type]
            "signalome workflow request config.clustering.clustering_engine",
        ),
        (
            lambda: SignalomeClusteringConfig(candidate_scoring_policy="invalid"),  # type: ignore[arg-type]
            "signalome workflow request config.clustering.candidate_scoring_policy",
        ),
        (
            lambda: SignalomePerformanceConfig(max_exact_tree_sites=0),
            "signalome workflow request config.performance.max_exact_tree_sites",
        ),
        (
            lambda: SignalomePerformanceConfig(max_full_candidate_scoring_sites=0),
            (
                "signalome workflow request config.performance."
                "max_full_candidate_scoring_sites"
            ),
        ),
        (
            lambda: SignalomeConfig(mode="invalid"),  # type: ignore[arg-type]
            "signalome workflow request config.mode",
        ),
        (
            lambda: SignalomeConfig(
                validation=SignalomeValidationConfig(
                    localisation_requirement=LocalisationRequirement()
                )
            ),
            "config.mode='production' requires",
        ),
        (
            lambda: SignalomeConfig(
                output=SignalomeOutputConfig(network_min_paired_finite_observations=3)
            ),
            "config.mode='production' requires",
        ),
    ],
)
def test_signalome_config_self_validates(factory: object, pattern: str) -> None:
    assert callable(factory)
    with pytest.raises(ContractValidationError, match=pattern):
        factory()


@pytest.mark.parametrize(
    ("factory", "invalid_pattern", "wrong_type_pattern"),
    [
        pytest.param(
            lambda value: KinaseActivityConfig(threshold=value),
            "activity_config.threshold must be between 0.0 and 1.0",
            "activity_config.threshold must be a float between 0.0 and 1.0",
            id="activity-threshold",
        ),
        pytest.param(
            lambda value: KinaseActivityConfig(ksea_evidence_threshold=value),
            "activity_config.ksea_evidence_threshold must be between 0.0 and 1.0",
            "activity_config.ksea_evidence_threshold must be a float between 0.0 and 1.0",
            id="activity-ksea-evidence-threshold",
        ),
        pytest.param(
            lambda value: SignalomeScientificConfig(substrate_support_cutoff=value),
            "signalome workflow request config.scientific.substrate_support_cutoff must be between 0.0 and 1.0",
            "signalome workflow request config.scientific.substrate_support_cutoff must be a float between 0.0 and 1.0",
            id="signalome-substrate-support-cutoff",
        ),
        pytest.param(
            lambda value: SignalomeOutputConfig(network_correlation_threshold=value),
            "signalome workflow request config.output.network_correlation_threshold must be between 0.0 and 1.0",
            "signalome workflow request config.output.network_correlation_threshold must be a float between 0.0 and 1.0",
            id="signalome-network-correlation-threshold",
        ),
        pytest.param(
            lambda value: SignalomeClusteringConfig(
                module_selection_primary_correlation_threshold=value
            ),
            "signalome workflow request config.clustering.module_selection_primary_correlation_threshold must be between 0.0 and 1.0",
            "signalome workflow request config.clustering.module_selection_primary_correlation_threshold must be a float between 0.0 and 1.0",
            id="signalome-module-selection-primary-correlation-threshold",
        ),
        pytest.param(
            lambda value: SignalomeClusteringConfig(
                module_selection_fallback_correlation_threshold=value
            ),
            "signalome workflow request config.clustering.module_selection_fallback_correlation_threshold must be between 0.0 and 1.0",
            "signalome workflow request config.clustering.module_selection_fallback_correlation_threshold must be a float between 0.0 and 1.0",
            id="signalome-module-selection-fallback-correlation-threshold",
        ),
    ],
)
@pytest.mark.parametrize(
    ("value", "expectation"),
    [
        pytest.param(-0.1, "invalid-range", id="below-min"),
        pytest.param(0.0, "valid", id="at-min"),
        pytest.param(0.5, "valid", id="inside-range"),
        pytest.param(1.0, "valid", id="at-max"),
        pytest.param(1.2, "invalid-range", id="above-max"),
        pytest.param(True, "wrong-type", id="wrong-type"),
        pytest.param(float("nan"), "invalid-range", id="nan"),
        pytest.param(float("inf"), "invalid-range", id="infinite"),
    ],
)
def test_probability_like_range_fields_self_validate(
    factory: object,
    value: object,
    expectation: str,
    invalid_pattern: str,
    wrong_type_pattern: str,
) -> None:
    assert callable(factory)
    if expectation == "valid":
        config = factory(value)  # type: ignore[misc]
        assert config is not None
        return

    pattern = wrong_type_pattern if expectation == "wrong-type" else invalid_pattern
    with pytest.raises(ContractValidationError, match=pattern):
        factory(value)  # type: ignore[misc]


def test_signalome_config_accepts_supported_clustering_engine_names() -> None:
    exact = SignalomeConfig(
        clustering=SignalomeClusteringConfig(
            clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
        )
    )
    scipy = SignalomeConfig(
        clustering=SignalomeClusteringConfig(
            clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
        )
    )
    assert (
        exact.clustering.clustering_engine == SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    )
    assert (
        scipy.clustering.clustering_engine
        == SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )


def test_signalome_output_config_accepts_network_minimum_observation_policy() -> None:
    config = SignalomeOutputConfig(
        network_policy="absolute_threshold",
        network_min_paired_finite_observations=3,
    )

    assert config.network_policy == "absolute_threshold"
    assert config.network_min_paired_finite_observations == 3


def test_signalome_config_presets_return_expected_values() -> None:
    default = SignalomeConfig()
    strict = SignalomeConfig.strict()
    permissive = SignalomeConfig.permissive_missing_scores()
    sampled = SignalomeConfig.sampled_candidate_scoring()
    compatibility = SignalomeConfig.compatibility()

    assert default.mode == SIGNALOME_MODE_PRODUCTION
    assert compatibility.mode == SIGNALOME_MODE_EXPLORATORY_COMPATIBILITY
    assert (
        strict.validation.score_preconditioning_policy
        == SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
    )
    assert (
        default.validation.score_preconditioning_policy
        == SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
    )
    assert (
        permissive.validation.score_preconditioning_policy
        == SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )
    assert (
        sampled.clustering.candidate_scoring_policy
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert (
        sampled.performance.max_exact_tree_sites
        == SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT
    )
    assert (
        sampled.performance.max_full_candidate_scoring_sites
        == SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT
    )
    assert (
        default.output.network_min_paired_finite_observations
        == SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT
    )
    assert (
        compatibility.output.network_min_paired_finite_observations
        < default.output.network_min_paired_finite_observations
    )
    assert (
        compatibility.validation.localisation_requirement == LocalisationRequirement()
    )
    assert SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT == 5
    assert (
        default.validation.reference_context_compatibility_policy
        == REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH
    )


def test_signalome_production_config_uses_strict_localisation() -> None:
    production = SignalomeConfig.production()

    assert production.validation.localisation_requirement == (
        LocalisationRequirement.production_site_level()
    )
    assert (
        production.validation.score_preconditioning_policy
        == SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
    )
    assert SignalomeConfig().validation.localisation_requirement == (
        LocalisationRequirement.production_site_level()
    )
    assert production.mode == SIGNALOME_MODE_PRODUCTION


def test_signalome_config_removes_large_dataset_preset() -> None:
    assert not hasattr(SignalomeConfig, "large_dataset")


def test_signalome_config_defaults_to_scipy_clustering_engine() -> None:
    assert (
        SignalomeConfig().clustering.clustering_engine
        == SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )


def test_signalome_config_rejects_removed_max_exact_clustering_sites_alias() -> None:
    with pytest.raises(
        TypeError,
        match="unexpected keyword argument 'max_exact_clustering_sites'",
    ):
        SignalomeConfig(max_exact_clustering_sites=1234)  # type: ignore[call-arg]


def test_signalome_config_accepts_engine_and_policy_names() -> None:
    config = SignalomeConfig(
        clustering=SignalomeClusteringConfig(
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
            clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
        )
    )
    assert (
        config.clustering.candidate_scoring_policy
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert (
        config.clustering.clustering_engine
        == SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )


def test_signalome_config_default_preset_matches_full_candidate_scoring() -> None:
    assert (
        SignalomeConfig.strict().clustering.candidate_scoring_policy
        == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    )


def test_signalome_clustering_config_rejects_removed_tree_engine_argument() -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument 'tree_engine'"):
        SignalomeClusteringConfig(tree_engine="exact")  # type: ignore[call-arg]


@pytest.mark.parametrize(
    "removed_name",
    [
        "cluster_tree_backend",
        "candidate_scoring_backend",
        "clustering_backend",
        "max_exact_cluster_tree_sites",
        "max_full_correlation_sites",
    ],
)
def test_signalome_config_rejects_removed_backend_style_names(
    removed_name: str,
) -> None:
    with pytest.raises(
        TypeError,
        match=f"unexpected keyword argument '{removed_name}'",
    ):
        SignalomeConfig(**{removed_name: "full"})  # type: ignore[arg-type]
