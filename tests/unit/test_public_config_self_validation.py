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
    KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT,
    SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT,
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
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
)
from phospy.errors import PhosPyInputError, WorkflowValidationError


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
        (
            {"policy": "log2", "pseudocount": "1.0"},
            "intensity_transform.pseudocount must be a float or int",
        ),
        (
            {"policy": "log2", "pseudocount": True},
            "intensity_transform.pseudocount must be a float or int",
        ),
        (
            {"policy": "log2", "pseudocount": float("nan")},
            "intensity_transform.pseudocount must be finite",
        ),
        (
            {"policy": "log2", "pseudocount": -0.1},
            "intensity_transform.pseudocount must be greater than or equal to 0",
        ),
    ],
)
def test_dataset_intensity_transform_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(PhosPyInputError, match=pattern):
        DatasetIntensityTransformConfig(**kwargs)  # type: ignore[arg-type]


def test_dataset_normalisation_config_rejects_unknown_policy() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="preprocessing_config.normalisation.policy must be one of",
    ):
        DatasetNormalisationConfig(policy="unsupported")  # type: ignore[arg-type]


def test_dataset_total_protein_correction_config_rejects_unknown_policy() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="preprocessing_config.total_protein_correction.policy must be one of",
    ):
        DatasetTotalProteinCorrectionConfig(policy="unsupported")  # type: ignore[arg-type]


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
            {"profile_missing_value_strategy": "invalid"},
            "scoring_config.profile_missing_value_strategy must be one of",
        ),
        (
            {"allow_mixed_total_protein_quantitative_meaning": "yes"},
            "scoring_config.allow_mixed_total_protein_quantitative_meaning must be a bool",
        ),
    ],
)
def test_kinase_scoring_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(WorkflowValidationError, match=pattern):
        KinaseScoringConfig(**kwargs)  # type: ignore[arg-type]


def test_kinase_scoring_config_presets_return_expected_values() -> None:
    default = KinaseScoringConfig.default()
    strict_missing = KinaseScoringConfig.strict_missing_values()

    assert (
        default.profile_missing_value_strategy
        == KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    )
    assert (
        strict_missing.profile_missing_value_strategy
        == KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
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
        (
            {"random_state": -1},
            "prediction_config.random_state must be greater than or equal to 0",
        ),
        (
            {"random_state": True},
            "prediction_config.random_state must be an int",
        ),
    ],
)
def test_kinase_prediction_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(WorkflowValidationError, match=pattern):
        KinasePredictionConfig(**kwargs)  # type: ignore[arg-type]


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
        WorkflowValidationError,
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
            {"threshold": 1.2},
            "activity_config.threshold must be between 0.0 and 1.0",
        ),
        (
            {"threshold": True},
            "activity_config.threshold must be a float between 0.0 and 1.0",
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
            {"ksea_evidence_threshold": 1.2},
            "activity_config.ksea_evidence_threshold must be between 0.0 and 1.0",
        ),
        (
            {"ksea_p_value_method": "invalid"},
            "activity_config.ksea_p_value_method must be one of",
        ),
        (
            {"ksea_adjust_p_values": "yes"},
            "activity_config.ksea_adjust_p_values must be a bool",
        ),
    ],
)
def test_kinase_activity_config_self_validates(
    kwargs: dict[str, object], pattern: str
) -> None:
    with pytest.raises(WorkflowValidationError, match=pattern):
        KinaseActivityConfig(**kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("factory", "pattern"),
    [
        (
            lambda: SignalomeScientificConfig(substrate_support_cutoff=1.5),
            "signalome workflow request config.scientific.substrate_support_cutoff",
        ),
        (
            lambda: SignalomeScientificConfig(substrate_support_cutoff=True),  # type: ignore[arg-type]
            "signalome workflow request config.scientific.substrate_support_cutoff",
        ),
        (
            lambda: SignalomeOutputConfig(network_correlation_threshold=-0.1),
            "signalome workflow request config.output.network_correlation_threshold",
        ),
        (
            lambda: SignalomeOutputConfig(network_policy="invalid"),  # type: ignore[arg-type]
            "signalome workflow request config.output.network_policy",
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
            lambda: SignalomeClusteringConfig(module_count=0),
            "signalome workflow request config.clustering.module_count",
        ),
        (
            lambda: SignalomeClusteringConfig(module_count=-1),
            "signalome workflow request config.clustering.module_count",
        ),
        (
            lambda: SignalomeClusteringConfig(module_count=True),  # type: ignore[arg-type]
            "signalome workflow request config.clustering.module_count must be an int",
        ),
        (
            lambda: SignalomeClusteringConfig(module_count=2.5),  # type: ignore[arg-type]
            "signalome workflow request config.clustering.module_count must be an int",
        ),
        (
            lambda: SignalomeClusteringConfig(
                module_selection_primary_correlation_threshold=1.2
            ),
            "module_selection_primary_correlation_threshold",
        ),
        (
            lambda: SignalomeClusteringConfig(
                module_selection_fallback_correlation_threshold=-0.1
            ),
            "module_selection_fallback_correlation_threshold",
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
    ],
)
def test_signalome_config_self_validates(factory: object, pattern: str) -> None:
    assert callable(factory)
    with pytest.raises(WorkflowValidationError, match=pattern):
        factory()


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


def test_signalome_config_presets_return_expected_values() -> None:
    default = SignalomeConfig()
    strict = SignalomeConfig.strict()
    permissive = SignalomeConfig.permissive_missing_scores()
    sampled = SignalomeConfig.sampled_candidate_scoring()

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
