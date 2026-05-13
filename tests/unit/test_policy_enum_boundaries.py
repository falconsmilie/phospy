from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteSequenceResolutionConfig,
    DatasetTotalProteinCorrectionConfig,
    KinaseScoringConfig,
    LocalisationRequirement,
    SignalomeValidationConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import WorkflowValidationError
from phospy.science.activities.threshold_membership import (
    threshold_membership_mask_array,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    IntensityTransformPolicy,
    MissingDataPolicy,
    NormalisationPolicy,
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixPolicy,
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    resolve_builder_provenance_stage_order,
)
from phospy.science.prediction.candidates import build_candidate_substrate_list
from phospy.science.prediction.scoring import select_downstream_score_matrix
from phospy.science.scoring.policy_models import DownstreamScoreSource, ThresholdMode


def test_preprocessing_plan_converts_public_strings_to_internal_policy_enums() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(policy="forbid"),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(policy="none"),
        )
    )

    assert isinstance(plan.missing_data_policy, MissingDataPolicy)
    assert plan.missing_data_policy is MissingDataPolicy.FORBID
    assert isinstance(
        plan.total_protein_correction_policy, TotalProteinCorrectionPolicy
    )
    assert plan.total_protein_correction_policy is TotalProteinCorrectionPolicy.NONE


@pytest.mark.parametrize(
    ("field_name", "value", "expected"),
    [
        ("intensity_transform_policy", "log2", IntensityTransformPolicy.LOG2),
        ("normalisation_policy", "quantile", NormalisationPolicy.QUANTILE),
        (
            "site_matrix_policy",
            "build_from_metadata",
            SiteMatrixPolicy.BUILD_FROM_METADATA,
        ),
        (
            "site_matrix_duplicate_site_policy",
            "aggregate_mean",
            SiteMatrixDuplicateSitePolicy.AGGREGATE_MEAN,
        ),
        (
            "site_sequence_resolution_mode",
            "fill_missing_only",
            SiteSequenceResolutionMode.FILL_MISSING_ONLY,
        ),
        (
            "site_sequence_resolution_conflict_policy",
            "replace_existing",
            SiteSequenceConflictPolicy.REPLACE_EXISTING,
        ),
        (
            "comparison_building_policy",
            "sample_metadata_pairs",
            ComparisonBuildingPolicy.SAMPLE_METADATA_PAIRS,
        ),
    ],
)
def test_preprocessing_plan_converts_selected_policy_strings(
    field_name: str,
    value: str,
    expected: object,
) -> None:
    plan = PreprocessingPlan(**{field_name: value})  # type: ignore[arg-type]
    assert getattr(plan, field_name) is expected


@pytest.mark.parametrize(
    "field_name",
    [
        "intensity_transform_policy",
        "normalisation_policy",
        "site_matrix_policy",
        "site_matrix_duplicate_site_policy",
        "site_sequence_resolution_mode",
        "site_sequence_resolution_conflict_policy",
        "comparison_building_policy",
    ],
)
def test_preprocessing_plan_rejects_unknown_selected_policy_strings(
    field_name: str,
) -> None:
    with pytest.raises(PhosPyInputError, match="must be one of:"):
        PreprocessingPlan(**{field_name: "invalid_policy"})  # type: ignore[arg-type]


def test_from_config_converts_selected_policies_to_internal_enums() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            normalisation=DatasetNormalisationConfig(policy="median_center"),
            site_matrix=DatasetSiteMatrixConfig(
                policy="build_from_metadata",
                duplicate_site_policy="aggregate_median",
            ),
            site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                mode="replace_existing",
                conflict_policy="error",
            ),
            comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs"),
            missing_data=DatasetMissingDataConfig(policy="forbid"),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(policy="none"),
        )
    )

    assert plan.intensity_transform_policy is IntensityTransformPolicy.LOG2
    assert plan.normalisation_policy is NormalisationPolicy.MEDIAN_CENTER
    assert plan.site_matrix_policy is SiteMatrixPolicy.BUILD_FROM_METADATA
    assert (
        plan.site_matrix_duplicate_site_policy
        is SiteMatrixDuplicateSitePolicy.AGGREGATE_MEDIAN
    )
    assert (
        plan.site_sequence_resolution_mode
        is SiteSequenceResolutionMode.REPLACE_EXISTING
    )
    assert (
        plan.site_sequence_resolution_conflict_policy
        is SiteSequenceConflictPolicy.ERROR
    )
    assert (
        plan.comparison_building_policy
        is ComparisonBuildingPolicy.SAMPLE_METADATA_PAIRS
    )


def test_stage_registry_serializes_selected_policy_values_as_stable_strings() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            normalisation=DatasetNormalisationConfig(policy="quantile"),
            site_matrix=DatasetSiteMatrixConfig(
                policy="build_from_metadata",
                duplicate_site_policy="first",
            ),
            site_sequence_resolution=DatasetSiteSequenceResolutionConfig(
                fasta_path="test.fasta",
                mode="fill_missing_only",
                conflict_policy="preserve_existing",
            ),
            comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs"),
        )
    )

    metadata_by_key = {
        item.stage_key: item for item in resolve_builder_provenance_stage_order(plan)
    }
    assert (
        metadata_by_key["intensity_transform"].operation_name(plan)
        == IntensityTransformPolicy.LOG2.value
    )
    assert (
        metadata_by_key["normalisation"].operation_name(plan)
        == NormalisationPolicy.QUANTILE.value
    )
    assert (
        metadata_by_key["site_matrix"].operation_name(plan)
        == SiteMatrixPolicy.BUILD_FROM_METADATA.value
    )
    assert (
        metadata_by_key["comparisons"].operation_name(plan)
        == ComparisonBuildingPolicy.SAMPLE_METADATA_PAIRS.value
    )
    site_sequence_parameters = metadata_by_key[
        "site_sequence_resolution"
    ].serialize_parameters(plan)
    assert (
        site_sequence_parameters["mode"]
        == SiteSequenceResolutionMode.FILL_MISSING_ONLY.value
    )
    assert (
        site_sequence_parameters["conflict_policy"]
        == SiteSequenceConflictPolicy.PRESERVE_EXISTING.value
    )


def test_selected_internal_policy_branches_do_not_compare_raw_string_literals() -> None:
    targets = (
        "src/phospy/science/datasets/preprocessing/stages/intensity_transform.py",
        "src/phospy/science/datasets/preprocessing/stages/normalisation.py",
        "src/phospy/science/datasets/preprocessing/stages/comparisons.py",
        "src/phospy/science/datasets/preprocessing/stages/site_matrix.py",
        "src/phospy/science/datasets/preprocessing/stages/site_sequence_resolution.py",
        "src/phospy/science/datasets/preprocessing/stages/missing_data/minprob.py",
    )
    forbidden_tokens = (
        '== "identity"',
        '!= "identity"',
        '== "log2"',
        '!= "log2"',
        '== "none"',
        '!= "none"',
        '== "median_center"',
        '== "quantile"',
        '== "as_input"',
        '== "build_from_metadata"',
        '== "sample_metadata_pairs"',
        '== "fill_missing_only"',
        '== "replace_existing"',
        '== "preserve_existing"',
    )
    for path in targets:
        content = Path(path).read_text(encoding="utf-8")
        for token in forbidden_tokens:
            assert token not in content, (
                f"raw-string policy branch found in {path}: {token}"
            )


def test_preprocessing_plan_rejects_unknown_policy_strings_without_defaulting() -> None:
    with pytest.raises(PhosPyInputError, match="must be one of:"):
        PreprocessingPlan(
            missing_data_policy="unknown_policy",  # type: ignore[arg-type]
            total_protein_correction_policy=TotalProteinCorrectionPolicy.NONE,
        )

    with pytest.raises(PhosPyInputError, match="must be one of:"):
        PreprocessingPlan(
            missing_data_policy=MissingDataPolicy.FORBID,
            total_protein_correction_policy="unknown_policy",  # type: ignore[arg-type]
        )


def test_downstream_score_selection_returns_enum_source() -> None:
    profile = pd.DataFrame({"K1": [0.1]}, index=["S1"])
    combined = pd.DataFrame({"K1": [0.2]}, index=["S1"])

    selected, source = select_downstream_score_matrix(
        profile_scores=profile,
        rank_weighted_fusion_scores=combined,
    )

    assert selected is combined
    assert source is DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES


def test_threshold_mode_drives_membership_logic_with_enum_or_string_input() -> None:
    scores = np.array([0.5, 0.4, np.nan])
    strict = threshold_membership_mask_array(
        scores,
        threshold=0.5,
        threshold_mode=ThresholdMode.GREATER_THAN,
    )
    inclusive = threshold_membership_mask_array(
        scores,
        threshold=0.5,
        threshold_mode="score >= threshold",
    )

    assert strict.tolist() == [False, False, False]
    assert inclusive.tolist() == [True, False, False]


def test_candidate_selection_rejects_unknown_threshold_mode() -> None:
    scores = pd.DataFrame({"K1": [0.7]}, index=["S1"])
    with pytest.raises(PhosPyInputError, match="must be one of:"):
        build_candidate_substrate_list(
            scores=scores,
            top=1,
            score_threshold=0.0,
            inclusion=1,
            threshold_mode="invalid",  # type: ignore[arg-type]
        )


def test_localisation_requirement_default_policy_is_allow_unknown() -> None:
    requirement = LocalisationRequirement()
    assert requirement.policy == "allow_unknown"
    assert requirement.requires_probability_column is False


def test_localisation_requirement_threshold_policy_requires_probability_column() -> (
    None
):
    requirement = LocalisationRequirement(minimum_probability=0.75)
    assert requirement.policy == "require_threshold"
    assert requirement.requires_probability_column is True


@pytest.mark.parametrize("invalid_threshold", [-0.1, 1.2, "high"])
def test_localisation_requirement_rejects_invalid_threshold(
    invalid_threshold: object,
) -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="localisation_requirement.minimum_probability",
    ):
        LocalisationRequirement(
            minimum_probability=invalid_threshold  # type: ignore[arg-type]
        )


def test_kinase_scoring_config_accepts_localisation_requirement() -> None:
    config = KinaseScoringConfig(
        min_substrates=2,
        localisation_requirement=LocalisationRequirement(require_present=True),
    )
    assert config.localisation_requirement.policy == "require_present"


def test_signalome_validation_config_accepts_localisation_requirement() -> None:
    config = SignalomeValidationConfig(
        localisation_requirement=LocalisationRequirement(minimum_probability=0.6)
    )
    assert config.localisation_requirement.policy == "require_threshold"
