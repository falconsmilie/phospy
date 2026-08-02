from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.configs import (
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetProteinAwarePreparationConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
    DatasetTotalProteinCorrectionIdentityConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
    PreprocessingPlan,
)
from phospy.science.datasets.preprocessing.plan_interpreter import (
    PreprocessingPlanInterpreter,
)
from phospy.science.datasets.preprocessing.plan_rules import (
    PreprocessingBatchCorrectionPlanRuleFamily,
    PreprocessingGroupCoveragePlanRuleFamily,
    PreprocessingLocalisationPlanRuleFamily,
)
from phospy.science.datasets.preprocessing.plan_stage_order import (
    PreprocessingStageOrderPlanner,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    IntensityTransformPolicy,
    LocalisationEligibilityMode,
    MissingDataPolicy,
    NormalisationPolicy,
    SiteMatrixPolicy,
    TotalProteinCorrectionPolicy,
)


def test_preprocessing_plan_interpreter_builds_default_plan() -> None:
    plan = PreprocessingPlanInterpreter().run(DatasetPreprocessingConfig())

    assert plan.intensity_transform_policy == "identity"
    assert plan.normalisation_policy == "none"
    assert plan.stage_order == ("localisation_confidence", "missing_data")
    assert "intensity_transform" not in plan.stage_order
    assert "normalisation" not in plan.stage_order


def test_preprocessing_plan_interpreter_orders_minprob_after_log2_transform() -> None:
    plan = PreprocessingPlanInterpreter().run(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            missing_data=DatasetMissingDataConfig(
                policy="impute_minprob",
                q=0.01,
                width=0.3,
                seed=12345,
                max_missing_fraction_per_row=0.5,
            ),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
        )
    )

    assert plan.stage_order == (
        "localisation_confidence",
        "intensity_transform",
        "missing_data",
        "site_matrix",
    )
    assert plan.stage_order_resolution[2].stage == "missing_data"
    assert (
        plan.stage_order_resolution[2].rationale
        == PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA
    )


def test_preprocessing_plan_interpreter_preserves_minprob_scale_error() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing_data.policy='impute_minprob' requires",
    ):
        PreprocessingPlanInterpreter().run(
            DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity"),
                missing_data=DatasetMissingDataConfig(
                    policy="impute_minprob",
                    q=0.01,
                    width=0.3,
                    seed=12345,
                    max_missing_fraction_per_row=0.5,
                ),
            )
        )


def test_preprocessing_plan_interpreter_resolves_total_mapping_table_identity() -> None:
    mapping_table = pd.DataFrame(
        {
            "phosphosite_id": [" MAPK14 ", pd.NA],
            "protein_id": ["MAPK14_TOTAL ", " AKT1_TOTAL"],
        }
    )

    plan = PreprocessingPlanInterpreter().run(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total",
                identity=DatasetTotalProteinCorrectionIdentityConfig(
                    mode="mapping_table",
                    mapping_table=mapping_table,
                    mapping_phosphosite_key="phosphosite_id",
                    mapping_total_protein_key="protein_id",
                ),
            ),
        )
    )

    identity = plan.total_protein_correction_identity_policy
    assert identity.mode == "mapping_table"
    assert identity.mapping_table == (
        ("MAPK14", "MAPK14_TOTAL"),
        ("", "AKT1_TOTAL"),
    )
    assert identity.mapping_phosphosite_key == "phosphosite_id"
    assert identity.mapping_total_protein_key == "protein_id"
    assert isinstance(identity.mapping_table_fingerprint, str)
    assert len(identity.mapping_table_fingerprint) == 64


def test_preprocessing_plan_interpreter_rejects_invalid_protein_aware_config_before_plan_construction() -> (
    None
):
    class PlanConstructionShouldNotRun:
        def __init__(self, *args: object, **kwargs: object) -> None:
            raise AssertionError("plan construction should not run")

    protein_aware_preparation = DatasetProteinAwarePreparationConfig()
    object.__setattr__(
        protein_aware_preparation,
        "policy",
        "joint_ptm_protein_model",
    )
    config = DatasetPreprocessingConfig(
        protein_aware_preparation=protein_aware_preparation
    )

    with pytest.raises(PhosPyInputError, match="protein_aware_preparation.policy"):
        PreprocessingPlanInterpreter(
            plan_type=PlanConstructionShouldNotRun,  # type: ignore[arg-type]
        ).run(config)


def test_preprocessing_plan_from_config_remains_compatibility_wrapper() -> None:
    config = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
        missing_data=DatasetMissingDataConfig(
            policy="impute_row_median",
            min_observed_values=1,
        ),
    )

    assert PreprocessingPlan.from_config(config) == PreprocessingPlanInterpreter().run(
        config
    )


def test_preprocessing_stage_order_planner_returns_structured_findings() -> None:
    stage_plan = PreprocessingStageOrderPlanner().run(
        site_sequence_resolution_enabled=False,
        intensity_transform_policy=IntensityTransformPolicy.LOG2,
        normalisation_policy=NormalisationPolicy.MEDIAN_CENTER,
        site_matrix_policy=SiteMatrixPolicy.BUILD_FROM_METADATA,
        comparison_building_policy=ComparisonBuildingPolicy.NONE,
        localisation_mode=LocalisationEligibilityMode.REQUIRE_THRESHOLD,
        missing_data_policy=MissingDataPolicy.IMPUTE_ROW_MEDIAN,
        batch_correction_method="linear_residualize_batch",
        total_correction_policy=TotalProteinCorrectionPolicy.NONE,
        group_coverage_filter_enabled=True,
    )

    assert stage_plan.stage_order == (
        "localisation_confidence",
        "group_coverage_filter",
        "missing_data",
        "intensity_transform",
        "batch_correction",
        "site_matrix",
        "normalisation",
    )
    assert [row.order_index for row in stage_plan.stage_order_resolution] == list(
        range(len(stage_plan.stage_order))
    )
    assert stage_plan.stage_order_resolution[1].rationale.startswith("group-aware")


def test_preprocessing_localisation_rule_family_normalises_waiver() -> None:
    resolved = PreprocessingLocalisationPlanRuleFamily().run(
        localisation_mode=LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER,
        localisation_min_confidence=0.75,
        localisation_confidence_column=" localisation_confidence ",
        localisation_waiver_reason=" documented import limitation ",
    )

    assert (
        resolved.localisation_mode
        is LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER
    )
    assert resolved.localisation_confidence_column == "localisation_confidence"
    assert resolved.localisation_waiver_reason == "documented import limitation"


def test_preprocessing_batch_rule_family_rejects_duplicate_condition_columns() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="batch_correction_condition_columns .* must not contain duplicates",
    ):
        PreprocessingBatchCorrectionPlanRuleFamily().run(
            batch_correction_method="linear_residualize_batch",
            batch_correction_batch_column="batch",
            batch_correction_condition_column="condition",
            batch_correction_condition_columns=("condition", "condition"),
            batch_correction_replicate_column=None,
            batch_correction_control_site_set=None,
            batch_correction_missingness_policy=None,
            batch_correction_internal_request=None,
            batch_correction_preserve_condition_effects=True,
        )


def test_preprocessing_group_coverage_rule_family_rejects_missing_stage() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="stage_order does not include 'group_coverage_filter'",
    ):
        PreprocessingGroupCoveragePlanRuleFamily().run(
            enabled=True,
            group_column="condition",
            min_finite_observations_per_group=1,
            min_finite_fraction_per_group=None,
            min_groups_passing_threshold=1,
            stage_order=("missing_data",),
        )
