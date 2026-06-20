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
