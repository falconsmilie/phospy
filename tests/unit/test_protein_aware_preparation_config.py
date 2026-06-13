from __future__ import annotations

from dataclasses import fields

import pytest

import phospy
import phospy.api as public_api
from phospy.api import (
    DatasetPreprocessingConfig,
    DatasetProteinAwarePreparationConfig,
    DatasetTotalProteinCorrectionConfig,
    PhosPyInputError,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
)


def test_total_protein_subtract_log_total_config_still_works() -> None:
    correction = DatasetTotalProteinCorrectionConfig(policy="subtract_log_total")

    preprocessing = DatasetPreprocessingConfig(
        total_protein_correction=correction,
    )

    assert preprocessing.total_protein_correction is correction
    assert preprocessing.total_protein_correction.policy == "subtract_log_total"
    assert preprocessing.protein_aware_preparation.policy == "disabled"


def test_protein_aware_preparation_config_can_be_declared() -> None:
    preparation = DatasetProteinAwarePreparationConfig(
        policy="prepare_model_inputs",
        protein_mapping_policy="allow_missing_with_report",
    )

    preprocessing = DatasetPreprocessingConfig(
        protein_aware_preparation=preparation,
    )

    assert preprocessing.protein_aware_preparation is preparation
    assert preprocessing.protein_aware_preparation.policy == "prepare_model_inputs"
    assert (
        preprocessing.protein_aware_preparation.protein_mapping_policy
        == "allow_missing_with_report"
    )
    assert preprocessing.total_protein_correction.policy == "none"


def test_protein_aware_preparation_is_disabled_by_default() -> None:
    assert DatasetProteinAwarePreparationConfig().policy == "disabled"
    assert DatasetPreprocessingConfig().protein_aware_preparation.policy == "disabled"


def test_protein_aware_preparation_unsupported_policy_value_is_rejected() -> None:
    with pytest.raises(PhosPyInputError, match="protein_aware_preparation.policy"):
        DatasetProteinAwarePreparationConfig(
            policy="joint_ptm_protein_model"  # type: ignore[arg-type]
        )


def test_protein_aware_preparation_declaration_does_not_execute_any_stage() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs"
            )
        )
    )

    assert "protein_aware_preparation" not in plan.stage_order
    assert DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION not in plan.stage_order
    assert plan.total_protein_correction_policy.value == "none"


def test_public_contract_protein_aware_preparation_is_typed_and_narrow() -> None:
    assert "DatasetProteinAwarePreparationConfig" in public_api.__all__
    assert "DatasetProteinAwarePreparationConfig" not in phospy.__all__
    assert (
        public_api.DatasetProteinAwarePreparationConfig
        is DatasetProteinAwarePreparationConfig
    )
    assert [field.name for field in fields(DatasetProteinAwarePreparationConfig)] == [
        "policy",
        "protein_mapping_policy",
    ]

    with pytest.raises(TypeError, match="unexpected keyword argument"):
        DatasetProteinAwarePreparationConfig(
            policy="prepare_model_inputs",
            modelling_policy="joint_ptm_protein_model",  # type: ignore[call-arg]
        )
