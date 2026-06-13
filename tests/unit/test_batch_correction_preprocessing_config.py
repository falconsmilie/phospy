from __future__ import annotations

import pytest

from phospy.api import DatasetBatchCorrectionConfig
from phospy.api.configs import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.errors import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    list_registered_preprocessing_stages,
)


def test_default_batch_correction_config_disables_correction() -> None:
    config = DatasetBatchCorrectionConfig()

    assert config.method == DATASET_BATCH_CORRECTION_METHOD_NONE
    assert config.batch_column == "batch"
    assert config.condition_column == "condition"
    assert config.preserve_condition_effects is True


def test_linear_residualize_batch_can_be_declared() -> None:
    config = DatasetBatchCorrectionConfig(
        method=DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
        batch_column="ms_run",
        condition_column="treatment",
        preserve_condition_effects=True,
    )

    assert config.method == "linear_residualize_batch"
    assert config.batch_column == "ms_run"
    assert config.condition_column == "treatment"
    assert config.preserve_condition_effects is True


def test_batch_correction_config_rejects_unsupported_method() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="preprocessing_config.batch_correction.method must be one of",
    ):
        DatasetBatchCorrectionConfig(
            method="unsupported_method"  # type: ignore[arg-type]
        )


def test_batch_correction_config_preserves_condition_effects_by_design() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="preserve_condition_effects must be True",
    ):
        DatasetBatchCorrectionConfig(
            method=DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
            preserve_condition_effects=False,  # type: ignore[arg-type]
        )


def test_batch_correction_config_appears_inside_dataset_preprocessing_config() -> None:
    config = DatasetPreprocessingConfig(
        batch_correction=DatasetBatchCorrectionConfig(
            method=DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH
        )
    )

    assert isinstance(config.batch_correction, DatasetBatchCorrectionConfig)
    assert config.batch_correction.method == "linear_residualize_batch"


def test_existing_preprocessing_config_construction_remains_backward_compatible() -> (
    None
):
    config = DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(policy="forbid")
    )

    assert config.missing_data.policy == "forbid"
    assert config.batch_correction.method == DATASET_BATCH_CORRECTION_METHOD_NONE
    assert (
        DatasetPreprocessingConfig.strict().batch_correction.method
        == DATASET_BATCH_CORRECTION_METHOD_NONE
    )
    assert (
        DatasetPreprocessingConfig.from_raw_phosphosite_table().batch_correction.method
        == DATASET_BATCH_CORRECTION_METHOD_NONE
    )


def test_declared_batch_correction_is_scheduled_after_intensity_transform() -> None:
    plan = PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            batch_correction=DatasetBatchCorrectionConfig(
                method=DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH
            ),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(policy="none"),
        )
    )

    assert DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION in plan.stage_order
    assert {
        stage.stage_key for stage in list_registered_preprocessing_stages()
    }.issuperset({DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION})
    assert plan.stage_order.index(DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION) > (
        plan.stage_order.index(DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM)
    )
    assert plan.stage_order.index(DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION) > (
        plan.stage_order.index(DATASET_PREPROCESSING_STAGE_MISSING_DATA)
    )
    assert DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION not in plan.stage_order
