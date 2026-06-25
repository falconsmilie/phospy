from __future__ import annotations

import inspect

import pytest

from phospy.api import (
    ControlSiteSet,
    CorrectionMissingnessPolicy,
    DatasetBatchCorrectionConfig,
    SpsRuvBatchCorrectionConfig,
)
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
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
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


def test_dataset_preprocessing_config_docstring_describes_executable_correction() -> (
    None
):
    docstring = inspect.getdoc(DatasetPreprocessingConfig) or ""
    normalized = " ".join(docstring.lower().split())

    assert "config-only" not in normalized
    assert "no correction is executed" not in normalized
    assert "optional executable preprocessing correction" in normalized
    assert "spsruvbatchcorrectionconfig" in normalized
    assert "readiness reporting" in normalized


def test_sps_ruv_config_docstring_marks_replicates_as_provenance_only() -> None:
    docstring = inspect.getdoc(SpsRuvBatchCorrectionConfig) or ""
    normalized = " ".join(docstring.lower().split())

    assert "replicate_column" in normalized
    assert "provenance and diagnostics only" in normalized
    assert "does not enable replicate-aware ruv-iii" in normalized


def test_sps_ruv_batch_correction_config_requires_explicit_public_contract() -> None:
    config = SpsRuvBatchCorrectionConfig(
        control_site_set=ControlSiteSet.from_site_keys(("site_a", "site_c")),
        batch_column="ms_run",
        condition_columns=("condition",),
        replicate_column="replicate",
        missingness_policy=CorrectionMissingnessPolicy(),
        n_unwanted_factors=1,
        diagnostics_enabled=True,
        provenance_enabled=True,
    )

    preprocessing = DatasetPreprocessingConfig(batch_correction=config)
    plan = PreprocessingPlan.from_config(preprocessing)

    assert preprocessing.batch_correction is config
    assert plan.batch_correction_method == "sps_ruv_style"
    assert plan.batch_correction_internal_request is not None
    assert plan.batch_correction_internal_request.n_unwanted_factors == 1
    assert plan.batch_correction_control_site_set is config.control_site_set
    assert plan.batch_correction_missingness_policy is config.missingness_policy


def test_sps_ruv_batch_correction_config_rejects_ruv_iii_style() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="replicate-aware RUV-III numerical semantics are not implemented",
    ):
        SpsRuvBatchCorrectionConfig(
            control_site_set=ControlSiteSet.from_site_keys(("site_a", "site_c")),
            batch_column="ms_run",
            condition_columns=("condition",),
            replicate_column="replicate",
            missingness_policy=CorrectionMissingnessPolicy(),
            n_unwanted_factors=1,
            method="ruv_iii_style",
        )


def test_sps_ruv_batch_correction_config_has_no_boolean_shortcut() -> None:
    assert "use_ruv" not in inspect.signature(SpsRuvBatchCorrectionConfig).parameters
    assert "use_ruv" not in inspect.signature(DatasetPreprocessingConfig).parameters


def test_sps_ruv_batch_correction_config_rejects_missing_provenance() -> None:
    with pytest.raises(PhosPyInputError, match="provenance_enabled must be True"):
        SpsRuvBatchCorrectionConfig(
            control_site_set=ControlSiteSet.from_site_keys(("site_a", "site_c")),
            batch_column="batch",
            condition_columns=("condition",),
            missingness_policy=CorrectionMissingnessPolicy(),
            n_unwanted_factors=1,
            provenance_enabled=False,  # type: ignore[arg-type]
        )


def test_sps_ruv_batch_correction_config_rejects_missing_controls() -> None:
    with pytest.raises(PhosPyInputError, match="control_site_set must contain"):
        SpsRuvBatchCorrectionConfig(
            control_site_set=ControlSiteSet(),
            batch_column="batch",
            condition_columns=("condition",),
            missingness_policy=CorrectionMissingnessPolicy(),
            n_unwanted_factors=1,
        )


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


def test_batch_correction_plan_rejects_correction_after_downstream_stage() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="unsupported stage_order.*downstream stages have consumed the matrix",
    ):
        PreprocessingPlan(
            batch_correction_method=DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
            stage_order=(
                DATASET_PREPROCESSING_STAGE_COMPARISONS,
                DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
            ),
        )


def test_batch_correction_plan_rejects_boundary_weakening_stage_order() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="weaken the analysis-ready dataset boundary",
    ):
        PreprocessingPlan(
            batch_correction_method=DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
            stage_order=(
                DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
                DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
            ),
        )


def test_batch_correction_plan_rejects_duplicate_stage_order_entries() -> None:
    with pytest.raises(PhosPyInputError, match="duplicate stages.*batch_correction"):
        PreprocessingPlan(
            batch_correction_method=DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
            stage_order=(
                DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
                DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
            ),
        )
