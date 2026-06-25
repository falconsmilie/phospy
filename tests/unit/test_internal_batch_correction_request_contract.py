from __future__ import annotations

import pytest

from phospy.api.configs import DATASET_BATCH_CORRECTION_METHODS
from phospy.contracts.configs.preprocessing import (
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
)
from phospy.errors import PhosPyInputError


def test_internal_batch_correction_request_accepts_typed_enum_values() -> None:
    request = InternalBatchCorrectionRequest(
        method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
        batch_column="ms_run",
        condition_columns=("condition", "timepoint"),
        replicate_column="bio_rep",
        control_site_source=(InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED),
        control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
        missing_value_policy=(InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING),
        imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
        n_unwanted_factors=2,
        stage_order=(
            InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
        ),
        diagnostics_enabled=True,
    )

    assert request.method is InternalBatchCorrectionMethod.SPS_RUV_STYLE
    assert request.batch_column == "ms_run"
    assert request.condition_columns == ("condition", "timepoint")
    assert request.n_unwanted_factors == 2
    assert request.diagnostics_enabled is True


def test_internal_batch_correction_request_coerces_supported_strings() -> None:
    request = InternalBatchCorrectionRequest(
        method="sps_ruv_style",  # type: ignore[arg-type]
        batch_column="batch",
        condition_columns=["condition"],  # type: ignore[arg-type]
        replicate_column=None,
        control_site_source="dataset_metadata",  # type: ignore[arg-type]
        control_site_mode="boolean_column",  # type: ignore[arg-type]
        missing_value_policy="allow_temporary_imputation",  # type: ignore[arg-type]
        imputation_policy="row_median_temporary",  # type: ignore[arg-type]
        n_unwanted_factors=None,
        stage_order=(  # type: ignore[arg-type]
            "after_intensity_transform_before_missing_data"
        ),
        diagnostics_enabled=False,
    )

    assert request.method is InternalBatchCorrectionMethod.SPS_RUV_STYLE
    assert request.condition_columns == ("condition",)
    assert (
        request.missing_value_policy
        is InternalBatchCorrectionMissingValuePolicy.ALLOW_TEMPORARY_IMPUTATION
    )
    assert (
        request.imputation_policy
        is InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY
    )


def test_internal_control_site_ruv_style_alias_normalizes_with_deprecation() -> None:
    with pytest.warns(
        DeprecationWarning,
        match="control_site_ruv_style is a deprecated internal compatibility alias",
    ):
        request = InternalBatchCorrectionRequest(
            method="control_site_ruv_style",  # type: ignore[arg-type]
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column=None,
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
            n_unwanted_factors=None,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        )

    assert request.method is InternalBatchCorrectionMethod.SPS_RUV_STYLE


def test_internal_batch_correction_request_rejects_invalid_method() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "internal batch-correction request.method must be one of: "
            "sps_ruv_style, control_site_ruv_style, ruv_iii_style; got 'combat'"
        ),
    ):
        InternalBatchCorrectionRequest(
            method="combat",  # type: ignore[arg-type]
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column=None,
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
            n_unwanted_factors=None,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        )


def test_internal_batch_correction_request_rejects_invalid_control_source() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="internal batch-correction request.control_site_source must be one of:",
    ):
        InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column=None,
            control_site_source="online_lookup",  # type: ignore[arg-type]
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
            n_unwanted_factors=None,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        )


def test_internal_batch_correction_request_rejects_bad_local_values() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="internal batch-correction request.condition_columns must not contain duplicates",
    ):
        InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
            batch_column="batch",
            condition_columns=("condition", "condition"),
            replicate_column=None,
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
            n_unwanted_factors=None,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "internal batch-correction request.n_unwanted_factors must be "
            "greater than or equal to 1"
        ),
    ):
        InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column=None,
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
            n_unwanted_factors=0,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        )


def test_internal_batch_correction_request_rejects_imputation_inconsistency() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "internal batch-correction request.imputation_policy must be 'none' "
            "when missing_value_policy='reject_missing'"
        ),
    ):
        InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column=None,
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=(
                InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY
            ),
            n_unwanted_factors=None,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        )


def test_internal_batch_correction_request_rejects_ruv_iii_style() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="replicate-aware RUV-III numerical semantics are not implemented",
    ):
        InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod.RUV_III_STYLE,
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column="replicate",
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
            n_unwanted_factors=None,
            stage_order=(
                InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
            ),
            diagnostics_enabled=True,
        )


def test_internal_batch_correction_request_rejects_unsupported_stage_order() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="internal batch-correction request.stage_order must be one of:",
    ):
        InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
            batch_column="batch",
            condition_columns=("condition",),
            replicate_column=None,
            control_site_source=(
                InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED
            ),
            control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
            n_unwanted_factors=None,
            stage_order="after_differential_workflow",  # type: ignore[arg-type]
            diagnostics_enabled=True,
        )


def test_public_batch_correction_methods_remain_unchanged() -> None:
    assert DATASET_BATCH_CORRECTION_METHODS == {
        "none",
        "linear_residualize_batch",
    }
