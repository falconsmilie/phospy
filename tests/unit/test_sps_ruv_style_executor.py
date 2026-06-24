from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.contracts.configs.preprocessing import (
    CorrectionMaskPolicy,
    CorrectionMissingnessPolicy,
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
    ObservationMask,
    OriginallyMissingCellTracking,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.science.batch_correction import (
    SPS_RUV_STYLE_EXECUTOR_ID,
    DeterministicSpsRuvStyleExecutor,
    SpsRuvStyleExecutorResult,
)
from phospy.science.datasets.preprocessing.control_sites import ControlSiteSet
from phospy.validation.datasets.batch_correction import ResolvedBatchDesignMetadata
from phospy.validation.workflows.batch_correction import ControlSiteEligibilityValidator
from phospy.workflows.batch_correction import BatchCorrectionPlanInterpreter


def test_sps_ruv_style_executor_produces_stable_expected_corrected_values() -> None:
    result = _run_executor(_phospho())

    expected = pd.DataFrame(
        {
            "sample_1": [12.0, 6.5, 24.0],
            "sample_2": [12.0, 10.5, 24.0],
            "sample_3": [12.0, 6.5, 24.0],
            "sample_4": [12.0, 10.5, 24.0],
        },
        index=_site_index(),
    )
    pdt.assert_frame_equal(
        result.corrected_matrix,
        expected,
        check_exact=False,
        atol=1e-10,
        rtol=0.0,
    )


def test_sps_ruv_style_executor_estimates_factors_from_control_rows_only() -> None:
    phospho = _phospho()
    altered_non_control = phospho.copy(deep=True)
    altered_non_control.loc["site_b", :] = [105.0, -91.0, 108.0, -88.0]

    first = _run_executor(phospho)
    second = _run_executor(altered_non_control)

    pdt.assert_frame_equal(
        first.estimated_unwanted_factors,
        second.estimated_unwanted_factors,
        check_exact=False,
        atol=1e-10,
        rtol=0.0,
    )
    assert first.diagnostics.control_site_count == 2


def test_sps_ruv_style_executor_preserves_protected_condition_effects() -> None:
    phospho = _phospho()

    result = _run_executor(phospho)

    assert _condition_effect(result.corrected_matrix).loc["site_b"] == pytest.approx(
        _condition_effect(phospho).loc["site_b"],
        abs=1e-10,
    )
    assert result.diagnostics.protected_design_terms == (
        "intercept",
        "condition[treated]",
    )


def test_sps_ruv_style_executor_restores_originally_missing_cells_and_status() -> None:
    phospho = _phospho()
    phospho.loc["site_b", "sample_2"] = pd.NA

    result = _run_executor(phospho, missing_cells=(("site_b", "sample_2"),))

    assert pd.isna(result.corrected_matrix.loc["site_b", "sample_2"])
    assert bool(result.output_observation_mask.loc["site_b", "sample_2"]) is False
    assert result.corrected_cell_status.loc["site_b", "sample_2"] == "restored_missing"
    assert result.withheld_cells == (("site_b", "sample_2"),)
    assert result.diagnostics.originally_missing_cell_count == 1
    assert result.diagnostics.withheld_cell_count == 1
    assert result.warnings


def test_sps_ruv_style_executor_returns_diagnostics_warnings_and_provenance() -> None:
    result = _run_executor(_phospho())

    assert isinstance(result, SpsRuvStyleExecutorResult)
    assert result.diagnostics.status == "applied"
    assert result.diagnostics.executor_id == SPS_RUV_STYLE_EXECUTOR_ID
    assert result.diagnostics.estimated_unwanted_factors == 1
    assert result.diagnostics.singular_values[0] > 0.0
    assert result.warnings == ()
    assert result.withheld_rows == ()
    assert result.rejected_rows == ()
    assert result.rejected_cells == ()
    assert result.provenance_payload["executor_id"] == SPS_RUV_STYLE_EXECUTOR_ID
    assert result.provenance_payload["method"] == "sps_ruv_style"
    assert "corrected_matrix_fingerprint" in result.provenance_payload


def _run_executor(
    phospho: pd.DataFrame,
    *,
    missing_cells: tuple[tuple[str, str], ...] = (),
) -> SpsRuvStyleExecutorResult:
    plan = BatchCorrectionPlanInterpreter().run(
        config=_config(),
        dataset_metadata=_metadata(),
        control_site_mapping=_control_mapping(),
        missingness_policy=_missingness_policy(missing_cells=missing_cells),
    )
    return DeterministicSpsRuvStyleExecutor().run(phospho=phospho, plan=plan)


def _config() -> InternalBatchCorrectionRequest:
    return InternalBatchCorrectionRequest(
        method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column="replicate",
        control_site_source=InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
        control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
        missing_value_policy=(
            InternalBatchCorrectionMissingValuePolicy.ALLOW_TEMPORARY_IMPUTATION
        ),
        imputation_policy=InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY,
        n_unwanted_factors=1,
        stage_order=InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM,
        diagnostics_enabled=True,
    )


def _metadata() -> ResolvedBatchDesignMetadata:
    return ResolvedBatchDesignMetadata(
        batch_by_sample={
            "sample_1": "run_1",
            "sample_2": "run_1",
            "sample_3": "run_2",
            "sample_4": "run_2",
        },
        condition_by_sample={
            "sample_1": "control",
            "sample_2": "treated",
            "sample_3": "control",
            "sample_4": "treated",
        },
        replicate_by_sample={
            "sample_1": "r1",
            "sample_2": "r1",
            "sample_3": "r2",
            "sample_4": "r2",
        },
        sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
    )


def _control_mapping():
    return ControlSiteEligibilityValidator().run(
        control_set=ControlSiteSet.from_site_keys(("site_a", "site_c")),
        site_keys=("site_a", "site_b", "site_c"),
        method="sps_ruv_style",
        min_eligible_controls=2,
        n_unwanted_factors=1,
    )


def _missingness_policy(
    *,
    missing_cells: tuple[tuple[str, str], ...],
) -> CorrectionMissingnessPolicy:
    return CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
        ),
        originally_missing_cells_tracked_by=OriginallyMissingCellTracking.OBSERVATION_MASK,
        correction_mask_policy=CorrectionMaskPolicy(),
        observation_mask=ObservationMask(
            feature_ids=tuple(str(value) for value in _site_index().tolist()),
            sample_ids=("sample_1", "sample_2", "sample_3", "sample_4"),
            originally_missing_cells=missing_cells,
        ),
    )


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [10.0, 5.0, 20.0],
            "sample_2": [10.0, 9.0, 20.0],
            "sample_3": [14.0, 8.0, 28.0],
            "sample_4": [14.0, 12.0, 28.0],
        },
        index=_site_index(),
    )


def _site_index() -> pd.Index:
    return pd.Index(["site_a", "site_b", "site_c"], name="site_key")


def _condition_effect(matrix: pd.DataFrame) -> pd.Series:
    treated_mean = matrix.loc[:, ["sample_2", "sample_4"]].mean(axis=1)
    control_mean = matrix.loc[:, ["sample_1", "sample_3"]].mean(axis=1)
    return treated_mean - control_mean
