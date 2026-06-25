from __future__ import annotations

from dataclasses import replace

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
from phospy.errors import PhosPyInputError
from phospy.science.batch_correction import (
    SPS_RUV_STYLE_EXECUTOR_ID,
    DeterministicSpsRuvStyleExecutor,
    SpsRuvStyleExecutorResult,
)
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteAnnotation,
    ControlSiteSet,
)
from phospy.validation.datasets.batch_correction import ResolvedBatchDesignMetadata
from phospy.validation.workflows.batch_correction import ControlSiteEligibilityValidator
from phospy.workflows.batch_correction import (
    BatchCorrectionPlanInterpreter,
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
)


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
    diagnostics_payload = result.diagnostics.to_payload()

    assert isinstance(result, SpsRuvStyleExecutorResult)
    assert result.diagnostics.status == "applied"
    assert result.diagnostics.executor_id == SPS_RUV_STYLE_EXECUTOR_ID
    assert result.diagnostics.estimated_unwanted_factors == 1
    assert result.diagnostics.singular_values[0] > 0.0
    assert diagnostics_payload["eligible_control_site_count"] == 2
    assert diagnostics_payload["rejected_control_site_count"] == 0
    assert diagnostics_payload["design_summary"]["number_of_batches"] == 2
    assert diagnostics_payload["design_summary"]["number_of_conditions"] == 2
    assert diagnostics_payload["batch_associated_variance"]["status"] == "computed"
    assert diagnostics_payload["missingness_imputation_summary"] == {
        "originally_missing_cell_count": 0,
        "withheld_cell_count": 0,
        "restored_missing_cell_count": 0,
        "temporary_imputation_applied": False,
        "temporary_imputation_allowed": True,
        "temporary_imputation_method": "row_median_temporary",
        "temporary_imputation_parameters": {"min_observed_values": 2},
        "output_policy": "no temporary imputation was needed",
    }
    assert result.warnings == ()
    assert result.withheld_rows == ()
    assert result.rejected_rows == ()
    assert result.rejected_cells == ()
    assert result.provenance_payload["executor_id"] == SPS_RUV_STYLE_EXECUTOR_ID
    assert result.provenance_payload["method"] == "sps_ruv_style"
    assert "corrected_matrix_fingerprint" in result.provenance_payload
    assert result.provenance_payload["diagnostics"] == diagnostics_payload
    assert result.corrected_preprocessing_output is not None
    assert result.corrected_preprocessing_output.batch_correction_report.status == (
        "applied"
    )


def test_sps_ruv_style_executor_defensively_rejects_non_estimable_factor_count() -> (
    None
):
    phospho = _phospho()
    phospho.loc["site_d", :] = [30.0, 30.0, 42.0, 42.0]

    with pytest.raises(
        PhosPyInputError,
        match="non-estimable n_unwanted_factors=2",
    ):
        _run_executor(
            phospho,
            control_site_set=ControlSiteSet.from_site_keys(
                ("site_a", "site_c", "site_d")
            ),
            n_unwanted_factors=2,
        )


def test_sps_ruv_style_executor_defensively_rejects_ruv_iii_style_plan() -> None:
    plan = BatchCorrectionPlanInterpreter().run(
        config=_config(),
        dataset_metadata=_metadata(),
        control_site_mapping=_control_mapping(),
        missingness_policy=_missingness_policy(missing_cells=()),
    )
    ruv_iii_plan = replace(plan, method="ruv_iii_style")

    with pytest.raises(
        PhosPyInputError,
        match="replicate-aware RUV-III numerical semantics are not implemented",
    ):
        DeterministicSpsRuvStyleExecutor().run(
            phospho=_phospho(),
            plan=ruv_iii_plan,
        )


def test_batch_correction_workflow_diagnostics_and_provenance_record_rejected_controls() -> (
    None
):
    control_site_set = ControlSiteSet(
        annotations=(
            ControlSiteAnnotation("site_a", control_status="control"),
            ControlSiteAnnotation(
                "site_b",
                control_status="excluded",
                exclusion_reason="unstable_reference_profile",
            ),
            ControlSiteAnnotation("site_c", control_status="control"),
        )
    )
    request = BatchCorrectionWorkflowRequest(
        phospho=_phospho(),
        config=_config(),
        sample_metadata=_sample_metadata(),
        control_site_set=control_site_set,
        missingness_policy=_missingness_policy(missing_cells=()),
    )

    result = BatchCorrectionWorkflow().run(request)

    executor_diagnostics = result.diagnostics["executor"]
    assert executor_diagnostics["eligible_control_site_count"] == 2
    assert executor_diagnostics["rejected_control_site_count"] == 1
    assert executor_diagnostics["rejected_control_sites"] == [
        {
            "site_key": "site_b",
            "scope": "row_eligibility",
            "control_status": "excluded",
            "row_position": 1,
            "annotation_indices": [1],
            "reasons": ["unstable_reference_profile"],
            "primary_reason": "unstable_reference_profile",
            "exclusion_reason": "unstable_reference_profile",
        }
    ]
    assert result.provenance.diagnostics["executor"] == executor_diagnostics
    assert result.provenance.preprocessing_stage_order == (
        "missing_data",
        "batch_correction",
        "downstream_workflows",
    )
    interpreter_plan = result.provenance.resolved_parameters["interpreter_plan"]
    assert interpreter_plan["executed_stage_order"] == list(
        result.provenance.preprocessing_stage_order
    )
    assert (
        interpreter_plan["requested_stage_order"]
        == "after_missing_data_before_downstream"
    )
    assert result.provenance.rejected_entities[0].entity_type == "site"
    assert result.provenance.rejected_entities[0].identifier == "site_b"
    assert result.provenance.rejected_entities[0].reason == (
        "unstable_reference_profile"
    )


def _run_executor(
    phospho: pd.DataFrame,
    *,
    missing_cells: tuple[tuple[str, str], ...] = (),
    control_site_set: ControlSiteSet | None = None,
    n_unwanted_factors: int = 1,
) -> SpsRuvStyleExecutorResult:
    plan = BatchCorrectionPlanInterpreter().run(
        config=_config(n_unwanted_factors=n_unwanted_factors),
        dataset_metadata=_metadata(),
        control_site_mapping=_control_mapping(
            site_keys=tuple(str(value) for value in phospho.index.tolist()),
            control_site_set=control_site_set,
            n_unwanted_factors=n_unwanted_factors,
        ),
        missingness_policy=_missingness_policy(
            missing_cells=missing_cells,
            feature_ids=tuple(str(value) for value in phospho.index.tolist()),
        ),
    )
    return DeterministicSpsRuvStyleExecutor().run(phospho=phospho, plan=plan)


def _config(*, n_unwanted_factors: int = 1) -> InternalBatchCorrectionRequest:
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
        n_unwanted_factors=n_unwanted_factors,
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


def _control_mapping(
    *,
    site_keys: tuple[str, ...] = ("site_a", "site_b", "site_c"),
    control_site_set: ControlSiteSet | None = None,
    n_unwanted_factors: int = 1,
):
    return ControlSiteEligibilityValidator().run(
        control_set=control_site_set
        if control_site_set is not None
        else ControlSiteSet.from_site_keys(("site_a", "site_c")),
        site_keys=site_keys,
        method="sps_ruv_style",
        min_eligible_controls=2,
        n_unwanted_factors=n_unwanted_factors,
    )


def _missingness_policy(
    *,
    missing_cells: tuple[tuple[str, str], ...],
    feature_ids: tuple[str, ...] | None = None,
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
            feature_ids=feature_ids
            if feature_ids is not None
            else tuple(str(value) for value in _site_index().tolist()),
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


def _sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("run_1", "run_1", "run_2", "run_2"),
            "condition": ("control", "treated", "control", "treated"),
            "replicate": ("r1", "r1", "r2", "r2"),
        },
        index=("sample_1", "sample_2", "sample_3", "sample_4"),
    )


def _condition_effect(matrix: pd.DataFrame) -> pd.Series:
    treated_mean = matrix.loc[:, ["sample_2", "sample_4"]].mean(axis=1)
    control_mean = matrix.loc[:, ["sample_1", "sample_3"]].mean(axis=1)
    return treated_mean - control_mean
