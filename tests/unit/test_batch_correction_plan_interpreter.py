from __future__ import annotations

import copy

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
from phospy.science.datasets.preprocessing.control_sites import ControlSiteSet
from phospy.validation.datasets.batch_correction import ResolvedBatchDesignMetadata
from phospy.validation.workflows.batch_correction import ControlSiteEligibilityValidator
from phospy.workflows.batch_correction import BatchCorrectionPlanInterpreter


def test_batch_correction_interpreter_resolves_controls_and_design() -> None:
    plan = BatchCorrectionPlanInterpreter().run(
        config=_config(),
        dataset_metadata=_metadata(),
        control_site_mapping=_control_mapping(),
        missingness_policy=_missingness_policy(),
    )

    assert plan.method == "sps_ruv_style"
    assert [row.site_key for row in plan.eligible_control_site_rows] == [
        "site_a",
        "site_c",
    ]
    assert [row.row_position for row in plan.eligible_control_site_rows] == [0, 2]
    assert plan.batch_terms == ("batch[run_2]",)
    assert plan.condition_terms_to_preserve == ("intercept", "condition[treated]")
    assert plan.resolved_design_matrix.loc["sample_2", "condition[treated]"] == 1.0
    assert plan.resolved_design_matrix.loc["sample_2", "batch[run_2]"] == 0.0
    assert plan.resolved_design_matrix.loc["sample_3", "batch[run_2]"] == 1.0


def test_batch_correction_interpreter_carries_replicates_and_missingness_policy() -> (
    None
):
    policy = _missingness_policy()

    plan = BatchCorrectionPlanInterpreter().run(
        config=_config(),
        dataset_metadata=_metadata(),
        control_site_mapping=_control_mapping(),
        missingness_policy=policy,
    )

    assert plan.replicate_structure.replicate_column == "replicate"
    assert plan.replicate_structure.replicate_labels == ("r1", "r1", "r2", "r2")
    assert plan.replicate_structure.replicate_groups == {
        "r1": ("sample_1", "sample_2"),
        "r2": ("sample_3", "sample_4"),
    }
    assert plan.observation_mask is policy.observation_mask
    assert plan.observation_mask.is_originally_missing("site_b", "sample_2")
    assert (
        plan.temporary_imputation_policy.method
        is TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY
    )


def test_batch_correction_interpreter_resolves_empty_observation_mask_when_complete() -> (
    None
):
    plan = BatchCorrectionPlanInterpreter().run(
        config=_config(
            missing_value_policy=(
                InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
            ),
            imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
        ),
        dataset_metadata=_metadata(),
        control_site_mapping=_control_mapping(),
        missingness_policy=CorrectionMissingnessPolicy(),
    )

    assert plan.observation_mask.feature_ids == ("site_a", "site_b", "site_c")
    assert plan.observation_mask.sample_ids == (
        "sample_1",
        "sample_2",
        "sample_3",
        "sample_4",
    )
    assert plan.observation_mask.originally_missing_cells == ()


def test_batch_correction_interpreter_resolves_stage_order_and_diagnostics() -> None:
    plan = BatchCorrectionPlanInterpreter().run(
        config=_config(
            diagnostics_enabled=True,
        ),
        dataset_metadata=_metadata(),
        control_site_mapping=_control_mapping(),
        missingness_policy=_missingness_policy(),
    )

    assert plan.stage_order == (
        "missing_data",
        "batch_correction",
        "downstream_workflows",
    )
    assert plan.stage_order_policy == ("after_missing_data_before_downstream")
    assert plan.to_payload()["executed_stage_order"] == list(plan.stage_order)
    assert plan.to_payload()["requested_stage_order"] == plan.stage_order_policy
    assert plan.diagnostic_requirements.diagnostics_enabled is True
    assert "diagnostic_tables" in plan.diagnostic_requirements.required_payloads


def test_batch_correction_interpreter_rejects_unexecutable_stage_order() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "stage_order='after_intensity_transform_before_missing_data' is "
            "unsupported.*supported stage order is missing_data -> "
            "batch_correction -> downstream_workflows.*provenance must match"
        ),
    ):
        BatchCorrectionPlanInterpreter().run(
            config=_config(
                stage_order=(
                    InternalBatchCorrectionStageOrder.AFTER_INTENSITY_TRANSFORM_BEFORE_MISSING_DATA
                ),
            ),
            dataset_metadata=_metadata(),
            control_site_mapping=_control_mapping(),
            missingness_policy=_missingness_policy(),
        )


def test_batch_correction_interpreter_provenance_seed_data_is_deterministic() -> None:
    kwargs = {
        "config": _config(),
        "dataset_metadata": _metadata(),
        "control_site_mapping": _control_mapping(),
        "missingness_policy": _missingness_policy(),
    }

    first = BatchCorrectionPlanInterpreter().run(**kwargs)
    second = BatchCorrectionPlanInterpreter().run(**copy.deepcopy(kwargs))

    assert first.provenance_seed_data == second.provenance_seed_data
    assert first.provenance_seed_data["method"] == "sps_ruv_style"
    assert first.provenance_seed_data["condition_columns"] == [
        "condition",
    ]
    assert first.provenance_seed_data["eligible_control_site_rows"] == [
        {
            "site_key": "site_a",
            "row_position": 0,
            "label": "caller_supplied_control",
            "weight": None,
            "group": None,
        },
        {
            "site_key": "site_c",
            "row_position": 2,
            "label": "caller_supplied_control",
            "weight": None,
            "group": None,
        },
    ]
    assert (
        first.provenance_seed_data["diagnostic_requirements"]
        == first.diagnostic_requirements.to_payload()
    )


def _config(
    *,
    missing_value_policy: InternalBatchCorrectionMissingValuePolicy = (
        InternalBatchCorrectionMissingValuePolicy.ALLOW_TEMPORARY_IMPUTATION
    ),
    imputation_policy: InternalBatchCorrectionImputationPolicy = (
        InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY
    ),
    stage_order: InternalBatchCorrectionStageOrder = (
        InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
    ),
    diagnostics_enabled: bool = False,
) -> InternalBatchCorrectionRequest:
    return InternalBatchCorrectionRequest(
        method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column="replicate",
        control_site_source=InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
        control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
        missing_value_policy=missing_value_policy,
        imputation_policy=imputation_policy,
        n_unwanted_factors=1,
        stage_order=stage_order,
        diagnostics_enabled=diagnostics_enabled,
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


def _missingness_policy() -> CorrectionMissingnessPolicy:
    return CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
        ),
        originally_missing_cells_tracked_by=OriginallyMissingCellTracking.OBSERVATION_MASK,
        correction_mask_policy=CorrectionMaskPolicy(),
        observation_mask=ObservationMask(
            feature_ids=("site_a", "site_b", "site_c"),
            sample_ids=("sample_1", "sample_2", "sample_3", "sample_4"),
            originally_missing_cells=(("site_b", "sample_2"),),
        ),
    )
