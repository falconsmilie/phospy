from __future__ import annotations

from dataclasses import replace

import numpy as np
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
from phospy.contracts.requests import DatasetBuildRequest
from phospy.errors import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.provenance import BatchCorrectionProvenance, fingerprint_matrix
from phospy.science.batch_correction import (
    SPS_RUV_STYLE_ALGORITHM_DESCRIPTION,
    SPS_RUV_STYLE_BATCH_TERM_ROLE,
    SPS_RUV_STYLE_EXECUTOR_ID,
    SPS_RUV_STYLE_REPLICATE_METADATA_ROLE,
    DeterministicSpsRuvStyleExecutor,
    SpsRuvStyleExecutorResult,
)
from phospy.science.datasets.builders.public import AnalysisReadyDatasetBuilder
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteAnnotation,
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.science.references.models import Organism
from phospy.validation.datasets.batch_correction import ResolvedBatchDesignMetadata
from phospy.validation.workflows.batch_correction import ControlSiteEligibilityValidator
from phospy.workflows.batch_correction import (
    BatchCorrectionPlanInterpreter,
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
    ResolvedBatchCorrectionPlan,
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


def test_sps_ruv_style_executor_batch_labels_affect_diagnostics_not_correction() -> (
    None
):
    phospho = _phospho()
    two_batch_metadata = _metadata()
    relabeled_batch_metadata = _metadata(
        batch_by_sample={
            "sample_1": "run_1",
            "sample_2": "run_2",
            "sample_3": "run_2",
            "sample_4": "run_1",
        }
    )

    first = _run_executor(phospho, dataset_metadata=two_batch_metadata)
    second = _run_executor(phospho, dataset_metadata=relabeled_batch_metadata)

    pdt.assert_frame_equal(
        first.corrected_matrix,
        second.corrected_matrix,
        check_exact=False,
        atol=1e-10,
        rtol=0.0,
    )
    pdt.assert_frame_equal(
        first.estimated_unwanted_factors,
        second.estimated_unwanted_factors,
        check_exact=False,
        atol=1e-10,
        rtol=0.0,
    )
    assert first.diagnostics.design_summary["number_of_batches"] == 2
    assert second.diagnostics.design_summary["number_of_batches"] == 2
    assert first.diagnostics.design_summary["design_matrix_shape"] == [4, 3]
    assert second.diagnostics.design_summary["design_matrix_shape"] == [4, 3]
    first_before = first.diagnostics.batch_associated_variance["before"]
    second_before = second.diagnostics.batch_associated_variance["before"]
    assert isinstance(first_before, dict)
    assert isinstance(second_before, dict)
    assert first_before["mean_r_squared"] != second_before["mean_r_squared"]


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
    assert result.corrected_preprocessing_output is None


def test_sps_ruv_style_executor_rejects_positive_infinity_in_observed_cell() -> None:
    phospho = _phospho()
    phospho.loc["site_b", "sample_2"] = np.inf

    with pytest.raises(
        PhosPyInputError,
        match=(
            "non-finite observed matrix value.*feature_id='site_b'.*"
            "sample_id='sample_2'"
        ),
    ):
        _run_executor(phospho)


def test_sps_ruv_style_executor_rejects_negative_infinity_in_observed_cell() -> None:
    phospho = _phospho()
    phospho.loc["site_b", "sample_2"] = -np.inf

    with pytest.raises(
        PhosPyInputError,
        match=(
            "non-finite observed matrix value.*feature_id='site_b'.*"
            "sample_id='sample_2'"
        ),
    ):
        _run_executor(phospho)


def test_sps_ruv_style_executor_rejects_duplicate_row_index_before_correction() -> None:
    phospho = _phospho()
    phospho.index = pd.Index(["site_a", "site_b", "site_a"], name="site_key")

    with pytest.raises(
        PhosPyInputError,
        match=r"feature/site labels.*unique.*'site_a'",
    ):
        _run_direct_executor_with_reference_plan(phospho)


def test_sps_ruv_style_executor_rejects_duplicate_column_label_before_correction() -> (
    None
):
    phospho = _phospho()
    phospho.columns = pd.Index(["sample_1", "sample_2", "sample_2", "sample_4"])

    with pytest.raises(
        PhosPyInputError,
        match=r"sample/column labels.*unique.*'sample_2'",
    ):
        _run_direct_executor_with_reference_plan(phospho)


def test_sps_ruv_style_executor_rejects_blank_row_index_before_correction() -> None:
    phospho = _phospho()
    phospho.index = pd.Index(["site_a", "   ", "site_c"], name="site_key")

    with pytest.raises(
        PhosPyInputError,
        match=r"feature/site labels.*nonblank.*blank label positions: 1",
    ):
        _run_direct_executor_with_reference_plan(phospho)


def test_sps_ruv_style_executor_rejects_blank_column_label_before_correction() -> None:
    phospho = _phospho()
    phospho.columns = pd.Index(["sample_1", "", "sample_3", "sample_4"])

    with pytest.raises(
        PhosPyInputError,
        match=r"sample/column labels.*nonblank.*blank label positions: 1",
    ):
        _run_direct_executor_with_reference_plan(phospho)


def test_sps_ruv_style_executor_returns_diagnostics_warnings_and_provenance() -> None:
    result = _run_executor(_phospho())
    diagnostics_payload = result.diagnostics.to_payload()

    assert isinstance(result, SpsRuvStyleExecutorResult)
    assert result.diagnostics.status == "applied"
    assert result.diagnostics.executor_id == SPS_RUV_STYLE_EXECUTOR_ID
    assert result.diagnostics.estimated_unwanted_factors == 1
    assert result.diagnostics.singular_values[0] > 0.0
    assert (
        diagnostics_payload["algorithm_description"]
        == SPS_RUV_STYLE_ALGORITHM_DESCRIPTION
    )
    assert diagnostics_payload["term_roles"]["batch_terms"] == (
        SPS_RUV_STYLE_BATCH_TERM_ROLE
    )
    assert diagnostics_payload["term_roles"]["replicate_metadata"] == (
        SPS_RUV_STYLE_REPLICATE_METADATA_ROLE
    )
    assert diagnostics_payload["eligible_control_site_count"] == 2
    assert diagnostics_payload["rejected_control_site_count"] == 0
    assert diagnostics_payload["design_summary"]["number_of_batches"] == 2
    assert diagnostics_payload["design_summary"]["number_of_conditions"] == 2
    assert diagnostics_payload["batch_associated_variance"]["status"] == "computed"
    assert diagnostics_payload["batch_associated_variance"]["batch_term_role"] == (
        SPS_RUV_STYLE_BATCH_TERM_ROLE
    )
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
    assert result.provenance_payload["resolved_plan"]["method"] == "sps_ruv_style"
    assert (
        result.provenance_payload["provenance_seed_data"]["method"] == "sps_ruv_style"
    )
    assert result.provenance_payload["algorithm_description"] == (
        SPS_RUV_STYLE_ALGORITHM_DESCRIPTION
    )
    assert result.provenance_payload["term_roles"]["batch_terms"] == (
        SPS_RUV_STYLE_BATCH_TERM_ROLE
    )
    assert result.provenance_payload["term_roles"]["replicate_metadata"] == (
        SPS_RUV_STYLE_REPLICATE_METADATA_ROLE
    )
    assert "corrected_matrix_fingerprint" in result.provenance_payload
    assert result.provenance_payload["diagnostics"] == diagnostics_payload
    assert result.corrected_preprocessing_output is not None
    report = result.corrected_preprocessing_output.batch_correction_report
    assert SPS_RUV_STYLE_ALGORITHM_DESCRIPTION in report.limitations
    assert SPS_RUV_STYLE_BATCH_TERM_ROLE in report.limitations
    assert result.corrected_preprocessing_output.batch_correction_report.status == (
        "applied"
    )


def test_sps_ruv_style_executor_rejects_internal_method_alias_plan() -> None:
    phospho = _phospho()
    plan = replace(_resolved_plan(phospho), method="control_site_ruv_style")

    with pytest.raises(
        PhosPyInputError,
        match=(
            "SPS/RUV-style executor requires a resolved SPS/RUV-style plan; "
            "got method='control_site_ruv_style'"
        ),
    ):
        DeterministicSpsRuvStyleExecutor().run(phospho=phospho, plan=plan)


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
                ("site_a", "site_c", "site_d"),
                source_metadata=_control_source_metadata(),
            ),
            n_unwanted_factors=2,
        )


def test_sps_ruv_style_executor_defensively_rejects_one_control_one_factor() -> None:
    plan = _resolved_plan(_phospho())
    underpowered_plan = replace(
        plan,
        eligible_control_site_rows=plan.eligible_control_site_rows[:1],
    )

    with pytest.raises(
        PhosPyInputError,
        match="n_unwanted_factors=1.*available eligible control count=1",
    ):
        DeterministicSpsRuvStyleExecutor().run(
            phospho=_phospho(),
            plan=underpowered_plan,
        )


def test_sps_ruv_style_executor_defensively_rejects_two_controls_two_factors() -> None:
    underpowered_plan = replace(_resolved_plan(_phospho()), n_unwanted_factors=2)

    with pytest.raises(
        PhosPyInputError,
        match="n_unwanted_factors=2.*available eligible control count=2",
    ):
        DeterministicSpsRuvStyleExecutor().run(
            phospho=_phospho(),
            plan=underpowered_plan,
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
        ),
        source_metadata=_control_source_metadata(),
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
    assert executor_diagnostics["algorithm_description"] == (
        SPS_RUV_STYLE_ALGORITHM_DESCRIPTION
    )
    assert executor_diagnostics["term_roles"]["batch_terms"] == (
        SPS_RUV_STYLE_BATCH_TERM_ROLE
    )
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
    provenance_executor = result.provenance.resolved_parameters["executor"]
    assert provenance_executor["algorithm_description"] == (
        SPS_RUV_STYLE_ALGORITHM_DESCRIPTION
    )
    assert provenance_executor["term_roles"]["batch_terms"] == (
        SPS_RUV_STYLE_BATCH_TERM_ROLE
    )
    assert result.provenance.preprocessing_stage_order == (
        "missing_data",
        "batch_correction",
        "downstream_workflows",
    )
    interpreter_plan = result.provenance.resolved_parameters["interpreter_plan"]
    assert interpreter_plan["executed_stage_order"] == tuple(
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


def test_batch_correction_workflow_rejects_actual_missing_cells_before_execution() -> (
    None
):
    phospho = _phospho()
    phospho.loc["site_b", "sample_2"] = pd.NA

    with pytest.raises(
        PhosPyInputError,
        match=(
            "native SPS/RUV-style correction.*actual missing values.*"
            "temporary imputation.*restored missing values.*cannot produce "
            "analysis-ready corrected output.*run missing-data preprocessing "
            "first or provide a complete upstream-imputed matrix with an "
            "observation mask"
        ),
    ):
        BatchCorrectionWorkflow().run(
            BatchCorrectionWorkflowRequest(
                phospho=phospho,
                config=_config(),
                sample_metadata=_sample_metadata(),
                control_site_set=ControlSiteSet.from_site_keys(
                    ("site_a", "site_c"),
                    source_metadata=_control_source_metadata(),
                ),
                missingness_policy=_missingness_policy(
                    missing_cells=(("site_b", "sample_2"),)
                ),
            )
        )


def test_batch_correction_workflow_complete_matrix_returns_applied_output() -> None:
    result = BatchCorrectionWorkflow().run(
        BatchCorrectionWorkflowRequest(
            phospho=_phospho(),
            config=_config(),
            sample_metadata=_sample_metadata(),
            control_site_set=ControlSiteSet.from_site_keys(
                ("site_a", "site_c"),
                source_metadata=_control_source_metadata(),
            ),
            missingness_policy=_missingness_policy(missing_cells=()),
        )
    )

    assert isinstance(
        result.corrected_preprocessing_output, CorrectedPreprocessingOutput
    )
    assert not bool(result.corrected_matrix.isna().to_numpy().any())
    pdt.assert_frame_equal(
        result.corrected_matrix,
        result.corrected_preprocessing_output.corrected_matrix,
    )


def test_dataset_builder_still_rejects_incomplete_corrected_output() -> None:
    phospho = _builder_phospho()
    corrected = phospho.copy(deep=True)
    corrected.loc[_builder_site_index()[1], "sample_2"] = pd.NA

    with pytest.raises(
        DatasetValidationError,
        match="dataset.phospho must not contain missing values",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_builder_site_metadata(),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_forged_incomplete_corrected_output(
                    corrected,
                    input_matrix=phospho,
                ),
            )
        )


def _run_executor(
    phospho: pd.DataFrame,
    *,
    missing_cells: tuple[tuple[str, str], ...] = (),
    control_site_set: ControlSiteSet | None = None,
    n_unwanted_factors: int = 1,
    dataset_metadata: ResolvedBatchDesignMetadata | None = None,
) -> SpsRuvStyleExecutorResult:
    plan = _resolved_plan(
        phospho,
        missing_cells=missing_cells,
        control_site_set=control_site_set,
        n_unwanted_factors=n_unwanted_factors,
        dataset_metadata=dataset_metadata,
    )
    return DeterministicSpsRuvStyleExecutor().run(phospho=phospho, plan=plan)


def _run_direct_executor_with_reference_plan(
    phospho: pd.DataFrame,
) -> SpsRuvStyleExecutorResult:
    plan = _resolved_plan(_phospho())
    return DeterministicSpsRuvStyleExecutor().run(phospho=phospho, plan=plan)


def _resolved_plan(
    phospho: pd.DataFrame,
    *,
    missing_cells: tuple[tuple[str, str], ...] = (),
    control_site_set: ControlSiteSet | None = None,
    n_unwanted_factors: int = 1,
    dataset_metadata: ResolvedBatchDesignMetadata | None = None,
) -> ResolvedBatchCorrectionPlan:
    return BatchCorrectionPlanInterpreter().run(
        config=_config(n_unwanted_factors=n_unwanted_factors),
        dataset_metadata=_metadata() if dataset_metadata is None else dataset_metadata,
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


def _metadata(
    *,
    batch_by_sample: dict[str, str] | None = None,
) -> ResolvedBatchDesignMetadata:
    return ResolvedBatchDesignMetadata(
        batch_by_sample=batch_by_sample
        if batch_by_sample is not None
        else {
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
            "sample_2": "r2",
            "sample_3": "r2",
            "sample_4": "r1",
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
        else ControlSiteSet.from_site_keys(
            ("site_a", "site_c"),
            source_metadata=_control_source_metadata(),
        ),
        site_keys=site_keys,
        method="sps_ruv_style",
        min_eligible_controls=2,
        n_unwanted_factors=n_unwanted_factors,
    )


def _control_source_metadata() -> ControlSiteSourceMetadata:
    return ControlSiteSourceMetadata(
        organism="rat",
        identifier_namespace="site_key",
        source_name="manual-curated-controls",
        source_version="manual-v1",
        license="caller local use",
        redistribution="not redistributed",
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


def _builder_site_index() -> pd.Index:
    return pd.Index(
        [
            (
                "phospy:v1|organism=rat|protein_namespace=protein_id|"
                "protein_identifier=MAPK14|residue=Y|position=182"
            ),
            (
                "phospy:v1|organism=rat|protein_namespace=protein_id|"
                "protein_identifier=AKT1|residue=T|position=308"
            ),
            (
                "phospy:v1|organism=rat|protein_namespace=protein_id|"
                "protein_identifier=GSK3B|residue=S|position=9"
            ),
        ],
        name="site_key",
    )


def _builder_phospho() -> pd.DataFrame:
    phospho = _phospho()
    phospho.index = _builder_site_index()
    return phospho


def _sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("run_1", "run_1", "run_2", "run_2"),
            "condition": ("control", "treated", "control", "treated"),
            "replicate": ("r1", "r2", "r2", "r1"),
        },
        index=("sample_1", "sample_2", "sample_3", "sample_4"),
    )


def _builder_site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ("MAPK14", "AKT1", "GSK3B"),
            "protein_id": ("MAPK14", "AKT1", "GSK3B"),
            "site": ("Y182", "T308", "S9"),
            "site_sequence": (
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "T" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ),
            "localisation_confidence": (0.95, 0.92, 0.9),
        },
        index=_builder_site_index(),
    )


def _forged_incomplete_corrected_output(
    corrected: pd.DataFrame,
    *,
    input_matrix: pd.DataFrame,
) -> CorrectedPreprocessingOutput:
    observation_mask = pd.DataFrame(
        True,
        index=corrected.index.copy(),
        columns=corrected.columns.copy(),
    )
    output = object.__new__(CorrectedPreprocessingOutput)
    object.__setattr__(output, "corrected_matrix", corrected)
    object.__setattr__(
        output,
        "batch_correction_report",
        BatchCorrectionReport(
            status="applied",
            policy=BatchCorrectionPolicy(
                method="sps_ruv_style",
                batch_column="batch",
                condition_column="condition",
                condition_columns=("condition",),
                preserve_condition_effects=True,
            ),
            diagnostics=BatchCorrectionDiagnostics(
                matrix_shape_before=(int(corrected.shape[0]), int(corrected.shape[1])),
                matrix_shape_after=(int(corrected.shape[0]), int(corrected.shape[1])),
            ),
        ),
    )
    object.__setattr__(
        output,
        "diagnostics",
        {"test": "forged incomplete corrected output"},
    )
    object.__setattr__(output, "output_observation_mask", observation_mask)
    object.__setattr__(output, "corrected_cell_status", None)
    object.__setattr__(
        output,
        "provenance",
        BatchCorrectionProvenance(
            requested_method="sps_ruv_style",
            resolved_parameters={
                "method": "sps_ruv_style",
                "n_unwanted_factors": 1,
                "source": "forged_unit_test",
            },
            preprocessing_stage_order=(
                "missing_data",
                "batch_correction",
                "downstream_workflows",
            ),
            control_site_source={
                "source_type": "caller_supplied",
                "organism": "rat",
                "identifier_namespace": "site_key",
                "source_name": "manual-curated-controls",
                "source_version": "manual-v1",
                "license": "caller local use",
                "redistribution": "not redistributed",
            },
            selected_site_key_rows=tuple(str(row) for row in corrected.index[:2]),
            batch_metadata={
                "column": "batch",
                "levels": ["run_1", "run_2"],
                "sample_order": list(corrected.columns.astype(str)),
            },
            replicate_metadata=None,
            design_metadata={
                "condition_columns": ["condition"],
                "preserve_condition_effects": True,
            },
            missing_value_policy={
                "policy": "reject_missing",
                "imputation_policy": "none",
            },
            observation_masks=(
                fingerprint_matrix(
                    observation_mask.astype("int8"),
                    name="batch_correction.native.observation_mask",
                ),
            ),
            input_matrix_fingerprint=fingerprint_matrix(
                input_matrix,
                name="batch_correction.native.input",
            ),
            output_matrix_fingerprint=fingerprint_matrix(
                corrected,
                name="batch_correction.native.corrected",
            ),
            diagnostics={"executor": {"status": "applied", "method": "sps_ruv_style"}},
            warnings=(),
            phospy_version="test",
            python_version="3.test",
            dependency_versions={
                "numpy": "test-numpy",
                "pandas": "test-pandas",
                "scipy": "test-scipy",
                "scikit-learn": "test-scikit-learn",
            },
        ),
    )
    object.__setattr__(
        output,
        "stage_order",
        ("missing_data", "batch_correction", "downstream_workflows"),
    )
    object.__setattr__(output, "consumed_by_downstream", False)
    return output


def _condition_effect(matrix: pd.DataFrame) -> pd.Series:
    treated_mean = matrix.loc[:, ["sample_2", "sample_4"]].mean(axis=1)
    control_mean = matrix.loc[:, ["sample_1", "sample_3"]].mean(axis=1)
    return treated_mean - control_mean
