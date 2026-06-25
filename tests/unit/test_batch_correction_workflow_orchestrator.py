from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd
import pytest

from phospy.contracts.configs.preprocessing import (
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
)
from phospy.errors import PhosPyInputError
from phospy.provenance import fingerprint_matrix
from phospy.provenance.models import BatchCorrectionProvenance
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteMapping,
    ControlSiteSet,
)
from phospy.validation.datasets.batch_correction import ResolvedBatchDesignMetadata
from phospy.workflows.batch_correction import (
    BatchCorrectionDiagnosticRequirements,
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
    EligibleControlSiteRow,
    ReplicateStructure,
    ResolvedBatchCorrectionPlan,
)
from phospy.workflows.batch_correction.provenance import (
    BatchCorrectionProvenanceRecorder,
)


def test_batch_correction_workflow_orders_validators_before_interpreter_and_executor() -> (
    None
):
    order: list[str] = []
    request = _request()
    corrected = _matrix() + 1.0

    workflow = BatchCorrectionWorkflow(
        request_validator=_RequestValidator(order),
        design_validator=_DesignValidator(order),
        stage_order_validator=_StageOrderValidator(order),
        control_site_validator=_ControlSiteValidator(order),
        missingness_validator=_MissingnessValidator(order),
        interpreter=_Interpreter(order),
        executor=_Executor(order, corrected=corrected),
        provenance_recorder=_ProvenanceRecorder(order, corrected=corrected),
    )

    result = workflow.run(request)

    assert order == [
        "request_validator",
        "design_validator",
        "stage_order_validator",
        "control_site_validator",
        "missingness_validator",
        "interpreter",
        "executor",
        "provenance_recorder",
    ]
    assert result.corrected_matrix.equals(corrected)
    assert result.diagnostics["executor"]["executor_id"] == "fake_executor"


def test_batch_correction_workflow_invalid_request_fails_before_executor() -> None:
    order: list[str] = []
    workflow = BatchCorrectionWorkflow(
        request_validator=_RejectingRequestValidator(order),
        design_validator=_DesignValidator(order),
        stage_order_validator=_StageOrderValidator(order),
        control_site_validator=_ControlSiteValidator(order),
        missingness_validator=_MissingnessValidator(order),
        interpreter=_Interpreter(order),
        executor=_Executor(order, corrected=_matrix()),
        provenance_recorder=_ProvenanceRecorder(order, corrected=_matrix()),
    )

    with pytest.raises(PhosPyInputError, match="invalid test request"):
        workflow.run(_request())

    assert order == ["request_validator"]
    assert "executor" not in order


def test_batch_correction_workflow_rejects_unexecutable_stage_order_before_executor() -> (
    None
):
    order: list[str] = []
    workflow = BatchCorrectionWorkflow(
        request_validator=_RequestValidator(order),
        design_validator=_DesignValidator(order),
        control_site_validator=_ControlSiteValidator(order),
        missingness_validator=_MissingnessValidator(order),
        interpreter=_Interpreter(order),
        executor=_Executor(order, corrected=_matrix()),
        provenance_recorder=_ProvenanceRecorder(order, corrected=_matrix()),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "unsupported by the current dataset preprocessing pipeline.*"
            "supported stage order is missing_data -> batch_correction -> "
            "downstream_workflows.*provenance must match"
        ),
    ):
        workflow.run(
            _request(
                config=_config(
                    stage_order=(
                        InternalBatchCorrectionStageOrder.AFTER_INTENSITY_TRANSFORM_BEFORE_MISSING_DATA
                    )
                )
            )
        )

    assert order == ["request_validator", "design_validator"]
    assert "executor" not in order


def test_batch_correction_workflow_rejects_ruv_iii_before_interpreter_and_executor() -> (
    None
):
    order: list[str] = []
    workflow = BatchCorrectionWorkflow(
        design_validator=_DesignValidator(order),
        stage_order_validator=_StageOrderValidator(order),
        control_site_validator=_ControlSiteValidator(order),
        missingness_validator=_MissingnessValidator(order),
        interpreter=_Interpreter(order),
        executor=_Executor(order, corrected=_matrix()),
        provenance_recorder=_ProvenanceRecorder(order, corrected=_matrix()),
    )
    request = _request()
    forged_request = BatchCorrectionWorkflowRequest(
        phospho=request.phospho,
        config=_forged_ruv_iii_config(),
        sample_metadata=request.sample_metadata,
        control_site_set=request.control_site_set,
        missingness_policy=request.missingness_policy,
    )

    with pytest.raises(
        PhosPyInputError,
        match="replicate-aware RUV-III numerical semantics are not implemented",
    ):
        workflow.run(forged_request)

    assert order == []


def test_batch_correction_workflow_default_provenance_recorder_assembles_sources() -> (
    None
):
    order: list[str] = []
    request = _request()
    corrected = _matrix() + 2.0
    workflow = BatchCorrectionWorkflow(
        request_validator=_RequestValidator(order),
        design_validator=_DesignValidator(order),
        stage_order_validator=_StageOrderValidator(order),
        control_site_validator=_ControlSiteValidator(order),
        missingness_validator=_MissingnessValidator(order),
        interpreter=_Interpreter(order),
        executor=_Executor(order, corrected=corrected),
    )

    result = workflow.run(request)
    provenance = result.provenance

    assert provenance.requested_method == "sps_ruv_style"
    assert provenance.preprocessing_stage_order == (
        "missing_data",
        "batch_correction",
        "downstream_workflows",
    )
    assert provenance.selected_site_key_rows == ("site_a", "site_b")
    assert (
        provenance.input_matrix_fingerprint.exact_hash_value
        == fingerprint_matrix(
            request.phospho,
            name="batch_correction.workflow.input",
        ).exact_hash_value
    )
    assert (
        provenance.output_matrix_fingerprint is not None
        and provenance.output_matrix_fingerprint.exact_hash_value
        == fingerprint_matrix(
            corrected,
            name="batch_correction.workflow.corrected",
        ).exact_hash_value
    )
    assert provenance.resolved_parameters["config"]["batch_column"] == "batch"
    assert (
        provenance.resolved_parameters["interpreter_seed_data"]["seed_source"]
        == "fake_interpreter"
    )
    assert provenance.resolved_parameters["executor"]["executor_id"] == "fake_executor"
    assert provenance.diagnostics["executor"]["executor_id"] == "fake_executor"
    assert provenance.warnings == ("executor warning",)
    assert provenance.phospy_version
    assert provenance.phospy_version != "unknown"
    assert provenance.python_version
    assert provenance.python_version != "unknown"
    assert {"numpy", "pandas"}.issubset(set(provenance.dependency_versions))


def test_batch_correction_provenance_recorder_populates_environment_fields() -> None:
    corrected = _matrix() + 3.0

    provenance = BatchCorrectionProvenanceRecorder().run(
        request=_request(),
        dataset_metadata=_metadata(),
        control_site_mapping=_control_mapping(),
        missingness_policy=_missingness_policy(),
        plan=_plan(),
        executor_result=_ExecutorResult(corrected_matrix=corrected),
    )

    assert provenance.phospy_version
    assert provenance.phospy_version != "unknown"
    assert provenance.python_version
    assert provenance.python_version != "unknown"
    assert provenance.dependency_versions
    assert {"numpy", "pandas", "scipy", "scikit-learn"}.issubset(
        set(provenance.dependency_versions)
    )


def test_batch_correction_workflow_combines_upstream_mask_with_executor_mask() -> None:
    order: list[str] = []
    request = _request(upstream_observation_mask=_upstream_observation_mask())
    recorder = _CapturingProvenanceRecorder(order, corrected=_matrix())
    workflow = BatchCorrectionWorkflow(
        request_validator=_RequestValidator(order),
        design_validator=_DesignValidator(order),
        stage_order_validator=_StageOrderValidator(order),
        control_site_validator=_ControlSiteValidator(order),
        missingness_validator=_MissingnessValidator(order),
        interpreter=_Interpreter(order),
        executor=_Executor(order, corrected=_matrix()),
        provenance_recorder=recorder,
    )

    workflow.run(request)

    assert recorder.output_observation_mask is not None
    assert bool(recorder.output_observation_mask.loc["site_b", "sample_2"]) is False


class _RequestValidator:
    def __init__(self, order: list[str]) -> None:
        self._order = order

    def run(self, request: object) -> BatchCorrectionWorkflowRequest:
        self._order.append("request_validator")
        assert isinstance(request, BatchCorrectionWorkflowRequest)
        return request


class _RejectingRequestValidator:
    def __init__(self, order: list[str]) -> None:
        self._order = order

    def run(self, request: object) -> BatchCorrectionWorkflowRequest:
        self._order.append("request_validator")
        raise PhosPyInputError("invalid test request")


class _DesignValidator:
    def __init__(self, order: list[str]) -> None:
        self._order = order

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
    ) -> ResolvedBatchDesignMetadata:
        self._order.append("design_validator")
        return _metadata()


class _StageOrderValidator:
    def __init__(self, order: list[str]) -> None:
        self._order = order

    def run(self, *, config: InternalBatchCorrectionRequest) -> None:
        self._order.append("stage_order_validator")


class _ControlSiteValidator:
    def __init__(self, order: list[str]) -> None:
        self._order = order

    def run(self, *, request: BatchCorrectionWorkflowRequest) -> ControlSiteMapping:
        self._order.append("control_site_validator")
        return _control_mapping()


class _MissingnessValidator:
    def __init__(self, order: list[str]) -> None:
        self._order = order

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
    ) -> CorrectionMissingnessPolicy:
        self._order.append("missingness_validator")
        return _missingness_policy()


class _Interpreter:
    def __init__(self, order: list[str]) -> None:
        self._order = order

    def run(self, **kwargs: Any) -> ResolvedBatchCorrectionPlan:
        self._order.append("interpreter")
        return _plan()


class _Executor:
    def __init__(self, order: list[str], *, corrected: pd.DataFrame) -> None:
        self._order = order
        self._corrected = corrected

    def run(self, **kwargs: Any) -> _ExecutorResult:
        self._order.append("executor")
        return _ExecutorResult(corrected_matrix=self._corrected)


class _ProvenanceRecorder:
    def __init__(self, order: list[str], *, corrected: pd.DataFrame) -> None:
        self._order = order
        self._corrected = corrected

    def run(self, **kwargs: Any) -> BatchCorrectionProvenance:
        self._order.append("provenance_recorder")
        return _provenance(corrected=self._corrected)


class _CapturingProvenanceRecorder(_ProvenanceRecorder):
    output_observation_mask: pd.DataFrame | None = None

    def run(self, **kwargs: Any) -> BatchCorrectionProvenance:
        executor_result = kwargs["executor_result"]
        self.output_observation_mask = executor_result.output_observation_mask
        return super().run(**kwargs)


@dataclass(frozen=True, slots=True)
class _Diagnostics:
    def to_payload(self) -> dict[str, object]:
        return {"executor_id": "fake_executor", "status": "applied"}


@dataclass(frozen=True, slots=True)
class _ExecutorResult:
    corrected_matrix: pd.DataFrame
    diagnostics: _Diagnostics = _Diagnostics()
    warnings: tuple[str, ...] = ("executor warning",)
    rejected_rows: tuple[str, ...] = ()
    rejected_cells: tuple[tuple[str, str], ...] = ()
    withheld_rows: tuple[str, ...] = ()
    withheld_cells: tuple[tuple[str, str], ...] = ()
    provenance_payload: dict[str, object] | None = None

    @property
    def output_observation_mask(self) -> pd.DataFrame:
        return pd.DataFrame(
            True,
            index=self.corrected_matrix.index.copy(),
            columns=self.corrected_matrix.columns.copy(),
        )

    def __post_init__(self) -> None:
        if self.provenance_payload is None:
            object.__setattr__(
                self,
                "provenance_payload",
                {"executor_id": "fake_executor", "payload_source": "fake"},
            )


def _request(
    *,
    config: InternalBatchCorrectionRequest | None = None,
    upstream_observation_mask: pd.DataFrame | None = None,
) -> BatchCorrectionWorkflowRequest:
    return BatchCorrectionWorkflowRequest(
        phospho=_matrix(),
        config=_config() if config is None else config,
        sample_metadata=pd.DataFrame(
            {
                "batch": ("run_1", "run_1", "run_2", "run_2"),
                "condition": ("control", "treated", "control", "treated"),
            },
            index=("sample_1", "sample_2", "sample_3", "sample_4"),
        ),
        control_site_set=ControlSiteSet.from_site_keys(("site_a", "site_b")),
        missingness_policy=_missingness_policy(),
        upstream_observation_mask=upstream_observation_mask,
    )


def _config(
    *,
    stage_order: InternalBatchCorrectionStageOrder = (
        InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
    ),
) -> InternalBatchCorrectionRequest:
    return InternalBatchCorrectionRequest(
        method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column=None,
        control_site_source=InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
        control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
        missing_value_policy=InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING,
        imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
        n_unwanted_factors=1,
        stage_order=stage_order,
        diagnostics_enabled=True,
    )


def _forged_ruv_iii_config() -> InternalBatchCorrectionRequest:
    config = object.__new__(InternalBatchCorrectionRequest)
    object.__setattr__(config, "method", InternalBatchCorrectionMethod.RUV_III_STYLE)
    object.__setattr__(config, "batch_column", "batch")
    object.__setattr__(config, "condition_columns", ("condition",))
    object.__setattr__(config, "replicate_column", "replicate")
    object.__setattr__(
        config,
        "control_site_source",
        InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
    )
    object.__setattr__(
        config,
        "control_site_mode",
        InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
    )
    object.__setattr__(
        config,
        "missing_value_policy",
        InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING,
    )
    object.__setattr__(
        config,
        "imputation_policy",
        InternalBatchCorrectionImputationPolicy.NONE,
    )
    object.__setattr__(config, "n_unwanted_factors", 1)
    object.__setattr__(
        config,
        "stage_order",
        InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM,
    )
    object.__setattr__(config, "diagnostics_enabled", True)
    return config


def _matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [
            (10.0, 11.0, 13.0, 14.0),
            (20.0, 22.0, 25.0, 27.0),
        ],
        index=pd.Index(("site_a", "site_b"), name="site_key"),
        columns=("sample_1", "sample_2", "sample_3", "sample_4"),
    )


def _upstream_observation_mask() -> pd.DataFrame:
    mask = pd.DataFrame(
        True,
        index=_matrix().index.copy(),
        columns=_matrix().columns.copy(),
    )
    mask.loc["site_b", "sample_2"] = False
    return mask


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
        sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
    )


def _control_mapping() -> ControlSiteMapping:
    return ControlSiteSet.from_site_keys(("site_a", "site_b")).map_to_site_keys(
        ("site_a", "site_b")
    )


def _missingness_policy() -> CorrectionMissingnessPolicy:
    return CorrectionMissingnessPolicy(
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.OBSERVATION_MASK
        ),
        observation_mask=ObservationMask(
            feature_ids=("site_a", "site_b"),
            sample_ids=("sample_1", "sample_2", "sample_3", "sample_4"),
            originally_missing_cells=(),
        ),
    )


def _plan() -> ResolvedBatchCorrectionPlan:
    design = pd.DataFrame(
        [
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, 1.0),
        ],
        index=pd.Index(("sample_1", "sample_2", "sample_3", "sample_4"), name="sample"),
        columns=pd.Index(
            ("intercept", "condition[treated]", "batch[run_2]"),
            name="term",
        ),
    )
    policy = _missingness_policy()
    return ResolvedBatchCorrectionPlan(
        method="sps_ruv_style",
        resolved_design_matrix=design,
        batch_terms=("batch[run_2]",),
        condition_terms_to_preserve=("intercept", "condition[treated]"),
        replicate_structure=ReplicateStructure(
            replicate_column=None,
            replicate_by_sample=None,
            replicate_labels=None,
            replicate_groups={},
        ),
        eligible_control_site_rows=(
            EligibleControlSiteRow(site_key="site_a", row_position=0),
            EligibleControlSiteRow(site_key="site_b", row_position=1),
        ),
        observation_mask=policy.observation_mask,  # type: ignore[arg-type]
        temporary_imputation_policy=policy.temporary_imputation,
        n_unwanted_factors=1,
        stage_order=("missing_data", "batch_correction", "downstream_workflows"),
        stage_order_policy="after_missing_data_before_downstream",
        diagnostic_requirements=BatchCorrectionDiagnosticRequirements(
            diagnostics_enabled=True,
            required_payloads=("resolved_design_matrix",),
        ),
        provenance_seed_data={"seed_source": "fake_interpreter"},
    )


def _provenance(*, corrected: pd.DataFrame) -> BatchCorrectionProvenance:
    return BatchCorrectionProvenance(
        requested_method="sps_ruv_style",
        resolved_parameters={},
        preprocessing_stage_order=("batch_correction",),
        control_site_source={},
        selected_site_key_rows=("site_a", "site_b"),
        batch_metadata={},
        replicate_metadata=None,
        design_metadata={},
        missing_value_policy={},
        observation_masks=(),
        input_matrix_fingerprint=fingerprint_matrix(_matrix(), name="fake.input"),
        output_matrix_fingerprint=fingerprint_matrix(corrected, name="fake.output"),
    )
