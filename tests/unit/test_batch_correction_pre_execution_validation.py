from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

import pandas as pd
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
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.workflows.batch_correction import (
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
)
from phospy.workflows.batch_correction.contracts import (
    BatchCorrectionExecutorContract,
    BatchCorrectionInterpreterContract,
)


def test_workflow_rejects_knn_temporary_imputation_before_interpreter_and_executor() -> (
    None
):
    executor = _SpyExecutor()
    interpreter = _SpyInterpreter()
    request = _request(
        phospho=_phospho_with_missing(),
        config=_config(
            imputation_policy=InternalBatchCorrectionImputationPolicy.KNN_TEMPORARY
        ),
        missingness_policy=_missingness_policy(
            method=TemporaryImputationMethod.KNN_TEMPORARY,
            parameters={
                "k": 3,
                "distance": "nan_euclidean",
                "max_missing_fraction_per_row": 0.5,
            },
            missing_cells=(("site_b", "sample_2"),),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="unsupported temporary imputation.*KNN temporary imputation is not implemented",
    ):
        BatchCorrectionWorkflow(
            interpreter=cast(BatchCorrectionInterpreterContract, interpreter),
            executor=cast(BatchCorrectionExecutorContract, executor),
        ).run(request)

    assert interpreter.call_count == 0
    assert executor.call_count == 0


def test_workflow_rejects_minprob_temporary_imputation_before_interpreter_and_executor() -> (
    None
):
    executor = _SpyExecutor()
    interpreter = _SpyInterpreter()
    request = _request(
        phospho=_phospho_with_missing(),
        config=_config(
            imputation_policy=InternalBatchCorrectionImputationPolicy.MINPROB_TEMPORARY
        ),
        missingness_policy=_missingness_policy(
            method=TemporaryImputationMethod.MINPROB_TEMPORARY,
            parameters={"q": 0.01, "width": 0.3},
            random_seed=101,
            missing_cells=(("site_b", "sample_2"),),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="unsupported temporary imputation.*MinProb temporary imputation is not implemented",
    ):
        BatchCorrectionWorkflow(
            interpreter=cast(BatchCorrectionInterpreterContract, interpreter),
            executor=cast(BatchCorrectionExecutorContract, executor),
        ).run(request)

    assert interpreter.call_count == 0
    assert executor.call_count == 0


def test_workflow_rejects_random_temporary_imputation_without_seed_before_executor() -> (
    None
):
    executor = _SpyExecutor()
    request = _request(
        phospho=_phospho_with_missing(),
        config=_config(
            imputation_policy=InternalBatchCorrectionImputationPolicy.MINPROB_TEMPORARY
        ),
        missingness_policy=_missingness_policy(
            temporary_imputation=_forged_minprob_policy_without_seed(),
            missing_cells=(("site_b", "sample_2"),),
        ),
    )

    with pytest.raises(PhosPyInputError, match="deterministic seed"):
        BatchCorrectionWorkflow(
            executor=cast(BatchCorrectionExecutorContract, executor)
        ).run(request)

    assert executor.call_count == 0


def test_workflow_rejects_factor_count_too_high_for_eligible_controls_before_executor() -> (
    None
):
    executor = _SpyExecutor()
    request = _request(config=_config(n_unwanted_factors=2))

    with pytest.raises(
        PhosPyInputError,
        match="too few eligible controls.*n_unwanted_factors=2",
    ):
        BatchCorrectionWorkflow(
            executor=cast(BatchCorrectionExecutorContract, executor)
        ).run(request)

    assert executor.call_count == 0


def test_workflow_rejects_factor_count_too_high_for_sample_design_rank_before_executor() -> (
    None
):
    executor = _SpyExecutor()
    request = _request(
        config=_config(n_unwanted_factors=2),
        control_site_set=ControlSiteSet.from_site_keys(
            ("site_a", "site_c", "site_d"),
            source_metadata=_control_source_metadata(),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="factor feasibility.*residual degrees of freedom",
    ):
        BatchCorrectionWorkflow(
            executor=cast(BatchCorrectionExecutorContract, executor)
        ).run(request)

    assert executor.call_count == 0


def test_workflow_rejects_zero_weight_control_before_executor() -> None:
    executor = _SpyExecutor()
    request = _request(
        control_site_set=ControlSiteSet.from_weighted_controls(
            {"site_a": 1.0, "site_c": 0.0},
            source_metadata=_control_source_metadata(),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="control weights must be positive finite values.*'site_c'",
    ):
        BatchCorrectionWorkflow(
            executor=cast(BatchCorrectionExecutorContract, executor)
        ).run(request)

    assert executor.call_count == 0


def test_workflow_design_validation_error_names_sps_ruv_not_linear_residualization() -> (
    None
):
    executor = _SpyExecutor()
    request = _request(sample_metadata=_confounded_sample_metadata())

    with pytest.raises(PhosPyInputError) as exc_info:
        BatchCorrectionWorkflow(
            executor=cast(BatchCorrectionExecutorContract, executor)
        ).run(request)

    message = str(exc_info.value)
    assert "SPS/RUV-style batch correction design validation" in message
    assert "perfectly confounded" in message
    assert "linear_residualize_batch" not in message
    assert executor.call_count == 0


def test_workflow_valid_factor_count_executes_exact_requested_count_without_downshift_warning() -> (
    None
):
    result = BatchCorrectionWorkflow().run(_request())

    executor_diagnostics = cast(Mapping[str, object], result.diagnostics["executor"])
    assert executor_diagnostics["requested_unwanted_factors"] == 1
    assert executor_diagnostics["estimated_unwanted_factors"] == 1
    assert not any(
        "proceeded with the estimated factor count" in warning
        for warning in result.warnings
    )


class _SpyExecutor:
    call_count = 0

    def run(self, **_: object) -> object:
        self.call_count += 1
        raise AssertionError("executor should not be called")


class _SpyInterpreter:
    call_count = 0

    def run(self, **_: object) -> object:
        self.call_count += 1
        raise AssertionError("interpreter should not be called")


@dataclass(frozen=True, slots=True)
class _Diagnostics:
    def to_payload(self) -> dict[str, object]:
        return {"executor_id": "spy_executor"}


@dataclass(frozen=True, slots=True)
class _ExecutorResult:
    corrected_matrix: pd.DataFrame
    diagnostics: _Diagnostics = _Diagnostics()
    warnings: tuple[str, ...] = ()
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
            object.__setattr__(self, "provenance_payload", {})


def _request(
    *,
    phospho: pd.DataFrame | None = None,
    config: InternalBatchCorrectionRequest | None = None,
    missingness_policy: CorrectionMissingnessPolicy | None = None,
    control_site_set: ControlSiteSet | None = None,
    sample_metadata: pd.DataFrame | None = None,
) -> BatchCorrectionWorkflowRequest:
    return BatchCorrectionWorkflowRequest(
        phospho=_phospho() if phospho is None else phospho,
        config=_config() if config is None else config,
        sample_metadata=_sample_metadata()
        if sample_metadata is None
        else sample_metadata,
        control_site_set=(
            ControlSiteSet.from_site_keys(
                ("site_a", "site_c"),
                source_metadata=_control_source_metadata(),
            )
            if control_site_set is None
            else control_site_set
        ),
        missingness_policy=(
            _missingness_policy(missing_cells=())
            if missingness_policy is None
            else missingness_policy
        ),
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


def _config(
    *,
    imputation_policy: InternalBatchCorrectionImputationPolicy = (
        InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY
    ),
    n_unwanted_factors: int = 1,
) -> InternalBatchCorrectionRequest:
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
        imputation_policy=imputation_policy,
        n_unwanted_factors=n_unwanted_factors,
        stage_order=InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM,
        diagnostics_enabled=True,
    )


def _missingness_policy(
    *,
    method: TemporaryImputationMethod = TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
    parameters: dict[str, object] | None = None,
    random_seed: int | None = None,
    missing_cells: tuple[tuple[str, str], ...],
    temporary_imputation: TemporaryImputationPolicy | None = None,
) -> CorrectionMissingnessPolicy:
    if temporary_imputation is None:
        temporary_imputation = TemporaryImputationPolicy(
            allowed=True,
            method=method,
            method_parameters=(
                {"min_observed_values": 2} if parameters is None else parameters
            ),  # type: ignore[arg-type]
            random_seed=random_seed,
        )
    return CorrectionMissingnessPolicy(
        temporary_imputation=temporary_imputation,
        originally_missing_cells_tracked_by=OriginallyMissingCellTracking.OBSERVATION_MASK,
        correction_mask_policy=CorrectionMaskPolicy(),
        observation_mask=ObservationMask(
            feature_ids=tuple(str(value) for value in _phospho().index.tolist()),
            sample_ids=tuple(str(value) for value in _phospho().columns.tolist()),
            originally_missing_cells=missing_cells,
        ),
    )


def _forged_minprob_policy_without_seed() -> TemporaryImputationPolicy:
    policy = object.__new__(TemporaryImputationPolicy)
    object.__setattr__(policy, "allowed", True)
    object.__setattr__(policy, "method", TemporaryImputationMethod.MINPROB_TEMPORARY)
    object.__setattr__(policy, "method_parameters", (("q", 0.01), ("width", 0.3)))
    object.__setattr__(policy, "random_seed", None)
    object.__setattr__(policy, "supported", True)
    object.__setattr__(policy, "unsupported_reason", None)
    object.__setattr__(policy, "imputed_values_are_observed_evidence", False)
    return policy


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [10.0, 5.0, 20.0, 30.0],
            "sample_2": [10.0, 9.0, 20.0, 31.0],
            "sample_3": [14.0, 8.0, 28.0, 42.0],
            "sample_4": [14.0, 12.0, 28.0, 43.0],
        },
        index=pd.Index(("site_a", "site_b", "site_c", "site_d"), name="site_key"),
    )


def _phospho_with_missing() -> pd.DataFrame:
    phospho = _phospho()
    phospho.loc["site_b", "sample_2"] = pd.NA
    return phospho


def _sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("run_1", "run_1", "run_2", "run_2"),
            "condition": ("control", "treated", "control", "treated"),
            "replicate": ("r1", "r1", "r2", "r2"),
        },
        index=("sample_1", "sample_2", "sample_3", "sample_4"),
    )


def _confounded_sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("run_1", "run_1", "run_2", "run_2"),
            "condition": ("control", "control", "treated", "treated"),
            "replicate": ("r1", "r1", "r2", "r2"),
        },
        index=("sample_1", "sample_2", "sample_3", "sample_4"),
    )
