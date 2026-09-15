from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from phospy.advanced.configs import (
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
)
from phospy.api.configs import DatasetPreprocessingConfig
from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import DeterminismKind
from phospy.science.datasets.builders.contracts import (
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.science.datasets.builders.provenance_assembler import (
    DatasetRunProvenanceAssembler,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.plan_interpreter import (
    PreprocessingPlanInterpreter,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ImputationInputScale,
    MissingDataPolicy,
)
from phospy.science.datasets.preprocessing.stages.missing_data import (
    stage as missing_data_stage_module,
)
from phospy.science.datasets.preprocessing.stages.missing_data.stage import (
    MISSING_DATA_STAGE_CONTRACT,
    MissingDataStage,
)
from phospy.science.transformations.models import IntensityScaleKind


def _valid_group_aware_kwargs() -> dict[str, Any]:
    return {
        "policy": "impute_group_aware",
        "group_column": "condition",
        "min_partial_observed_fraction": 0.5,
        "min_reference_observed_fraction": 0.75,
        "q": 0.01,
        "width": 0.3,
        "seed": 12345,
        "k": 3,
        "distance": "nan_euclidean",
    }


def _valid_direct_plan_kwargs() -> dict[str, Any]:
    return {
        "intensity_transform_policy": "log2",
        "missing_data_policy": "impute_group_aware",
        "missing_data_group_column": "condition",
        "missing_data_min_partial_observed_fraction": 0.5,
        "missing_data_min_reference_observed_fraction": 0.75,
        "missing_data_q": 0.01,
        "missing_data_width": 0.3,
        "missing_data_seed": 12345,
        "missing_data_k": 3,
        "missing_data_distance": "nan_euclidean",
        "stage_order": ("intensity_transform", "missing_data"),
    }


def test_group_aware_public_config_is_accepted() -> None:
    config = DatasetMissingDataConfig(**_valid_group_aware_kwargs())

    assert config.policy == "impute_group_aware"
    assert config.group_column == "condition"
    assert config.no_overlap_policy == "error"
    assert config.max_missing_fraction_per_row is None


@pytest.mark.parametrize(
    "field_name",
    [
        "group_column",
        "min_partial_observed_fraction",
        "min_reference_observed_fraction",
        "q",
        "width",
        "seed",
        "k",
        "distance",
    ],
)
def test_group_aware_required_fields_are_enforced(field_name: str) -> None:
    kwargs = _valid_group_aware_kwargs()
    kwargs[field_name] = None

    with pytest.raises(PhosPyInputError, match=field_name):
        DatasetMissingDataConfig(**kwargs)


@pytest.mark.parametrize(
    ("field_name", "value"),
    [
        ("min_partial_observed_fraction", 0.0),
        ("min_partial_observed_fraction", 1.01),
        ("min_reference_observed_fraction", -0.1),
        ("min_reference_observed_fraction", 2.0),
    ],
)
def test_group_aware_thresholds_require_open_zero_closed_one(
    field_name: str,
    value: float,
) -> None:
    kwargs = _valid_group_aware_kwargs()
    kwargs[field_name] = value

    with pytest.raises(PhosPyInputError, match=r"0 < value <= 1"):
        DatasetMissingDataConfig(**kwargs)


def test_group_aware_rejects_non_nan_euclidean_distance() -> None:
    kwargs = _valid_group_aware_kwargs()
    kwargs["distance"] = "euclidean"

    with pytest.raises(PhosPyInputError, match="nan_euclidean"):
        DatasetMissingDataConfig(**kwargs)


def test_group_aware_rejects_linear_input_scale() -> None:
    kwargs = _valid_group_aware_kwargs()
    kwargs["input_scale"] = "linear"

    with pytest.raises(PhosPyInputError, match="input_scale.*log2"):
        DatasetMissingDataConfig(**kwargs)


def test_group_aware_rejects_standalone_row_eligibility_control() -> None:
    kwargs = _valid_group_aware_kwargs()
    kwargs["max_missing_fraction_per_row"] = 0.5

    with pytest.raises(PhosPyInputError, match="max_missing_fraction_per_row"):
        DatasetMissingDataConfig(**kwargs)


@pytest.mark.parametrize(
    "no_overlap_policy",
    ["column_mean_with_caveat", "unsupported_fallback"],
)
def test_group_aware_rejects_unsupported_knn_fallback(
    no_overlap_policy: str,
) -> None:
    kwargs = _valid_group_aware_kwargs()
    kwargs["no_overlap_policy"] = no_overlap_policy

    with pytest.raises(PhosPyInputError, match="no_overlap_policy"):
        DatasetMissingDataConfig(**kwargs)


def test_standalone_knn_and_minprob_defaults_remain_unchanged() -> None:
    knn = DatasetMissingDataConfig(
        policy="impute_knn",
        k=2,
        distance="nan_euclidean",
        max_missing_fraction_per_row=0.5,
        input_scale="linear",
    )
    minprob = DatasetMissingDataConfig(
        policy="impute_minprob",
        q=0.01,
        width=0.3,
        seed=12345,
        max_missing_fraction_per_row=0.5,
    )

    assert knn.no_overlap_policy == "column_mean_with_caveat"
    assert minprob.input_scale is None


def test_group_aware_plan_resolves_log2_order_and_seeded_contract() -> None:
    plan = PreprocessingPlanInterpreter().run(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            missing_data=DatasetMissingDataConfig(**_valid_group_aware_kwargs()),
        )
    )

    assert plan.missing_data_policy is MissingDataPolicy.IMPUTE_GROUP_AWARE
    assert plan.missing_data_input_scale is ImputationInputScale.LOG2
    assert plan.missing_data_input_scale_source == "method_required"
    assert plan.stage_order.index("intensity_transform") < plan.stage_order.index(
        "missing_data"
    )
    interpreted = MISSING_DATA_STAGE_CONTRACT.interpret(plan)
    assert interpreted.determinism_kind is DeterminismKind.SEEDED_STOCHASTIC
    assert interpreted.quantitative_contract.accepted_input_scale_kinds == frozenset(
        {IntensityScaleKind.LOG2}
    )
    assert interpreted.parameters["missing_data_group_column"] == "condition"


def test_group_aware_declared_log2_input_does_not_require_transform() -> None:
    plan = PreprocessingPlanInterpreter().run(
        DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(**_valid_group_aware_kwargs()),
        ),
        declared_input_scale_kind=IntensityScaleKind.LOG2,
    )

    assert "intensity_transform" not in plan.stage_order
    assert plan.stage_order[-1] == "missing_data"


def test_group_aware_direct_plan_defaults_no_overlap_and_resolves_log2_order() -> None:
    plan = PreprocessingPlan(**_valid_direct_plan_kwargs())

    assert plan.missing_data_policy is MissingDataPolicy.IMPUTE_GROUP_AWARE
    assert plan.missing_data_no_overlap_policy == "error"
    assert plan.missing_data_input_scale is ImputationInputScale.LOG2
    assert plan.missing_data_input_scale_source == "method_required"
    assert plan.missing_data_imputation_operation_order == "after_intensity_transform"
    assert plan.stage_order == ("intensity_transform", "missing_data")
    assert tuple(item.stage for item in plan.stage_order_resolution) == plan.stage_order


@pytest.mark.parametrize(
    "field_name",
    [
        "missing_data_group_column",
        "missing_data_min_partial_observed_fraction",
        "missing_data_min_reference_observed_fraction",
        "missing_data_q",
        "missing_data_width",
        "missing_data_seed",
        "missing_data_k",
        "missing_data_distance",
    ],
)
def test_group_aware_direct_plan_requires_contract_fields(field_name: str) -> None:
    kwargs = _valid_direct_plan_kwargs()
    kwargs[field_name] = None

    with pytest.raises(PhosPyInputError, match=field_name):
        PreprocessingPlan(**kwargs)


@pytest.mark.parametrize(
    ("field_name", "value", "message"),
    [
        ("missing_data_min_partial_observed_fraction", 0.0, r"0 < value <= 1"),
        ("missing_data_min_reference_observed_fraction", 1.1, r"0 < value <= 1"),
        ("missing_data_input_scale", "linear", "requires.*log2"),
        ("missing_data_max_missing_fraction_per_row", 0.5, "row eligibility"),
        ("missing_data_distance", "euclidean", "nan_euclidean"),
        (
            "missing_data_no_overlap_policy",
            "column_mean_with_caveat",
            "column_mean_with_caveat",
        ),
    ],
)
def test_group_aware_direct_plan_rejects_incompatible_fields(
    field_name: str,
    value: object,
    message: str,
) -> None:
    kwargs = _valid_direct_plan_kwargs()
    kwargs[field_name] = value

    with pytest.raises(PhosPyInputError, match=message):
        PreprocessingPlan(**kwargs)


def test_group_aware_execution_is_guarded_without_mutating_state(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mechanism_calls: list[str] = []

    def _record_knn_call(*args: object, **kwargs: object) -> None:
        del args, kwargs
        mechanism_calls.append("knn")

    def _record_minprob_call(*args: object, **kwargs: object) -> None:
        del args, kwargs
        mechanism_calls.append("minprob")

    monkeypatch.setattr(missing_data_stage_module, "run_knn_policy", _record_knn_call)
    monkeypatch.setattr(
        missing_data_stage_module,
        "run_minprob_policy",
        _record_minprob_call,
    )
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [float("nan")]},
        index=pd.Index(["row_a"], name="site_key"),
    )
    site_metadata = pd.DataFrame(
        {"site": ["S1"]},
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"condition": ["A", "B"]},
        index=pd.Index(phospho.columns, name="sample"),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        plan=PreprocessingPlan(**_valid_direct_plan_kwargs()),
    )
    original_phospho = phospho.copy(deep=True)
    original_site_metadata = site_metadata.copy(deep=True)
    original_sample_metadata = sample_metadata.copy(deep=True)

    with pytest.raises(
        PhosPyInputError,
        match=r"planning/contract-only.*not implemented",
    ) as exc_info:
        MissingDataStage().run(state)

    assert "unsupported" not in str(exc_info.value).lower()
    assert mechanism_calls == []
    pd.testing.assert_frame_equal(state.phospho, original_phospho)
    pd.testing.assert_frame_equal(state.site_metadata, original_site_metadata)
    pd.testing.assert_frame_equal(state.sample_metadata, original_sample_metadata)
    assert state.row_audit is None
    assert state.imputation_observation_mask is None
    assert state.report_rows == ()


def test_group_aware_plan_is_serialized_through_provenance_assembler() -> None:
    plan = PreprocessingPlanInterpreter().run(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            missing_data=DatasetMissingDataConfig(**_valid_group_aware_kwargs()),
        )
    )
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["row_a"], name="site_key"),
    )
    site_metadata = pd.DataFrame(
        {"site": ["S1"]},
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"condition": ["A", "B"]},
        index=pd.Index(phospho.columns, name="sample"),
    )
    request = InterpretedDatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        organism=None,
        preprocessing_plan=plan,
    )
    preprocessed = PreprocessedDatasetBuildTables(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=None,
        preprocessing_trace=None,
    )

    provenance = DatasetRunProvenanceAssembler().run(
        request=request,
        preprocessed=preprocessed,
        validated_site_metadata=site_metadata,
        resolved_phospho=phospho,
        resolved_total=None,
        preprocessing_trace=None,
        intensity_scale_label="log2",
        intensity_scale_establishment={"source": "test"},
        quantitative_meaning_provenance={
            "schema_version": 1,
            "source_quantity": None,
            "target_quantity": "phosphosite_abundance",
            "operation_id": "test.group_aware",
            "producer_id": "test",
            "evidence_mode": "test",
            "parameters": {},
            "input_table_fingerprints": {},
            "output_table_fingerprint": None,
            "trace_id": None,
            "diagnostic_caveat_codes": [],
        },
        quantitative_meaning="phosphosite_abundance",
        allow_opaque_site_values=False,
    )

    payload = provenance.workflow_parameters["preprocessing_plan"]
    assert payload["missing_data_group_column"] == "condition"
    assert payload["missing_data_min_partial_observed_fraction"] == 0.5
    assert payload["missing_data_min_reference_observed_fraction"] == 0.75
    assert payload["missing_data_q"] == 0.01
    assert payload["missing_data_width"] == 0.3
    assert payload["missing_data_seed"] == 12345
    assert payload["missing_data_k"] == 3
    assert payload["missing_data_distance"] == "nan_euclidean"
    assert payload["missing_data_no_overlap_policy"] == "error"
    assert payload["missing_data_input_scale"] == "log2"
    assert payload["missing_data_input_scale_source"] == "method_required"
    assert (
        payload["missing_data_imputation_operation_order"]
        == "after_intensity_transform"
    )
    assert payload["stage_order"] == list(plan.stage_order)
    assert [item["stage"] for item in payload["resolved_stage_order"]] == list(
        plan.stage_order
    )
