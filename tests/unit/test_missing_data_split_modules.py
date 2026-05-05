from __future__ import annotations

import pandas as pd
import pytest

from phospy.datasets.preprocessing.models import PreprocessingPlan, PreprocessingState
from phospy.datasets.preprocessing.stages.missing_data.audit import (
    build_knn_audit_records,
    build_minprob_audit_records,
    build_row_median_audit_records,
)
from phospy.datasets.preprocessing.stages.missing_data.diagnostics import (
    build_input_profile,
    build_missing_data_diagnostics,
)
from phospy.datasets.preprocessing.stages.missing_data.forbid import (
    fail_if_forbid_policy_has_missing_values,
)
from phospy.datasets.preprocessing.stages.missing_data.knn import run_knn_policy
from phospy.datasets.preprocessing.stages.missing_data.minprob import run_minprob_policy
from phospy.datasets.preprocessing.stages.missing_data.row_median import (
    run_row_median_policy,
)
from phospy.datasets.preprocessing.stages.missing_data.stage import MissingDataStage
from phospy.errors.input import PhosPyInputError


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": [f"G{i}" for i in range(len(index))],
            "site": [f"S{i}" for i in range(len(index))],
            "site_sequence": [f"SEQ_{i}" for i in range(len(index))],
        },
        index=index.copy(),
    )


def test_forbid_policy_rejects_missing_values_without_full_pipeline() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, float("nan")], "sample_b": [2.0, 3.0]},
        index=pd.Index(["row_a", "row_b"]),
    )

    with pytest.raises(PhosPyInputError, match="missing_data.policy='forbid'"):
        fail_if_forbid_policy_has_missing_values(phospho)


def test_row_median_policy_is_independently_testable() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), float("nan")],
            "sample_b": [2.0, 10.0, float("nan")],
            "sample_c": [3.0, 20.0, 9.0],
            "sample_d": [4.0, float("nan"), float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute", "row_drop"]),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_row_median",
            missing_data_min_observed_values=2,
            stage_order=("missing_data",),
        ),
    )

    outcome = run_row_median_policy(state)

    assert outcome.phospho.index.tolist() == ["row_keep", "row_impute"]
    assert float(outcome.phospho.loc["row_impute", "sample_a"]) == 15.0
    assert float(outcome.phospho.loc["row_impute", "sample_d"]) == 15.0
    assert outcome.dropped_row_ids == ("row_drop",)


def test_knn_policy_is_independently_testable() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 1.0, 2.0, 10.0],
            "sample_b": [1.0, 2.0, 2.0, float("nan")],
            "sample_c": [float("nan"), 3.0, 3.0, float("nan")],
        },
        index=pd.Index(["row_impute", "row_ref_1", "row_ref_2", "row_drop"]),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            missing_data_policy="impute_knn",
            missing_data_k=1,
            missing_data_distance="nan_euclidean",
            missing_data_max_missing_fraction_per_row=0.5,
            stage_order=("missing_data",),
        ),
    )

    outcome = run_knn_policy(state)

    assert float(outcome.phospho.loc["row_impute", "sample_c"]) == pytest.approx(3.0)
    assert int(outcome.phospho.isna().to_numpy().sum()) == 0
    assert outcome.dropped_row_ids == ("row_drop",)


def test_minprob_policy_is_independently_testable() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), float("nan"), 4.0],
            "sample_b": [9.0, 8.0, float("nan"), 6.0],
            "sample_c": [11.0, 7.0, 5.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute_a", "row_drop", "row_impute_c"]),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            intensity_transform_policy="log2",
            missing_data_policy="impute_minprob",
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=12345,
            missing_data_max_missing_fraction_per_row=0.5,
            stage_order=("intensity_transform", "missing_data"),
        ),
    )

    outcome = run_minprob_policy(state)

    assert int(outcome.phospho.isna().to_numpy().sum()) == 0
    assert outcome.dropped_row_ids == ("row_drop",)
    assert "sample_a" in outcome.per_column_distribution_parameters


def test_diagnostics_builder_supports_each_policy_shape() -> None:
    diagnostics = build_missing_data_diagnostics(
        missing_data_policy="impute_row_median",
        imputation_method_id="row_median",
        imputation_method_family="deterministic_row_statistic",
        input_missing_cell_count=3,
        output_missing_cell_count=0,
        imputed_cell_count=3,
        affected_row_ids=("r1", "r2"),
        affected_column_ids=("c1",),
        imputed_row_ids=("r1",),
        imputed_column_ids=("c1",),
        dropped_row_ids=("r2",),
        random_seed=None,
        method_parameters={"min_observed_values": 1},
        matrix_scale_requirement=None,
        stage_order=("missing_data",),
        missingness_mask_hash="abc",
        left_censored_assumption=False,
        rows_not_imputable=(),
        row_medians_used={"r1": 2.0},
        per_column_distribution_parameters=None,
        dropped_rows_above_max_missing_fraction=(),
        neighbour_count=None,
        distance_metric=None,
    )

    assert diagnostics["missing_data_policy"] == "impute_row_median"
    assert diagnostics["imputation_method_id"] == "row_median"
    assert diagnostics["row_medians_used"] == {"r1": 2.0}


def test_row_audit_builders_are_independently_testable() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), float("nan"), 4.0],
            "sample_b": [2.0, 2.0, float("nan"), 5.0],
            "sample_c": [3.0, float("nan"), 9.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute", "row_drop", "row_impute_2"]),
    )
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            intensity_transform_policy="log2",
            missing_data_policy="impute_minprob",
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=7,
            missing_data_max_missing_fraction_per_row=0.5,
            stage_order=("intensity_transform", "missing_data"),
        ),
    )
    input_profile = build_input_profile(state.phospho)

    row_median_outcome = run_row_median_policy(
        replace_plan(state, missing_data_policy="impute_row_median", min_obs=2)
    )
    knn_outcome = run_knn_policy(
        replace_plan(
            state, missing_data_policy="impute_knn", k=1, distance="nan_euclidean"
        )
    )
    minprob_outcome = run_minprob_policy(state)

    row_median_records = build_row_median_audit_records(
        plan=replace_plan(
            state, missing_data_policy="impute_row_median", min_obs=2
        ).plan,
        input_profile=input_profile,
        outcome=row_median_outcome,
    )
    knn_records = build_knn_audit_records(
        plan=replace_plan(
            state, missing_data_policy="impute_knn", k=1, distance="nan_euclidean"
        ).plan,
        input_profile=input_profile,
        outcome=knn_outcome,
    )
    minprob_records = build_minprob_audit_records(
        plan=state.plan,
        input_profile=input_profile,
        outcome=minprob_outcome,
    )

    assert row_median_records
    assert knn_records
    assert minprob_records
    assert all(record.stage == "missing_data" for record in row_median_records)
    assert all(record.stage == "missing_data" for record in knn_records)
    assert all(record.stage == "missing_data" for record in minprob_records)


def test_unsupported_policy_fails_clearly() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["row_a"]),
    )
    plan = PreprocessingPlan(
        missing_data_policy="forbid",
        stage_order=("missing_data",),
    )
    object.__setattr__(plan, "missing_data_policy", "unsupported_policy")
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=plan,
    )

    with pytest.raises(PhosPyInputError, match="unsupported missing_data.policy"):
        MissingDataStage().run(state)


def replace_plan(
    state: PreprocessingState,
    *,
    missing_data_policy: str,
    min_obs: int | None = None,
    k: int | None = None,
    distance: str | None = None,
) -> PreprocessingState:
    return PreprocessingState(
        phospho=state.phospho.copy(deep=True),
        site_metadata=state.site_metadata.copy(deep=True),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            intensity_transform_policy=state.plan.intensity_transform_policy,
            missing_data_policy=missing_data_policy,
            missing_data_min_observed_values=min_obs,
            missing_data_q=state.plan.missing_data_q,
            missing_data_width=state.plan.missing_data_width,
            missing_data_seed=state.plan.missing_data_seed,
            missing_data_max_missing_fraction_per_row=state.plan.missing_data_max_missing_fraction_per_row,
            missing_data_k=k,
            missing_data_distance=distance,
            stage_order=state.plan.stage_order,
        ),
    )
