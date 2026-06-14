from __future__ import annotations

import inspect

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.provenance_adapter import (
    PreprocessingProvenanceAdapter,
)
from phospy.science.datasets.preprocessing.stages.missing_data import (
    minprob as minprob_module,
)
from phospy.science.datasets.preprocessing.stages.missing_data.audit import (
    build_knn_audit_records,
    build_minprob_audit_records,
    build_row_median_audit_records,
)
from phospy.science.datasets.preprocessing.stages.missing_data.diagnostics import (
    build_input_profile,
    build_missing_data_diagnostics,
)
from phospy.science.datasets.preprocessing.stages.missing_data.forbid import (
    fail_if_forbid_policy_has_missing_values,
)
from phospy.science.datasets.preprocessing.stages.missing_data.knn import run_knn_policy
from phospy.science.datasets.preprocessing.stages.missing_data.minprob import (
    run_minprob_policy,
)
from phospy.science.datasets.preprocessing.stages.missing_data.row_median import (
    run_row_median_policy,
)
from phospy.science.datasets.preprocessing.stages.missing_data.stage import (
    MissingDataStage,
)


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": [f"G{i}" for i in range(len(index))],
            "site": [f"S{i}" for i in range(len(index))],
            "site_sequence": [f"SEQ_{i}" for i in range(len(index))],
        },
        index=index.copy(),
    )


def _minprob_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [10.0, float("nan"), 12.0, 4.0],
            "sample_b": [9.0, 8.0, float("nan"), 6.0],
            "sample_c": [11.0, 7.0, 5.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute_a", "row_impute_b", "row_impute_c"]),
    )


def _minprob_state(phospho: pd.DataFrame, *, seed: int = 12345) -> PreprocessingState:
    return PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            intensity_transform_policy="log2",
            missing_data_policy="impute_minprob",
            missing_data_q=0.01,
            missing_data_width=0.3,
            missing_data_seed=seed,
            missing_data_max_missing_fraction_per_row=1.0,
            stage_order=("intensity_transform", "missing_data"),
        ),
    )


def _knn_tie_phospho(
    row_order: tuple[str, ...] = ("target", "b_ref", "a_ref", "row_drop"),
) -> pd.DataFrame:
    row_values = {
        "target": {"sample_a": 0.0, "sample_b": float("nan"), "sample_c": 0.0},
        "a_ref": {"sample_a": 0.0, "sample_b": 10.0, "sample_c": 0.0},
        "b_ref": {"sample_a": 0.0, "sample_b": 20.0, "sample_c": 0.0},
        "row_drop": {
            "sample_a": 10.0,
            "sample_b": float("nan"),
            "sample_c": float("nan"),
        },
    }
    return pd.DataFrame(
        [row_values[row_id] for row_id in row_order],
        index=pd.Index(row_order),
        columns=["sample_a", "sample_b", "sample_c"],
    )


def _knn_state(phospho: pd.DataFrame) -> PreprocessingState:
    return PreprocessingState(
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


def _knn_provenance_tables(
    phospho: pd.DataFrame,
) -> tuple[tuple[object, ...], pd.DataFrame, pd.DataFrame]:
    state = _knn_state(phospho)
    final_state, trace = PreprocessingPipeline().run_with_trace(state)
    row_counts, operations = PreprocessingProvenanceAdapter().build_tables(
        plan=state.plan,
        input_row_count=int(len(state.phospho.index)),
        output_row_count=int(len(final_state.phospho.index)),
        trace=trace,
    )
    return trace, row_counts, operations


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


@pytest.mark.reproducibility
def test_knn_imputation_is_reproducible_across_repeated_runs() -> None:
    first = run_knn_policy(_knn_state(_knn_tie_phospho().copy(deep=True)))
    second = run_knn_policy(_knn_state(_knn_tie_phospho().copy(deep=True)))

    pdt.assert_frame_equal(first.phospho, second.phospho)
    pdt.assert_frame_equal(first.imputed_mask, second.imputed_mask)
    assert first.imputed_rows == second.imputed_rows
    assert first.dropped_rows_missing_fraction == second.dropped_rows_missing_fraction


@pytest.mark.reproducibility
def test_knn_imputation_tie_fixture_is_deterministic() -> None:
    b_ref_first = run_knn_policy(
        _knn_state(_knn_tie_phospho(("target", "b_ref", "a_ref", "row_drop")))
    )
    a_ref_first = run_knn_policy(
        _knn_state(_knn_tie_phospho(("target", "a_ref", "b_ref", "row_drop")))
    )

    assert float(b_ref_first.phospho.loc["target", "sample_b"]) == pytest.approx(10.0)
    assert float(a_ref_first.phospho.loc["target", "sample_b"]) == pytest.approx(10.0)
    pdt.assert_frame_equal(
        b_ref_first.phospho.sort_index(),
        a_ref_first.phospho.sort_index(),
    )


@pytest.mark.reproducibility
def test_knn_imputation_diagnostics_are_reproducible() -> None:
    first = MissingDataStage().run(_knn_state(_knn_tie_phospho().copy(deep=True)))
    second = MissingDataStage().run(_knn_state(_knn_tie_phospho().copy(deep=True)))

    assert first.diagnostics == second.diagnostics
    diagnostics = first.diagnostics["diagnostics"]
    assert diagnostics["imputation_method_id"] == "knn"
    assert diagnostics["random_seed"] is None
    assert diagnostics["method_parameters"] == {
        "k": 1,
        "distance": "nan_euclidean",
        "max_missing_fraction_per_row": 0.5,
    }


@pytest.mark.reproducibility
def test_knn_imputation_provenance_is_reproducible() -> None:
    first_trace, first_row_counts, first_operations = _knn_provenance_tables(
        _knn_tie_phospho().copy(deep=True)
    )
    second_trace, second_row_counts, second_operations = _knn_provenance_tables(
        _knn_tie_phospho().copy(deep=True)
    )

    assert first_trace == second_trace
    pdt.assert_frame_equal(first_row_counts, second_row_counts)
    pdt.assert_frame_equal(first_operations, second_operations)
    missing_data_row = first_operations.loc[
        first_operations.loc[:, "stage"] == "missing_data"
    ].iloc[0]
    parameters = missing_data_row["parameters"]
    assert isinstance(parameters, dict)
    assert parameters["missing_data_policy"] == "impute_knn"
    assert parameters["missing_data_k"] == 1
    assert parameters["missing_data_distance"] == "nan_euclidean"
    assert parameters["missing_data_max_missing_fraction_per_row"] == 0.5
    summary = parameters["execution_summary"]
    assert isinstance(summary, dict)
    assert summary["imputation_scope"] == "global_matrix"
    assert summary["diagnostic_summary"]["imputation_method_id"] == "knn"


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


def test_minprob_policy_is_deterministic_for_fixed_seed() -> None:
    first = run_minprob_policy(_minprob_state(_minprob_phospho().copy(deep=True)))
    second = run_minprob_policy(_minprob_state(_minprob_phospho().copy(deep=True)))

    pdt.assert_frame_equal(first.phospho, second.phospho)


def test_minprob_policy_is_stable_under_column_reordering_after_realignment() -> None:
    baseline = run_minprob_policy(_minprob_state(_minprob_phospho().copy(deep=True)))

    reordered_input = (
        _minprob_phospho().loc[:, ["sample_c", "sample_a", "sample_b"]].copy(deep=True)
    )
    reordered = run_minprob_policy(_minprob_state(reordered_input))
    reordered_realigned = reordered.phospho.loc[:, baseline.phospho.columns]

    pdt.assert_frame_equal(baseline.phospho, reordered_realigned)


def test_minprob_policy_is_stable_for_existing_columns_when_unrelated_column_is_inserted() -> (
    None
):
    baseline = run_minprob_policy(_minprob_state(_minprob_phospho().copy(deep=True)))

    with_extra = _minprob_phospho().copy(deep=True)
    with_extra.loc[:, "sample_extra"] = [13.0, 12.0, float("nan"), 11.0]
    with_extra = with_extra.loc[:, ["sample_a", "sample_extra", "sample_b", "sample_c"]]
    with_extra_outcome = run_minprob_policy(_minprob_state(with_extra))
    existing_only = with_extra_outcome.phospho.drop(columns=["sample_extra"]).loc[
        :, baseline.phospho.columns
    ]

    pdt.assert_frame_equal(baseline.phospho, existing_only)


def test_minprob_policy_changes_imputed_values_when_seed_changes() -> None:
    first = run_minprob_policy(
        _minprob_state(_minprob_phospho().copy(deep=True), seed=12345)
    )
    second = run_minprob_policy(
        _minprob_state(_minprob_phospho().copy(deep=True), seed=54321)
    )

    assert float(first.phospho.loc["row_impute_a", "sample_a"]) != pytest.approx(
        float(second.phospho.loc["row_impute_a", "sample_a"])
    )


def test_minprob_policy_changes_column_rng_stream_when_column_label_changes() -> None:
    baseline = run_minprob_policy(_minprob_state(_minprob_phospho().copy(deep=True)))
    renamed_input = _minprob_phospho().rename(columns={"sample_a": "sample_a_renamed"})
    renamed = run_minprob_policy(_minprob_state(renamed_input))

    assert float(baseline.phospho.loc["row_impute_a", "sample_a"]) != pytest.approx(
        float(renamed.phospho.loc["row_impute_a", "sample_a_renamed"])
    )


def test_minprob_stage_diagnostics_include_configured_seed() -> None:
    result = MissingDataStage().run(_minprob_state(_minprob_phospho().copy(deep=True)))
    diagnostics = result.diagnostics["diagnostics"]

    assert diagnostics["random_seed"] == 12345
    assert diagnostics["method_parameters"]["seed"] == 12345


def test_minprob_implementation_does_not_call_python_builtin_hash() -> None:
    source = inspect.getsource(minprob_module)
    assert "hash(" not in source


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
        imputation_mask_hash=None,
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
