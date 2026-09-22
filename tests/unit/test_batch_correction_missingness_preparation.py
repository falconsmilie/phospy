from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.provenance import fingerprint_matrix
from phospy.science.batch_correction.executor import (
    _prepare_matrix as prepare_sps_ruv_matrix,
)
from phospy.science.batch_correction.ruv_iii_executor import (
    _prepare_matrix as prepare_ruv_iii_matrix,
)
from phospy.science.configs.preprocessing import TemporaryImputationMethod


@dataclass(frozen=True)
class _ObservationMask:
    feature_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    originally_missing_cells: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class _TemporaryPolicy:
    allowed: bool = True
    method: TemporaryImputationMethod = TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY
    method_parameters: tuple[tuple[str, int], ...] = (("min_observed_values", 2),)


@dataclass(frozen=True)
class _Plan:
    observation_mask: _ObservationMask
    temporary_imputation_policy: _TemporaryPolicy = _TemporaryPolicy()


@dataclass(frozen=True)
class _ReferencePrepared:
    working: pd.DataFrame
    originally_missing: pd.DataFrame
    actual_missing: pd.DataFrame
    missing_cells: tuple[tuple[str, str], ...]
    completion_applied: bool


def _reference_prepare(
    matrix: pd.DataFrame,
    plan: _Plan,
) -> _ReferencePrepared:
    """Retain the scalar 1.7.5 mechanics as a test-only equivalence oracle."""
    values = matrix.astype("float64").copy(deep=True)
    originally_missing = pd.DataFrame(
        False,
        index=matrix.index.copy(),
        columns=matrix.columns.copy(),
    )
    for feature_id, sample_id in plan.observation_mask.originally_missing_cells:
        originally_missing.loc[feature_id, sample_id] = True
    actual_missing = values.isna()
    if bool(actual_missing.to_numpy().any()):
        minimum = dict(plan.temporary_imputation_policy.method_parameters)[
            "min_observed_values"
        ]
        for row_id in values.index.tolist():
            row = values.loc[row_id, :]
            observed = row.dropna()
            if int(observed.shape[0]) < minimum:
                raise PhosPyInputError(str(row_id))
            missing = row.isna()
            if bool(missing.any()):
                values.loc[row_id, missing] = float(observed.median())
    return _ReferencePrepared(
        working=values,
        originally_missing=originally_missing,
        actual_missing=actual_missing,
        missing_cells=plan.observation_mask.originally_missing_cells,
        completion_applied=bool(actual_missing.to_numpy().any()),
    )


def _matrix_and_plan(
    case: str,
    *,
    reordered: bool = False,
) -> tuple[pd.DataFrame, _Plan]:
    matrix = pd.DataFrame(
        np.arange(1.0, 21.0).reshape(5, 4),
        index=pd.Index(["f1", "f2", "f3", "f4", "f5"], name="feature"),
        columns=pd.Index(["s1", "s2", "s3", "s4"], name="sample"),
    )
    governed_by_case = {
        "complete": (),
        "one": (("f3", "s2"),),
        "scattered": (("f1", "s4"), ("f3", "s2"), ("f5", "s1")),
        "dense": (
            ("f1", "s1"),
            ("f1", "s3"),
            ("f2", "s2"),
            ("f2", "s4"),
            ("f3", "s1"),
            ("f3", "s2"),
            ("f4", "s3"),
            ("f4", "s4"),
            ("f5", "s1"),
            ("f5", "s4"),
        ),
        "upstream": (("f2", "s2"), ("f4", "s3")),
        "mixed": (("f1", "s4"), ("f2", "s2"), ("f4", "s3")),
    }
    governed = governed_by_case[case]
    actual = governed
    if case == "upstream":
        actual = ()
    elif case == "mixed":
        actual = (governed[0], governed[2])
    for feature_id, sample_id in actual:
        matrix.loc[feature_id, sample_id] = np.nan
    if reordered:
        matrix = matrix.loc[
            ["f5", "f2", "f4", "f1", "f3"],
            ["s3", "s1", "s4", "s2"],
        ]
    plan = _Plan(
        _ObservationMask(
            tuple(str(value) for value in matrix.index),
            tuple(str(value) for value in matrix.columns),
            governed,
        )
    )
    return matrix, plan


@pytest.mark.parametrize(
    ("case", "reordered"),
    [
        ("complete", False),
        ("one", False),
        ("scattered", False),
        ("dense", False),
        ("upstream", False),
        ("mixed", False),
        ("mixed", True),
    ],
)
def test_vectorized_preparation_exactly_matches_scalar_reference(
    case: str,
    reordered: bool,
) -> None:
    matrix, plan = _matrix_and_plan(case, reordered=reordered)
    reference = _reference_prepare(matrix, plan)

    ruv_iii = prepare_ruv_iii_matrix(phospho=matrix, plan=plan)  # type: ignore[arg-type]
    sps_ruv = prepare_sps_ruv_matrix(phospho=matrix, plan=plan)  # type: ignore[arg-type]

    for prepared in (ruv_iii, sps_ruv):
        pd.testing.assert_frame_equal(
            prepared.working,
            reference.working,
            check_exact=True,
        )
        pd.testing.assert_frame_equal(
            prepared.originally_missing,
            reference.originally_missing,
            check_exact=True,
        )
        pd.testing.assert_frame_equal(
            prepared.actual_missing,
            reference.actual_missing,
            check_exact=True,
        )
        assert prepared.missing_cells == reference.missing_cells
        assert (
            fingerprint_matrix(
                prepared.working,
                name="working",
            ).exact_hash_value
            == fingerprint_matrix(
                reference.working,
                name="working",
            ).exact_hash_value
        )
        assert (
            fingerprint_matrix(
                prepared.originally_missing.astype("int8"),
                name="originally-missing",
            ).exact_hash_value
            == fingerprint_matrix(
                reference.originally_missing.astype("int8"),
                name="originally-missing",
            ).exact_hash_value
        )
    assert ruv_iii.temporary_completion_applied is reference.completion_applied
    assert bool(sps_ruv.warnings) is reference.completion_applied


@pytest.mark.parametrize("axis", ["feature", "sample"])
@pytest.mark.parametrize(
    ("prepare", "message"),
    [
        (prepare_ruv_iii_matrix, "observation mask axes"),
        (prepare_sps_ruv_matrix, "observation mask .* must match"),
    ],
)
def test_preparation_rejects_mismatched_axes(
    axis: str,
    prepare: object,
    message: str,
) -> None:
    matrix, plan = _matrix_and_plan("complete")
    mask = plan.observation_mask
    if axis == "feature":
        mask = _ObservationMask(tuple(reversed(mask.feature_ids)), mask.sample_ids, ())
    else:
        mask = _ObservationMask(mask.feature_ids, tuple(reversed(mask.sample_ids)), ())

    with pytest.raises(PhosPyInputError, match=message):
        prepare(phospho=matrix, plan=_Plan(mask))  # type: ignore[operator]


@pytest.mark.parametrize("prepare", [prepare_ruv_iii_matrix, prepare_sps_ruv_matrix])
def test_preparation_rejects_actual_missing_cell_absent_from_mask(
    prepare: object,
) -> None:
    matrix, plan = _matrix_and_plan("complete")
    matrix.loc["f3", "s2"] = np.nan

    with pytest.raises(
        PhosPyInputError, match="missing cells|missing cells not governed"
    ):
        prepare(phospho=matrix, plan=plan)  # type: ignore[operator]


@pytest.mark.parametrize("prepare", [prepare_ruv_iii_matrix, prepare_sps_ruv_matrix])
@pytest.mark.parametrize("non_finite", [np.inf, -np.inf])
def test_preparation_rejects_observed_infinity(
    prepare: object,
    non_finite: float,
) -> None:
    matrix, plan = _matrix_and_plan("complete")
    matrix.loc["f2", "s3"] = non_finite

    with pytest.raises(PhosPyInputError, match="finite"):
        prepare(phospho=matrix, plan=plan)  # type: ignore[operator]


def test_sps_ruv_preparation_reports_first_non_finite_cell_in_axis_order() -> None:
    matrix, plan = _matrix_and_plan("complete")
    matrix.loc["f4", "s1"] = -np.inf
    matrix.loc["f2", "s4"] = np.inf

    with pytest.raises(
        PhosPyInputError,
        match="feature_id='f2', sample_id='s4'",
    ):
        prepare_sps_ruv_matrix(phospho=matrix, plan=plan)  # type: ignore[arg-type]


@pytest.mark.parametrize("prepare", [prepare_ruv_iii_matrix, prepare_sps_ruv_matrix])
def test_preparation_reports_first_insufficient_row_in_axis_order(
    prepare: object,
) -> None:
    matrix, _ = _matrix_and_plan("complete")
    governed = (
        ("f2", "s1"),
        ("f2", "s2"),
        ("f2", "s3"),
        ("f2", "s4"),
        ("f4", "s1"),
        ("f4", "s2"),
        ("f4", "s3"),
        ("f4", "s4"),
    )
    for feature_id, sample_id in governed:
        matrix.loc[feature_id, sample_id] = np.nan
    plan = _Plan(_ObservationMask(tuple(matrix.index), tuple(matrix.columns), governed))

    with pytest.raises(PhosPyInputError, match="row 'f2'"):
        prepare(phospho=matrix, plan=plan)  # type: ignore[operator]
