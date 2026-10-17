from __future__ import annotations

import statistics
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest

from benchmarks.measure_ruv_iii_missingness_preparation import _source_provenance
from phospy.preprocessing.batch_correction.contracts.missingness import (
    TemporaryImputationMethod,
)
from phospy.preprocessing.batch_correction.execution.ruv_iii_executor import (
    _materialize_observation_mask,
    _prepare_matrix,
    _temporary_row_median_completion,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]


@dataclass(frozen=True)
class _Mask:
    feature_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    originally_missing_cells: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class _Policy:
    allowed: bool = True
    method: TemporaryImputationMethod = TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY
    method_parameters: tuple[tuple[str, int], ...] = (("min_observed_values", 2),)


@dataclass(frozen=True)
class _Plan:
    observation_mask: _Mask
    temporary_imputation_policy: _Policy = _Policy()


def _fixture(*, case: str) -> tuple[pd.DataFrame, _Plan]:
    rows, columns = 1_000, 20
    rng = np.random.default_rng(72_101)
    feature_ids = tuple(f"feature_{position:05d}" for position in range(rows))
    sample_ids = tuple(f"sample_{position:03d}" for position in range(columns))
    values = rng.normal(size=(rows, columns))
    governed_count = 0 if case == "complete" else 2_000
    flat = np.sort(rng.choice(values.size, size=governed_count, replace=False))
    row_positions, column_positions = np.unravel_index(flat, values.shape)
    governed = tuple(
        (feature_ids[int(row)], sample_ids[int(column)])
        for row, column in zip(row_positions, column_positions, strict=True)
    )
    if case == "actual_missing":
        values[row_positions, column_positions] = np.nan
    matrix = pd.DataFrame(values, index=feature_ids, columns=sample_ids)
    return matrix, _Plan(_Mask(feature_ids, sample_ids, governed))


def _scalar_reference(matrix: pd.DataFrame, plan: _Plan) -> None:
    working = matrix.copy(deep=True)
    governed = pd.DataFrame(False, index=matrix.index, columns=matrix.columns)
    for feature_id, sample_id in plan.observation_mask.originally_missing_cells:
        governed.loc[feature_id, sample_id] = True
    if not bool(working.isna().to_numpy().any()):
        return
    for row_id in working.index:
        row = working.loc[row_id, :]
        missing = row.isna()
        if bool(missing.any()):
            working.loc[row_id, missing] = float(row.dropna().median())


def _scalar_mask(matrix: pd.DataFrame, plan: _Plan) -> pd.DataFrame:
    result = pd.DataFrame(False, index=matrix.index, columns=matrix.columns)
    for feature_id, sample_id in plan.observation_mask.originally_missing_cells:
        result.loc[feature_id, sample_id] = True
    return result


def _scalar_row_median(values: np.ndarray, actual_missing: np.ndarray) -> np.ndarray:
    result = values.copy()
    for row_position in range(result.shape[0]):
        missing = actual_missing[row_position]
        if bool(missing.any()):
            result[row_position, missing] = float(np.nanmedian(result[row_position, :]))
    return result


def _median_seconds(callable_: object, *, repeats: int) -> float:
    timings: list[float] = []
    for _ in range(repeats):
        started = time.perf_counter()
        callable_()  # type: ignore[operator]
        timings.append(time.perf_counter() - started)
    return float(statistics.median(timings))


@pytest.mark.parametrize("case", ["upstream_imputed", "actual_missing"])
def test_preparation_remains_materially_faster_than_scalar_reference(
    case: str,
    record_property: object,
) -> None:
    matrix, plan = _fixture(case=case)
    _prepare_matrix(phospho=matrix, plan=plan)  # type: ignore[arg-type]
    vectorized = _median_seconds(
        lambda: _prepare_matrix(phospho=matrix, plan=plan),  # type: ignore[arg-type]
        repeats=3,
    )
    scalar = _median_seconds(lambda: _scalar_reference(matrix, plan), repeats=2)
    speedup = scalar / vectorized
    record_property("ruv_iii_preparation_vectorized_seconds", f"{vectorized:.6f}")  # type: ignore[operator]
    record_property("ruv_iii_preparation_scalar_seconds", f"{scalar:.6f}")  # type: ignore[operator]
    record_property("ruv_iii_preparation_speedup", f"{speedup:.3f}")  # type: ignore[operator]

    assert speedup >= 3.0


def test_complete_input_keeps_the_empty_mask_fast_path(
    record_property: object,
) -> None:
    complete_matrix, complete_plan = _fixture(case="complete")
    upstream_matrix, upstream_plan = _fixture(case="upstream_imputed")
    complete = _median_seconds(
        lambda: _prepare_matrix(  # type: ignore[arg-type]
            phospho=complete_matrix,
            plan=complete_plan,
        ),
        repeats=5,
    )
    governed = _median_seconds(
        lambda: _prepare_matrix(  # type: ignore[arg-type]
            phospho=upstream_matrix,
            plan=upstream_plan,
        ),
        repeats=5,
    )
    ratio = complete / governed
    record_property("ruv_iii_complete_seconds", f"{complete:.6f}")  # type: ignore[operator]
    record_property("ruv_iii_governed_numeric_seconds", f"{governed:.6f}")  # type: ignore[operator]
    record_property("ruv_iii_complete_to_governed_ratio", f"{ratio:.3f}")  # type: ignore[operator]

    assert ratio <= 1.25


def test_mask_materialisation_is_measured_independently(
    record_property: object,
) -> None:
    matrix, plan = _fixture(case="upstream_imputed")
    expected = _scalar_mask(matrix, plan)
    actual = _materialize_observation_mask(phospho=matrix, plan=plan)  # type: ignore[arg-type]
    pd.testing.assert_frame_equal(actual, expected)

    vectorized = _median_seconds(
        lambda: _materialize_observation_mask(  # type: ignore[arg-type]
            phospho=matrix,
            plan=plan,
        ),
        repeats=5,
    )
    scalar = _median_seconds(lambda: _scalar_mask(matrix, plan), repeats=3)
    speedup = scalar / vectorized
    record_property("ruv_iii_mask_materialisation_seconds", f"{vectorized:.6f}")  # type: ignore[operator]
    record_property("ruv_iii_mask_scalar_seconds", f"{scalar:.6f}")  # type: ignore[operator]
    record_property("ruv_iii_mask_speedup", f"{speedup:.3f}")  # type: ignore[operator]

    assert speedup >= 3.0


def test_row_median_completion_is_measured_independently(
    record_property: object,
) -> None:
    matrix, plan = _fixture(case="actual_missing")
    original = matrix.to_numpy(dtype="float64", copy=True)
    missing = np.isnan(original)
    expected = _scalar_row_median(original, missing)
    actual = original.copy()
    assert _temporary_row_median_completion(
        values=actual,
        actual_missing=missing,
        index=matrix.index,
        policy=plan.temporary_imputation_policy,
    )
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(matrix.to_numpy(), original)

    vectorized = _median_seconds(
        lambda: _temporary_row_median_completion(
            values=original.copy(),
            actual_missing=missing,
            index=matrix.index,
            policy=plan.temporary_imputation_policy,
        ),
        repeats=5,
    )
    scalar = _median_seconds(
        lambda: _scalar_row_median(original, missing),
        repeats=3,
    )
    speedup = scalar / vectorized
    record_property("ruv_iii_row_median_completion_seconds", f"{vectorized:.6f}")  # type: ignore[operator]
    record_property("ruv_iii_row_median_scalar_seconds", f"{scalar:.6f}")  # type: ignore[operator]
    record_property("ruv_iii_row_median_speedup", f"{speedup:.3f}")  # type: ignore[operator]

    assert speedup >= 3.0


def test_benchmark_provenance_hashes_current_owner_sources() -> None:
    payload = _source_provenance()

    source_hashes = payload["production_source_sha256"]
    assert isinstance(source_hashes, dict)
    assert tuple(source_hashes) == (
        "src/phospy/preprocessing/batch_correction/execution/missingness.py",
        "src/phospy/preprocessing/batch_correction/execution/ruv_iii.py",
        "src/phospy/preprocessing/batch_correction/execution/ruv_iii_executor.py",
    )
    assert all(
        isinstance(digest, str) and len(digest) == 64
        for digest in source_hashes.values()
    )
