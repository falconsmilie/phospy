from __future__ import annotations

import numpy as np
import pandas as pd
import pandas.testing as pdt

from phospy.advanced import DatasetBatchCorrectionConfig
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionResult,
    LinearResidualizeBatchCorrectionEngine,
)


def test_linear_residualize_batch_removes_simple_additive_batch_effect() -> None:
    result = _run_engine(_phospho_with_batch_effect())

    expected = pd.DataFrame(
        {
            "sample_1": [10.0, 2.0],
            "sample_2": [14.0, 1.0],
            "sample_3": [10.0, 2.0],
            "sample_4": [14.0, 1.0],
        },
        index=_site_index(),
        columns=_sample_columns(),
    )
    pdt.assert_frame_equal(
        result.corrected_matrix,
        expected,
        check_exact=False,
        atol=1e-10,
        rtol=0.0,
    )


def test_linear_residualize_batch_preserves_condition_effect() -> None:
    phospho = _phospho_with_batch_effect()

    result = _run_engine(phospho)

    pdt.assert_series_equal(
        _condition_effect(result.corrected_matrix),
        _condition_effect(phospho),
        check_exact=False,
        atol=1e-10,
        rtol=0.0,
    )


def test_linear_residualize_batch_reduces_batch_effect() -> None:
    phospho = _phospho_with_batch_effect()

    result = _run_engine(phospho)

    assert _batch_effect_magnitude(result.corrected_matrix) < (
        _batch_effect_magnitude(phospho) * 1e-10
    )


def test_linear_residualize_batch_does_not_mutate_input_matrix() -> None:
    phospho = _phospho_with_batch_effect()
    original = phospho.copy(deep=True)

    result = _run_engine(phospho)

    pdt.assert_frame_equal(phospho, original)
    assert result.corrected_matrix is not phospho


def test_linear_residualize_batch_output_is_deterministic() -> None:
    phospho = _phospho_with_batch_effect()

    first = _run_engine(phospho)
    second = _run_engine(phospho)

    pdt.assert_frame_equal(first.corrected_matrix, second.corrected_matrix)
    assert first.report.to_payload() == second.report.to_payload()
    assert dict(first.diagnostics) == dict(second.diagnostics)


def test_linear_residualize_batch_report_status_is_applied() -> None:
    result = _run_engine(_phospho_with_batch_effect())

    assert isinstance(result, BatchCorrectionResult)
    assert result.report.status == "applied"
    assert result.report.method == "linear_residualize_batch"
    assert result.report.confounding_check_status == "passed"
    assert result.report.batch_levels == ("run_1", "run_2")
    assert result.report.condition_levels == ("control", "treated")
    assert result.report.matrix_shape_before == (2, 4)
    assert result.report.matrix_shape_after == (2, 4)
    assert result.diagnostics["status"] == "applied"


def test_linear_residualize_batch_preserves_matrix_shape_index_and_columns() -> None:
    phospho = _phospho_with_batch_effect()

    result = _run_engine(phospho)

    assert result.corrected_matrix.shape == phospho.shape
    assert result.corrected_matrix.index.equals(phospho.index)
    assert result.corrected_matrix.columns.equals(phospho.columns)
    assert result.corrected_matrix.index.name == phospho.index.name
    assert result.corrected_matrix.columns.name == phospho.columns.name


def _run_engine(phospho: pd.DataFrame) -> BatchCorrectionResult:
    return LinearResidualizeBatchCorrectionEngine().run(
        phospho=phospho,
        batch_labels=("run_1", "run_1", "run_2", "run_2"),
        condition_labels=("control", "treated", "control", "treated"),
        config=DatasetBatchCorrectionConfig(method="linear_residualize_batch"),
    )


def _phospho_with_batch_effect() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [10.0, 2.0],
            "sample_2": [14.0, 1.0],
            "sample_3": [15.0, -1.0],
            "sample_4": [19.0, -2.0],
        },
        index=_site_index(),
        columns=_sample_columns(),
    )


def _site_index() -> pd.Index:
    return pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id")


def _sample_columns() -> pd.Index:
    return pd.Index(
        ["sample_1", "sample_2", "sample_3", "sample_4"],
        name="sample_id",
    )


def _condition_effect(matrix: pd.DataFrame) -> pd.Series:
    treated_mean = matrix.loc[:, ["sample_2", "sample_4"]].mean(axis=1)
    control_mean = matrix.loc[:, ["sample_1", "sample_3"]].mean(axis=1)
    return treated_mean - control_mean


def _batch_effect_magnitude(matrix: pd.DataFrame) -> float:
    control_batch_gap = matrix.loc[:, "sample_3"] - matrix.loc[:, "sample_1"]
    treated_batch_gap = matrix.loc[:, "sample_4"] - matrix.loc[:, "sample_2"]
    batch_gaps = pd.concat([control_batch_gap, treated_batch_gap], axis=1)
    return float(np.abs(batch_gaps.to_numpy()).mean())
