from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.science.differential.empirical_bayes import (
    EmpiricalBayesFit,
    fit_empirical_bayes,
)
from phospy.science.differential.executor import DifferentialAnalysisExecutor
from phospy.science.differential.models import (
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
    EmpiricalBayesConfig,
    EmpiricalBayesPriorDiagnostics,
)

_FEATURE_INDEX = pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id")


def _series(values: list[float], *, name: str) -> pd.Series:
    return pd.Series(values, index=_FEATURE_INDEX.copy(), name=name, dtype=float)


def _prior_diagnostics() -> EmpiricalBayesPriorDiagnostics:
    return EmpiricalBayesPriorDiagnostics(
        method="standard",
        robust=False,
        trend=False,
        winsor_tail_p=(0.05, 0.1),
        base_prior_variance=1.0,
        base_prior_degrees_of_freedom=5.0,
        robust_outlier_count=0,
        robust_outlier_fraction=0.0,
        winsorized_low_count=0,
        winsorized_high_count=0,
        prior_variance=_series([1.0, 1.2], name="prior_residual_variance"),
        prior_degrees_of_freedom=_series([5.0, 5.0], name="prior_degrees_of_freedom"),
    )


def _valid_result_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "logFC": [0.3, -0.1],
            "t": [3.0, -1.2],
            "P.Value": [0.01, 0.24],
            "adj.P.Val": [0.02, 0.24],
        },
        index=_FEATURE_INDEX.copy(),
    )


def _build_result(*, table: pd.DataFrame) -> DifferentialAnalysisResult:
    return DifferentialAnalysisResult(
        residual_variance=_series([0.1, 0.2], name="residual_variance"),
        posterior_residual_variance=_series(
            [0.11, 0.21], name="posterior_residual_variance"
        ),
        prior_residual_variance=_series([0.12, 0.22], name="prior_residual_variance"),
        prior_degrees_of_freedom_series_value=_series(
            [5.0, 5.0], name="prior_degrees_of_freedom"
        ),
        prior_variance=0.17,
        prior_degrees_of_freedom=5.0,
        residual_degrees_of_freedom=2.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=_prior_diagnostics(),
        mean_variance_trend_diagnostics=None,
        contrast_tables={"B_vs_A": table},
    )


@pytest.mark.parametrize(
    ("column", "value"),
    (
        ("P.Value", -1e-6),
        ("P.Value", 1.000001),
        ("adj.P.Val", -1e-6),
        ("adj.P.Val", 1.000001),
    ),
)
def test_result_model_rejects_out_of_range_probability_columns(
    column: str, value: float
) -> None:
    table = _valid_result_table()
    table.iloc[0, table.columns.get_loc(column)] = value

    with pytest.raises(PhosPyInputError) as error:
        _build_result(table=table)

    message = str(error.value)
    assert f".{column} must be within [0, 1]" in message
    assert "invalid values" in message


@pytest.mark.parametrize("column", ("P.Value", "adj.P.Val"))
@pytest.mark.parametrize("value", (float("nan"), float("inf"), float("-inf")))
def test_result_model_rejects_non_finite_probability_values(
    column: str, value: float
) -> None:
    table = _valid_result_table()
    table.iloc[1, table.columns.get_loc(column)] = value

    with pytest.raises(PhosPyInputError) as error:
        _build_result(table=table)

    assert column in str(error.value)


def _base_request(*, matrix: pd.DataFrame, empirical_bayes: EmpiricalBayesConfig):
    design = pd.DataFrame(
        {
            "A": [1.0, 1.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0],
        },
        index=pd.Index(["A_1", "A_2", "B_1", "B_2"], name="sample"),
    )
    contrasts = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0]},
        index=pd.Index(["A", "B"], name="coefficient"),
    )
    return DifferentialAnalysisRequest(
        matrix=matrix,
        design=design,
        contrasts=contrasts,
        empirical_bayes=empirical_bayes,
    )


def test_executor_supports_zero_and_near_zero_variance_rows() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [5.0, 1.0, 3.0],
            "A_2": [5.0, 1.0 + 1e-12, 3.1],
            "B_1": [5.0, 1.0 + 2e-12, 4.0],
            "B_2": [5.0, 1.0 + 3e-12, 3.9],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
            name="site_id",
        ),
    )
    result = DifferentialAnalysisExecutor().run(
        _base_request(matrix=matrix, empirical_bayes=EmpiricalBayesConfig())
    )
    table = result.table_for("B_vs_A")

    assert table.at["MAPK14;Y182;", "logFC"] == pytest.approx(0.0)
    assert table.at["MAPK14;Y182;", "P.Value"] == pytest.approx(1.0)
    assert np.isfinite(table.loc[:, "t"]).all()
    assert np.isfinite(table.loc[:, "P.Value"]).all()
    assert (table.loc[:, "P.Value"] >= 0.0).all()
    assert (table.loc[:, "P.Value"] <= 1.0).all()


def test_executor_remains_stable_for_tiny_group_sizes() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 0.8, 1.2],
            "A_2": [1.1, 2.2, 0.9, 1.3],
            "B_1": [1.7, 2.1, 1.1, 1.9],
            "B_2": [1.8, 2.3, 1.0, 2.0],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;", "RPS6KB1;T389;"],
            name="site_id",
        ),
    )
    result = DifferentialAnalysisExecutor().run(
        _base_request(
            matrix=matrix,
            empirical_bayes=EmpiricalBayesConfig(
                method="robust",
                trend=True,
                winsor_tail_p=(0.05, 0.1),
            ),
        )
    )
    table = result.table_for("B_vs_A")

    assert np.isfinite(table.loc[:, "t"]).all()
    assert np.isfinite(table.loc[:, "P.Value"]).all()
    assert (table.loc[:, "P.Value"] >= 0.0).all()
    assert (table.loc[:, "P.Value"] <= 1.0).all()


def test_fit_empirical_bayes_rejects_invalid_residual_degrees_of_freedom() -> None:
    with pytest.raises(
        ValueError,
        match="residual_dof must be finite and > 0.0",
    ):
        fit_empirical_bayes(
            variances=np.array([0.2, 0.4], dtype=float),
            residual_dof=0.0,
            method="standard",
            trend=False,
            winsor_tail_p=(0.05, 0.1),
        )


def test_executor_rejects_invalid_moderated_degrees_of_freedom(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _bad_fit(**_: object) -> EmpiricalBayesFit:
        return EmpiricalBayesFit(
            prior_variance=np.array([0.5, 0.6], dtype=float),
            prior_degrees_of_freedom=np.array([-5.0, -5.0], dtype=float),
            base_prior_variance=0.5,
            base_prior_degrees_of_freedom=-5.0,
            robust_outlier_count=0,
            robust_outlier_fraction=0.0,
            winsorized_low_count=0,
            winsorized_high_count=0,
        )

    monkeypatch.setattr(
        "phospy.science.differential.executor.fit_empirical_bayes", _bad_fit
    )
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.1],
            "B_1": [1.6, 2.3],
            "B_2": [1.7, 2.4],
        },
        index=_FEATURE_INDEX.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="invalid moderated degrees of freedom",
    ):
        DifferentialAnalysisExecutor().run(
            _base_request(matrix=matrix, empirical_bayes=EmpiricalBayesConfig())
        )


def test_executor_rejects_unstable_standard_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _zero_variance_fit(**_: object) -> EmpiricalBayesFit:
        return EmpiricalBayesFit(
            prior_variance=np.array([0.0, 0.0], dtype=float),
            prior_degrees_of_freedom=np.array(
                [float("inf"), float("inf")], dtype=float
            ),
            base_prior_variance=0.0,
            base_prior_degrees_of_freedom=float("inf"),
            robust_outlier_count=0,
            robust_outlier_fraction=0.0,
            winsorized_low_count=0,
            winsorized_high_count=0,
        )

    monkeypatch.setattr(
        "phospy.science.differential.executor.fit_empirical_bayes",
        _zero_variance_fit,
    )
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.0, 2.1],
            "B_1": [1.5, 2.2],
            "B_2": [1.5, 2.3],
        },
        index=_FEATURE_INDEX.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="unstable standard errors",
    ):
        DifferentialAnalysisExecutor().run(
            _base_request(matrix=matrix, empirical_bayes=EmpiricalBayesConfig())
        )
