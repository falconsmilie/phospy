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
from phospy.science.differential.linear_model import (
    DIFFERENTIAL_LINEAR_MODEL_DECOMPOSITION_METHOD,
    DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER,
    DIFFERENTIAL_LINEAR_MODEL_SOLVER,
    DifferentialDesignDecompositionError,
    decompose_differential_design,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
    EmpiricalBayesConfig,
    EmpiricalBayesPriorDiagnostics,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_FEATURE_INDEX = pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id")
_STRICT_FEATURE_INDEX = protein_site_key_index(
    protein_identifiers=["MAPK14", "GSK3B"],
    sites=["Y182", "S9"],
)
_DISPLAY_IDS = ["MAPK14;Y182;", "GSK3B;S9;"]


def _series(
    values: list[float],
    *,
    name: str,
    index: pd.Index | None = None,
) -> pd.Series:
    resolved_index = _STRICT_FEATURE_INDEX if index is None else index
    return pd.Series(values, index=resolved_index.copy(), name=name, dtype=float)


def _prior_diagnostics(index: pd.Index | None = None) -> EmpiricalBayesPriorDiagnostics:
    resolved_index = _STRICT_FEATURE_INDEX if index is None else index
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
        prior_variance=_series(
            [1.0, 1.2],
            name="prior_residual_variance",
            index=resolved_index,
        ),
        prior_degrees_of_freedom=_series(
            [5.0, 5.0],
            name="prior_degrees_of_freedom",
            index=resolved_index,
        ),
    )


def _valid_result_table() -> pd.DataFrame:
    site_keys = _STRICT_FEATURE_INDEX.astype(str).tolist()
    return pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": _DISPLAY_IDS,
            **site_key_context_columns(_STRICT_FEATURE_INDEX),
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "logFC": [0.3, -0.1],
            "t": [3.0, -1.2],
            "P.Value": [0.01, 0.24],
            "adj.P.Val": [0.02, 0.24],
        },
        index=_STRICT_FEATURE_INDEX.copy(),
    )


def _build_result(*, table: pd.DataFrame) -> DifferentialAnalysisResult:
    index = table.index.copy()
    return DifferentialAnalysisResult(
        residual_variance=_series(
            [0.1, 0.2],
            name="residual_variance",
            index=index,
        ),
        posterior_residual_variance=_series(
            [0.11, 0.21],
            name="posterior_residual_variance",
            index=index,
        ),
        prior_residual_variance=_series(
            [0.12, 0.22],
            name="prior_residual_variance",
            index=index,
        ),
        prior_degrees_of_freedom_series_value=_series(
            [5.0, 5.0],
            name="prior_degrees_of_freedom",
            index=index,
        ),
        prior_variance=0.17,
        prior_degrees_of_freedom=5.0,
        residual_degrees_of_freedom=2.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=_prior_diagnostics(index),
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


def _manual_contrast_effects(
    *,
    matrix: pd.DataFrame,
    design: pd.DataFrame,
    contrasts: pd.DataFrame,
) -> pd.DataFrame:
    response = matrix.loc[:, design.index].to_numpy(dtype=float).T
    design_values = design.to_numpy(dtype=float)
    contrast_values = contrasts.loc[design.columns].to_numpy(dtype=float)
    coefficients, *_ = np.linalg.lstsq(design_values, response, rcond=None)
    return pd.DataFrame(
        coefficients.T @ contrast_values,
        index=matrix.index.copy(),
        columns=contrasts.columns.copy(),
    )


def _continuous_covariate_inputs() -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    samples = pd.Index(
        ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"],
        name="sample",
    )
    dose = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], dtype=float)
    noise = np.array([0.10, -0.12, 0.06, -0.08, 0.09, -0.04], dtype=float)
    condition_b = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], dtype=float)
    feature_values = {
        "site1": 10.0 + 2.0 * dose + noise,
        "site2": 5.0 + condition_b + 0.5 * dose + noise[::-1],
        "site3": 1.0 + 0.2 * np.arange(dose.size, dtype=float) + noise * 0.3,
    }
    matrix = pd.DataFrame(
        {
            sample: [values[position] for values in feature_values.values()]
            for position, sample in enumerate(samples)
        },
        index=pd.Index(tuple(feature_values), name="site_id"),
    )
    condition_design = pd.DataFrame(
        {
            "A": [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
            "B": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        },
        index=samples.copy(),
        dtype=float,
    )
    adjusted_design = condition_design.copy(deep=True)
    adjusted_design.loc[:, "dose"] = dose
    adjusted_contrast = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0, 0.0]},
        index=pd.Index(["A", "B", "dose"], name="coefficient"),
    )
    return matrix, condition_design, adjusted_design, adjusted_contrast


def _categorical_covariate_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    samples = pd.Index(
        ["A_F1", "A_F2", "A_M1", "B_F1", "B_M1", "B_M2"],
        name="sample",
    )
    sex_m = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 1.0], dtype=float)
    condition_b = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], dtype=float)
    noise = np.array([0.05, -0.04, 0.03, -0.02, 0.04, -0.03], dtype=float)
    feature_values = {
        "site1": 10.0 + condition_b + 4.0 * sex_m + noise,
        "site2": 8.0 - 2.0 * sex_m + noise[::-1],
        "site3": 4.0 + 0.5 * condition_b + 0.7 * sex_m + noise * 2.0,
    }
    matrix = pd.DataFrame(
        {
            sample: [values[position] for values in feature_values.values()]
            for position, sample in enumerate(samples)
        },
        index=pd.Index(tuple(feature_values), name="site_id"),
    )
    design = pd.DataFrame(
        {
            "A": [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
            "B": [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            "sex[M]": sex_m,
        },
        index=samples.copy(),
        dtype=float,
    )
    contrast = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0, 0.0]},
        index=pd.Index(["A", "B", "sex[M]"], name="coefficient"),
    )
    return matrix, design, contrast


def test_condition_only_executor_preserves_group_mean_contrast() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 0.8],
            "A_2": [1.2, 2.1, 0.9],
            "B_1": [1.8, 1.9, 1.4],
            "B_2": [2.0, 2.2, 1.6],
        },
        index=pd.Index(["site1", "site2", "site3"], name="site_id"),
    )
    request = _base_request(matrix=matrix, empirical_bayes=EmpiricalBayesConfig())

    result = DifferentialAnalysisExecutor().run(request)
    table = result.table_for("B_vs_A")

    expected = matrix.loc[:, ["B_1", "B_2"]].mean(axis=1) - matrix.loc[
        :, ["A_1", "A_2"]
    ].mean(axis=1)
    np.testing.assert_allclose(
        table.loc[:, "logFC"].to_numpy(dtype=float),
        expected.to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    assert result.residual_degrees_of_freedom == pytest.approx(2.0)


def test_adjusted_contrast_differs_when_covariate_explains_signal() -> None:
    matrix, condition_design, adjusted_design, adjusted_contrast = (
        _continuous_covariate_inputs()
    )
    condition_contrast = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0]},
        index=pd.Index(["A", "B"], name="coefficient"),
    )

    unadjusted = DifferentialAnalysisExecutor().run(
        DifferentialAnalysisRequest(
            matrix=matrix,
            design=condition_design,
            contrasts=condition_contrast,
        )
    )
    adjusted = DifferentialAnalysisExecutor().run(
        DifferentialAnalysisRequest(
            matrix=matrix,
            design=adjusted_design,
            contrasts=adjusted_contrast,
        )
    )

    unadjusted_log_fc = float(unadjusted.table_for("B_vs_A").at["site1", "logFC"])
    adjusted_log_fc = float(adjusted.table_for("B_vs_A").at["site1", "logFC"])
    assert abs(unadjusted_log_fc) > 5.0
    assert abs(adjusted_log_fc) < 0.1
    assert abs(unadjusted_log_fc - adjusted_log_fc) > 5.0


def test_categorical_covariate_adjusted_model_uses_supplied_design() -> None:
    matrix, design, contrast = _categorical_covariate_inputs()

    result = DifferentialAnalysisExecutor().run(
        DifferentialAnalysisRequest(matrix=matrix, design=design, contrasts=contrast)
    )
    expected = _manual_contrast_effects(
        matrix=matrix,
        design=design,
        contrasts=contrast,
    )

    np.testing.assert_allclose(
        result.table_for("B_vs_A").loc[:, "logFC"].to_numpy(dtype=float),
        expected.loc[:, "B_vs_A"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    assert result.table_for("B_vs_A").at["site1", "logFC"] == pytest.approx(0.975)


def test_continuous_covariate_adjusted_model_uses_supplied_design() -> None:
    matrix, _, adjusted_design, adjusted_contrast = _continuous_covariate_inputs()

    result = DifferentialAnalysisExecutor().run(
        DifferentialAnalysisRequest(
            matrix=matrix,
            design=adjusted_design,
            contrasts=adjusted_contrast,
        )
    )
    expected = _manual_contrast_effects(
        matrix=matrix,
        design=adjusted_design,
        contrasts=adjusted_contrast,
    )

    np.testing.assert_allclose(
        result.table_for("B_vs_A").loc[:, "logFC"].to_numpy(dtype=float),
        expected.loc[:, "B_vs_A"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    assert result.residual_degrees_of_freedom == pytest.approx(3.0)


def test_full_rank_near_confounded_covariate_fit_is_stable() -> None:
    samples = pd.Index(
        ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"],
        name="sample",
    )
    condition_b = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], dtype=float)
    almost_condition_b = condition_b + np.array(
        [-1e-5, 0.0, 1e-5, -1e-5, 0.0, 1e-5],
        dtype=float,
    )
    design = pd.DataFrame(
        {
            "A": 1.0 - condition_b,
            "B": condition_b,
            "almost_B": almost_condition_b,
        },
        index=samples,
        dtype=float,
    )
    contrasts = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0, 0.0]},
        index=pd.Index(["A", "B", "almost_B"], name="coefficient"),
    )
    matrix = pd.DataFrame(
        {
            sample: [
                2.0
                + 0.4 * condition_b[position]
                + 0.7 * almost_condition_b[position]
                + 0.01 * float(position),
                -1.0
                + 0.2 * condition_b[position]
                - 0.3 * almost_condition_b[position]
                + 0.02 * float(position),
            ]
            for position, sample in enumerate(samples)
        },
        index=pd.Index(["site1", "site2"], name="site_id"),
    )

    decomposition = decompose_differential_design(design.to_numpy(dtype=float))
    assert decomposition.rank == 3
    assert (
        decomposition.condition_number < DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER
    )

    result = DifferentialAnalysisExecutor().run(
        DifferentialAnalysisRequest(matrix=matrix, design=design, contrasts=contrasts)
    )
    table = result.table_for("B_vs_A")
    expected = _manual_contrast_effects(
        matrix=matrix,
        design=design,
        contrasts=contrasts,
    )
    np.testing.assert_allclose(
        table.loc[:, "logFC"].to_numpy(dtype=float),
        expected.loc[:, "B_vs_A"].to_numpy(dtype=float),
        rtol=1e-8,
        atol=1e-8,
    )
    assert np.isfinite(table.loc[:, "t"]).all()


def test_executor_rejects_exactly_rank_deficient_design() -> None:
    matrix, _, adjusted_design, adjusted_contrast = _continuous_covariate_inputs()
    deficient_design = adjusted_design.copy(deep=True)
    deficient_design.loc[:, "duplicate_B"] = deficient_design.loc[:, "B"]
    deficient_contrast = adjusted_contrast.copy(deep=True)
    deficient_contrast.loc["duplicate_B"] = 0.0

    with pytest.raises(
        PhosPyInputError,
        match="rank deficient under scaled SVD decomposition",
    ):
        DifferentialAnalysisExecutor().run(
            DifferentialAnalysisRequest(
                matrix=matrix,
                design=deficient_design,
                contrasts=deficient_contrast,
            )
        )


def test_rescaled_covariate_preserves_estimable_condition_contrast() -> None:
    matrix, _, design, contrast = _continuous_covariate_inputs()
    rescaled_design = design.copy(deep=True)
    rescaled_design.loc[:, "dose"] = rescaled_design.loc[:, "dose"] * 1.0e6

    original = DifferentialAnalysisExecutor().run(
        DifferentialAnalysisRequest(matrix=matrix, design=design, contrasts=contrast)
    )
    rescaled = DifferentialAnalysisExecutor().run(
        DifferentialAnalysisRequest(
            matrix=matrix,
            design=rescaled_design,
            contrasts=contrast,
        )
    )

    pd.testing.assert_series_equal(
        original.table_for("B_vs_A").loc[:, "logFC"],
        rescaled.table_for("B_vs_A").loc[:, "logFC"],
        rtol=1e-9,
        atol=1e-9,
        check_names=False,
    )
    pd.testing.assert_series_equal(
        original.table_for("B_vs_A").loc[:, "t"],
        rescaled.table_for("B_vs_A").loc[:, "t"],
        rtol=1e-7,
        atol=1e-7,
        check_names=False,
    )


def test_scaled_svd_decomposition_matches_known_coefficient_covariance_fixture() -> (
    None
):
    design = np.array(
        [
            [1.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 1.0],
        ],
        dtype=float,
    )
    response = np.array([[1.0], [3.0], [3.0], [5.0]], dtype=float)

    decomposition = decompose_differential_design(design)
    fit = decomposition.fit(response)
    covariance = decomposition.coefficient_covariance
    contrast_covariance = decomposition.contrast_covariance(
        np.array([[-1.0], [1.0]], dtype=float)
    )

    np.testing.assert_allclose(fit.coefficients[:, 0], np.array([2.0, 4.0]))
    assert fit.residual_variance[0] == pytest.approx(2.0)
    np.testing.assert_allclose(covariance, np.diag([0.5, 0.5]), atol=1e-12)
    assert contrast_covariance[0, 0] == pytest.approx(1.0)
    assert (
        decomposition.decomposition_method
        == DIFFERENTIAL_LINEAR_MODEL_DECOMPOSITION_METHOD
    )
    assert decomposition.solver == DIFFERENTIAL_LINEAR_MODEL_SOLVER


def test_contrast_covariance_regression_for_adjusted_design() -> None:
    design = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=float,
    )
    contrasts = np.array(
        [
            [-1.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        dtype=float,
    )

    covariance = decompose_differential_design(design).contrast_covariance(contrasts)

    np.testing.assert_allclose(
        covariance,
        np.array([[1.0, 0.0], [0.0, 1.0]], dtype=float),
        atol=1e-12,
    )


def test_scaled_svd_condition_threshold_boundary() -> None:
    def _design_for_condition(target_condition: float) -> np.ndarray:
        rho = (target_condition**2 - 1.0) / (target_condition**2 + 1.0)
        return np.array(
            [
                [1.0, rho],
                [0.0, np.sqrt(1.0 - rho**2)],
                [0.0, 0.0],
            ],
            dtype=float,
        )

    at_threshold = decompose_differential_design(
        _design_for_condition(10.0),
        max_condition_number=10.0,
    )
    assert at_threshold.condition_number == pytest.approx(10.0)

    with pytest.raises(
        DifferentialDesignDecompositionError,
        match="too ill-conditioned",
    ):
        decompose_differential_design(
            _design_for_condition(10.1),
            max_condition_number=10.0,
        )


def test_executor_rejects_misaligned_contrast_vector_clearly() -> None:
    matrix, _, adjusted_design, _ = _continuous_covariate_inputs()
    invalid_contrast = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0, 0.0]},
        index=pd.Index(["A", "B", "unknown"], name="coefficient"),
    )

    with pytest.raises(
        PhosPyInputError,
        match=r"differential\.contrasts\.index must match "
        r"differential\.design\.columns",
    ):
        DifferentialAnalysisExecutor().run(
            DifferentialAnalysisRequest(
                matrix=matrix,
                design=adjusted_design,
                contrasts=invalid_contrast,
            )
        )


def test_empirical_bayes_runs_with_adjusted_design() -> None:
    matrix, _, adjusted_design, adjusted_contrast = _continuous_covariate_inputs()

    result = DifferentialAnalysisExecutor().run(
        DifferentialAnalysisRequest(
            matrix=matrix,
            design=adjusted_design,
            contrasts=adjusted_contrast,
            empirical_bayes=EmpiricalBayesConfig(method="robust", trend=True),
        )
    )

    assert result.empirical_bayes_method == "robust"
    assert result.empirical_bayes_robust is True
    assert result.empirical_bayes_trend is True
    assert result.mean_variance_trend_diagnostics is not None
    assert np.isfinite(result.table_for("B_vs_A").loc[:, "P.Value"]).all()


def test_request_rejects_all_constant_site_intensities() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [5.0, 1.0],
            "A_2": [5.0, 1.0 + 1e-12],
            "B_1": [5.0, 1.0 + 2e-12],
            "B_2": [5.0, 1.0 + 3e-12],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;"],
            name="site_id",
        ),
    )

    with pytest.raises(PhosPyInputError, match="all-constant site intensities"):
        _base_request(matrix=matrix, empirical_bayes=EmpiricalBayesConfig())


def test_executor_supports_near_zero_variance_rows() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 3.0],
            "A_2": [1.0 + 1e-12, 3.1],
            "B_1": [1.0 + 2e-12, 4.0],
            "B_2": [1.0 + 3e-12, 3.9],
        },
        index=pd.Index(
            ["GSK3B;S9;", "AKT1;T308;"],
            name="site_id",
        ),
    )
    result = DifferentialAnalysisExecutor().run(
        _base_request(matrix=matrix, empirical_bayes=EmpiricalBayesConfig())
    )
    table = result.table_for("B_vs_A")

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


def test_differential_rejects_or_documents_non_finite_p_values_before_correction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _NonFiniteTDistribution:
        @staticmethod
        def sf(values: np.ndarray, *, df: np.ndarray) -> np.ndarray:
            del df
            p_values = np.full(np.asarray(values, dtype=float).shape, 0.05)
            p_values[0] = np.nan
            return p_values

    def _unexpected_correction_call(_: np.ndarray, *, method: str) -> np.ndarray:
        del method
        raise AssertionError(
            "p-value correction must not run after non-finite P.Value generation"
        )

    monkeypatch.setattr(
        "phospy.science.differential.executor.stats.t",
        _NonFiniteTDistribution(),
    )
    monkeypatch.setattr(
        "phospy.science.differential.executor.adjust_p_values",
        _unexpected_correction_call,
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
        match=r"P\.Value must be finite and within \[0, 1\]",
    ):
        DifferentialAnalysisExecutor().run(
            _base_request(matrix=matrix, empirical_bayes=EmpiricalBayesConfig())
        )
