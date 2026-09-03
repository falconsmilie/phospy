from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest

import phospy.science.differential.protein_covariate_adjusted as kernel_module
from phospy.errors import PhosPyInputError
from phospy.science.differential.empirical_bayes import fit_empirical_bayes
from phospy.science.differential.models import EmpiricalBayesConfig
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationRequest,
)
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_POSITIVE,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID,
)
from phospy.science.differential.protein_covariate_adjusted import (
    PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
    ProteinCovariateAdjustedDifferentialKernel,
)
from phospy.science.statistics.multiple_testing import (
    MULTIPLE_TESTING_CORRECTION_BONFERRONI,
    adjust_p_values,
)

SAMPLES = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
CONDITION_A = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=float)
CONDITION_B = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], dtype=float)
NOISE_SEED = np.array([0.8, -1.3, 0.4, -0.7, 1.1, -0.2], dtype=float)


def test_one_site_one_protein_matches_independent_oracle() -> None:
    protein = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    values = _site_values(beta_a=5.0, beta_b=7.0, gamma=0.4, protein=protein)
    request = _request(
        matrix=_matrix({"site_a": values}),
        proteins=_proteins({"protein_a": protein}),
        pairs=_pairs({"site_a": ("MAPK14", "protein_a")}),
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)
    oracle = _oracle(values=values, protein=protein)

    assert result.tested_site_ids == ("site_a",)
    assert result.prior_degrees_of_freedom == pytest.approx(0.0)
    assert result.augmented_design_diagnostics_dataframe().loc[
        "protein_a",
        "protein_covariate_contrast_weights",
    ] == (0.0,)
    np.testing.assert_allclose(
        result.coefficient_dataframe().loc["site_a"].to_numpy(dtype=float),
        oracle["coefficients"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        result.residuals_dataframe().loc["site_a"].to_numpy(dtype=float),
        oracle["residuals"],
        rtol=1e-12,
        atol=1e-12,
    )
    assert result.residual_variance_series().loc["site_a"] == pytest.approx(
        float(oracle["residual_variance"]),
        rel=1e-12,
        abs=1e-12,
    )
    assert result.table_for("B_vs_A").loc["site_a", "logFC"] == pytest.approx(
        float(oracle["contrast_effects"][0]),
        rel=1e-12,
        abs=1e-12,
    )
    assert result.contrast_standard_error_scale_dataframe().loc[
        "site_a",
        "B_vs_A",
    ] == pytest.approx(float(oracle["contrast_scales"][0]), rel=1e-12, abs=1e-12)


def test_grouped_fit_restores_site_order_and_decomposes_once_per_protein(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    protein_a = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    protein_b = np.array([2.0, 1.1, 0.3, 1.7, 0.4, 2.2], dtype=float)
    values = {
        "site_b": _site_values(
            beta_a=3.0,
            beta_b=2.7,
            gamma=-0.3,
            protein=protein_b,
            seed=NOISE_SEED[::-1],
        ),
        "site_a": _site_values(beta_a=5.0, beta_b=7.0, gamma=0.4, protein=protein_a),
        "site_c": _site_values(
            beta_a=1.0,
            beta_b=1.6,
            gamma=1.2,
            protein=protein_a,
            seed=np.roll(NOISE_SEED, 1),
        ),
    }
    calls: list[np.ndarray] = []
    real_decompose = kernel_module.decompose_differential_design

    def _counting_decompose(
        design: np.ndarray,
        *,
        max_condition_number: float = 1.0e10,
    ):
        calls.append(np.asarray(design, dtype=float).copy())
        return real_decompose(
            design,
            max_condition_number=max_condition_number,
        )

    monkeypatch.setattr(
        kernel_module,
        "decompose_differential_design",
        _counting_decompose,
    )
    request = _request(
        matrix=_matrix(values),
        proteins=_proteins({"protein_a": protein_a, "protein_b": protein_b}),
        pairs=_pairs(
            {
                "site_b": ("GSK3B", "protein_b"),
                "site_a": ("MAPK14", "protein_a"),
                "site_c": ("AKT1", "protein_a"),
            }
        ),
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert result.tested_site_ids == ("site_b", "site_a", "site_c")
    assert len(calls) == 2
    assert result.augmented_design_diagnostics_dataframe().index.tolist() == [
        "protein_b",
        "protein_a",
    ]
    for site_id, protein in (
        ("site_b", protein_b),
        ("site_a", protein_a),
        ("site_c", protein_a),
    ):
        oracle = _oracle(values=values[site_id], protein=protein)
        np.testing.assert_allclose(
            result.coefficient_dataframe().loc[site_id].to_numpy(dtype=float),
            oracle["coefficients"],
            rtol=1e-12,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            result.contrast_standard_error_scale_dataframe()
            .loc[site_id]
            .to_numpy(dtype=float),
            oracle["contrast_scales"],
            rtol=1e-12,
            atol=1e-12,
        )
    scales = result.contrast_standard_error_scale_dataframe()["B_vs_A"]
    assert scales.loc["site_a"] == pytest.approx(scales.loc["site_c"])
    assert not math.isclose(float(scales.loc["site_a"]), float(scales.loc["site_b"]))


def test_ordered_grouping_reuses_one_vectorized_fit_per_protein_row(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    row_key_sequence = (
        "protein_b",
        "protein_a",
        "protein_b",
        "protein_c",
        "protein_a",
    )
    groups = kernel_module._positions_by_total_protein_row_key(row_key_sequence)

    assert list(groups) == ["protein_b", "protein_a", "protein_c"]
    assert {row_key: list(positions) for row_key, positions in groups.items()} == {
        "protein_b": [0, 2],
        "protein_a": [1, 4],
        "protein_c": [3],
    }

    protein_a = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    protein_b = np.array([2.0, 1.1, 0.3, 1.7, 0.4, 2.2], dtype=float)
    protein_c = np.array([1.4, 0.2, 2.1, 0.7, 1.8, 2.6], dtype=float)
    values = {
        "site_0": _site_values(
            beta_a=4.0,
            beta_b=4.8,
            gamma=0.3,
            protein=protein_b,
            seed=NOISE_SEED[::-1],
        ),
        "site_1": _site_values(
            beta_a=5.0,
            beta_b=6.0,
            gamma=-0.2,
            protein=protein_a,
        ),
        "site_2": _site_values(
            beta_a=3.0,
            beta_b=4.1,
            gamma=0.7,
            protein=protein_b,
            seed=np.roll(NOISE_SEED, 1),
        ),
        "site_3": _site_values(
            beta_a=2.5,
            beta_b=3.4,
            gamma=-0.4,
            protein=protein_c,
            seed=np.roll(NOISE_SEED, 2),
        ),
        "site_4": _site_values(
            beta_a=6.2,
            beta_b=6.9,
            gamma=0.5,
            protein=protein_a,
            seed=np.roll(NOISE_SEED, 3),
        ),
    }
    decompose_calls: list[np.ndarray] = []
    fit_calls: list[tuple[int, ...]] = []
    real_decompose = kernel_module.decompose_differential_design
    real_fit_group = kernel_module._fit_group

    def _counting_decompose(
        design: np.ndarray,
        *,
        max_condition_number: float = 1.0e10,
    ):
        decompose_calls.append(np.asarray(design, dtype=float).copy())
        return real_decompose(
            design,
            max_condition_number=max_condition_number,
        )

    def _counting_fit_group(**kwargs: Any) -> None:
        fit_calls.append(tuple(int(value) for value in kwargs["positions"]))
        real_fit_group(**kwargs)

    monkeypatch.setattr(
        kernel_module,
        "decompose_differential_design",
        _counting_decompose,
    )
    monkeypatch.setattr(kernel_module, "_fit_group", _counting_fit_group)

    request = _request(
        matrix=_matrix(values),
        proteins=_proteins(
            {
                "protein_b": protein_b,
                "protein_a": protein_a,
                "protein_c": protein_c,
            }
        ),
        pairs=_pairs(
            {
                "site_0": ("P_B", "protein_b"),
                "site_1": ("P_A", "protein_a"),
                "site_2": ("P_B", "protein_b"),
                "site_3": ("P_C", "protein_c"),
                "site_4": ("P_A", "protein_a"),
            }
        ),
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert len(decompose_calls) == 3
    assert fit_calls == [(0, 2), (1, 4), (3,)]
    assert result.tested_site_ids == (
        "site_0",
        "site_1",
        "site_2",
        "site_3",
        "site_4",
    )
    assert result.augmented_design_diagnostics_dataframe().index.tolist() == [
        "protein_b",
        "protein_a",
        "protein_c",
    ]
    assert result.site_diagnostics_dataframe()["protein_identifier"].tolist() == [
        "P_B",
        "P_A",
        "P_B",
        "P_C",
        "P_A",
    ]


def test_invalid_protein_covariate_and_augmented_design_groups_are_typed() -> None:
    valid = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    constant = np.full(6, 4.2, dtype=float)
    non_finite = np.array([0.2, 0.4, np.nan, 0.6, 0.8, 1.0], dtype=float)
    collinear = CONDITION_B.copy()
    ill_conditioned = CONDITION_B + np.array(
        [-1.0e-10, 0.0, 1.0e-10, -1.0e-10, 0.0, 1.0e-10],
        dtype=float,
    )
    request = _request(
        matrix=_matrix(
            {
                "valid_site": _site_values(5.0, 6.0, 0.25, valid),
                "constant_site": _site_values(2.0, 2.5, 0.1, valid),
                "nonfinite_site": _site_values(2.0, 2.7, -0.1, valid),
                "collinear_site": _site_values(3.0, 3.7, 0.2, valid),
                "ill_conditioned_site": _site_values(4.0, 4.9, 0.3, valid),
            }
        ),
        proteins=_proteins(
            {
                "valid_protein": valid,
                "constant_protein": constant,
                "nonfinite_protein": non_finite,
                "collinear_protein": collinear,
                "ill_conditioned_protein": ill_conditioned,
            }
        ),
        pairs=_pairs(
            {
                "valid_site": ("MAPK14", "valid_protein"),
                "constant_site": ("AKT1", "constant_protein"),
                "nonfinite_site": ("GSK3B", "nonfinite_protein"),
                "collinear_site": ("RPS6", "collinear_protein"),
                "ill_conditioned_site": ("EIF4EBP1", "ill_conditioned_protein"),
            }
        ),
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert result.tested_site_ids == ("valid_site",)
    failures = result.site_failure_diagnostics_dataframe()
    assert failures.index.tolist() == [
        "constant_site",
        "nonfinite_site",
        "collinear_site",
        "ill_conditioned_site",
    ]
    assert failures.loc["constant_site", DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID
    )
    assert failures.loc["constant_site", DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE
    )
    assert failures.loc["nonfinite_site", DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE
    )
    group_failures = result.augmented_design_failure_diagnostics_dataframe()
    assert (
        group_failures.loc[
            "collinear_protein",
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        ]
        == DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
    )
    assert (
        group_failures.loc[
            "ill_conditioned_protein",
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        ]
        == DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED
    )
    assert (
        group_failures.loc[
            "ill_conditioned_protein",
            DIFFERENTIAL_RESULT_STATUS_COLUMN,
        ]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
    )


def test_non_estimable_augmented_contrast_is_typed_in_failure_report() -> None:
    protein = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    zero_contrast = pd.DataFrame(
        {"zero": [0.0, 0.0]},
        index=pd.Index(["A", "B"], name="coefficient"),
    )
    request = _request(
        matrix=_matrix({"site_a": _site_values(5.0, 6.0, 0.25, protein)}),
        proteins=_proteins({"protein_a": protein}),
        pairs=_pairs({"site_a": ("MAPK14", "protein_a")}),
        contrasts=zero_contrast,
    )

    with pytest.raises(PhosPyInputError, match="no successfully tested sites") as exc:
        ProteinCovariateAdjustedDifferentialKernel().run(request)

    report = exc.value.diagnostics
    assert isinstance(report, dict)
    failures = report["site_failure_diagnostics"]
    assert isinstance(failures, pd.DataFrame)
    assert failures.loc["site_a", DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE
    )
    assert failures.loc["site_a", DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE
    )


def test_non_positive_augmented_residual_df_is_typed_in_failure_report() -> None:
    samples = ("A_1", "A_2", "B_1")
    design = pd.DataFrame(
        {"A": [1.0, 1.0, 0.0], "B": [0.0, 0.0, 1.0]},
        index=pd.Index(samples, name="sample_id"),
    )
    contrasts = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0]},
        index=pd.Index(["A", "B"], name="coefficient"),
    )
    protein = np.array([0.0, 1.0, 0.0], dtype=float)
    matrix = pd.DataFrame(
        {"A_1": [2.0], "A_2": [2.3], "B_1": [3.1]},
        index=pd.Index(["site_a"], name="site_key"),
    )
    request = _request(
        matrix=matrix,
        proteins=_proteins({"protein_a": protein}, samples=samples),
        pairs=_pairs({"site_a": ("MAPK14", "protein_a")}),
        design=design,
        contrasts=contrasts,
        sample_order=samples,
    )

    with pytest.raises(PhosPyInputError, match="no successfully tested sites") as exc:
        ProteinCovariateAdjustedDifferentialKernel().run(request)

    report = exc.value.diagnostics
    assert isinstance(report, dict)
    group_failures = report["augmented_design_failure_diagnostics"]
    assert isinstance(group_failures, pd.DataFrame)
    assert (
        group_failures.loc[
            "protein_a",
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        ]
        == DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF
    )


@pytest.mark.parametrize(
    "empirical_bayes",
    [
        EmpiricalBayesConfig(method="standard", trend=False),
        EmpiricalBayesConfig(method="robust", trend=False),
        EmpiricalBayesConfig(method="standard", trend=True),
        EmpiricalBayesConfig(method="robust", trend=True),
    ],
)
def test_empirical_bayes_modes_match_existing_direct_helper(
    empirical_bayes: EmpiricalBayesConfig,
) -> None:
    protein_a = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    protein_b = np.array([2.0, 1.1, 0.3, 1.7, 0.4, 2.2], dtype=float)
    request = _request(
        matrix=_matrix(
            {
                "site_a": _site_values(5.0, 6.0, 0.25, protein_a),
                "site_b": _site_values(4.0, 4.3, -0.2, protein_a),
                "site_c": _site_values(1.5, 2.1, 0.8, protein_b),
                "site_d": _site_values(3.0, 2.2, -0.5, protein_b),
            }
        ),
        proteins=_proteins({"protein_a": protein_a, "protein_b": protein_b}),
        pairs=_pairs(
            {
                "site_a": ("MAPK14", "protein_a"),
                "site_b": ("AKT1", "protein_a"),
                "site_c": ("GSK3B", "protein_b"),
                "site_d": ("RPS6", "protein_b"),
            }
        ),
        empirical_bayes=empirical_bayes,
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)
    expected = fit_empirical_bayes(
        variances=result.residual_variance_series().to_numpy(dtype=float),
        residual_dof=result.residual_degrees_of_freedom,
        method=empirical_bayes.method,
        trend=empirical_bayes.trend,
        winsor_tail_p=empirical_bayes.winsor_tail_p,
        mean_intensity=result.site_diagnostics_dataframe()
        .loc[:, "mean_intensity"]
        .to_numpy(dtype=float),
    )

    assert result.empirical_bayes_method == empirical_bayes.method
    assert result.empirical_bayes_robust is (empirical_bayes.method == "robust")
    assert result.empirical_bayes_trend is empirical_bayes.trend
    assert (result.mean_variance_trend_diagnostics is not None) is (
        empirical_bayes.trend
    )
    np.testing.assert_allclose(
        result.prior_residual_variance_series().to_numpy(dtype=float),
        expected.prior_variance,
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        result.prior_degrees_of_freedom_series().to_numpy(dtype=float),
        expected.prior_degrees_of_freedom,
        rtol=1e-12,
        atol=1e-12,
    )


def test_multiple_testing_uses_successfully_tested_sites_only_per_contrast() -> None:
    protein_a = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    protein_b = np.array([2.0, 1.1, 0.3, 1.7, 0.4, 2.2], dtype=float)
    constant = np.full(6, 3.0, dtype=float)
    request = _request(
        matrix=_matrix(
            {
                "site_a": _site_values(5.0, 6.0, 0.2, protein_a),
                "site_b": _site_values(4.0, 4.8, -0.3, protein_b),
                "site_c": _site_values(2.0, 2.4, 0.1, protein_a),
            }
        ),
        proteins=_proteins(
            {
                "protein_a": protein_a,
                "protein_b": protein_b,
                "constant_protein": constant,
            }
        ),
        pairs=_pairs(
            {
                "site_a": ("MAPK14", "protein_a"),
                "site_b": ("GSK3B", "protein_b"),
                "site_c": ("AKT1", "constant_protein"),
            }
        ),
        contrasts=_contrasts(multiple=True),
        multiple_testing_method=MULTIPLE_TESTING_CORRECTION_BONFERRONI,
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert result.tested_site_ids == ("site_a", "site_b")
    for contrast_name in ("B_vs_A", "A_level"):
        table = result.table_for(contrast_name)
        p_values = table["P.Value"].to_numpy(dtype=float)
        expected = adjust_p_values(
            p_values,
            method=MULTIPLE_TESTING_CORRECTION_BONFERRONI,
        )
        np.testing.assert_allclose(table["adj.P.Val"].to_numpy(dtype=float), expected)
        np.testing.assert_allclose(
            table["adj.P.Val"].to_numpy(dtype=float),
            np.clip(p_values * 2.0, 0.0, 1.0),
        )
        assert not np.array_equal(
            table["adj.P.Val"].to_numpy(dtype=float),
            np.clip(p_values * 3.0, 0.0, 1.0),
        )


def test_post_fit_non_finite_site_does_not_abort_valid_shared_group() -> None:
    protein = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    overflowing = np.array(
        [1.0e308, -1.0e308, 1.0e308, -1.0e308, 1.0e308, -1.0e308],
        dtype=float,
    )
    request = _request(
        matrix=_matrix(
            {
                "valid_site": _site_values(5.0, 6.0, 0.2, protein),
                "overflowing_site": overflowing,
            }
        ),
        proteins=_proteins({"protein_a": protein}),
        pairs=_pairs(
            {
                "valid_site": ("MAPK14", "protein_a"),
                "overflowing_site": ("AKT1", "protein_a"),
            }
        ),
    )

    with np.errstate(over="ignore"):
        result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert result.tested_site_ids == ("valid_site",)
    assert result.table_for("B_vs_A").index.tolist() == ["valid_site"]
    assert result.prior_diagnostics.prior_variance.index.tolist() == ["valid_site"]
    assert result.augmented_design_diagnostics_dataframe().index.tolist() == [
        "protein_a"
    ]
    failures = result.site_failure_diagnostics_dataframe()
    assert failures.index.tolist() == ["overflowing_site"]
    assert failures.loc["overflowing_site", DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID
    )
    assert (
        failures.loc[
            "overflowing_site",
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        ]
        == DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_FINITE
    )


def test_zero_residual_variance_site_is_withheld_before_moderation() -> None:
    protein = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    request = _request(
        matrix=_matrix(
            {
                "exact_site": _site_values(
                    5.0,
                    6.0,
                    0.2,
                    protein,
                    noise_scale=0.0,
                ),
                "noisy_site_a": _site_values(4.0, 4.8, -0.3, protein),
                "noisy_site_b": _site_values(
                    2.0,
                    2.5,
                    0.4,
                    protein,
                    seed=np.roll(NOISE_SEED, 1),
                ),
            }
        ),
        proteins=_proteins({"protein_a": protein}),
        pairs=_pairs(
            {
                "exact_site": ("MAPK14", "protein_a"),
                "noisy_site_a": ("AKT1", "protein_a"),
                "noisy_site_b": ("GSK3B", "protein_a"),
            }
        ),
        multiple_testing_method=MULTIPLE_TESTING_CORRECTION_BONFERRONI,
    )

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)

    assert result.tested_site_ids == ("noisy_site_a", "noisy_site_b")
    residual_variance = result.residual_variance_series()
    assert np.isfinite(residual_variance.to_numpy(dtype=float)).all()
    assert (residual_variance.to_numpy(dtype=float) > 0.0).all()
    failures = result.site_failure_diagnostics_dataframe()
    assert failures.loc["exact_site", DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID
    )
    assert failures.loc["exact_site", DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_POSITIVE
    )
    table = result.table_for("B_vs_A")
    p_values = table["P.Value"].to_numpy(dtype=float)
    np.testing.assert_allclose(
        table["adj.P.Val"].to_numpy(dtype=float),
        adjust_p_values(p_values, method=MULTIPLE_TESTING_CORRECTION_BONFERRONI),
    )
    np.testing.assert_allclose(
        table["adj.P.Val"].to_numpy(dtype=float),
        np.clip(p_values * 2.0, 0.0, 1.0),
    )
    assert not np.array_equal(
        table["adj.P.Val"].to_numpy(dtype=float),
        np.clip(p_values * 3.0, 0.0, 1.0),
    )
    exported_failures = result.site_failure_diagnostics_dataframe()
    exported_failures.loc["exact_site", DIFFERENTIAL_RESULT_STATUS_COLUMN] = "changed"
    assert (
        result.site_failure_diagnostics_dataframe().loc[
            "exact_site",
            DIFFERENTIAL_RESULT_STATUS_COLUMN,
        ]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID
    )


def test_all_post_fit_failed_sites_raise_typed_global_error() -> None:
    protein = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    request = _request(
        matrix=_matrix(
            {
                "exact_site_a": _site_values(
                    5.0,
                    6.0,
                    0.2,
                    protein,
                    noise_scale=0.0,
                ),
                "exact_site_b": _site_values(
                    4.0,
                    4.5,
                    -0.1,
                    protein,
                    noise_scale=0.0,
                ),
            }
        ),
        proteins=_proteins({"protein_a": protein}),
        pairs=_pairs(
            {
                "exact_site_a": ("MAPK14", "protein_a"),
                "exact_site_b": ("AKT1", "protein_a"),
            }
        ),
    )

    with pytest.raises(PhosPyInputError, match="post-fit numerical eligibility") as exc:
        ProteinCovariateAdjustedDifferentialKernel().run(request)

    report = exc.value.diagnostics
    assert isinstance(report, dict)
    assert report["status_counts"] == {
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID: 2
    }
    assert report["reason_counts"] == {
        DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_POSITIVE: 2
    }
    failures = report["site_failure_diagnostics"]
    assert isinstance(failures, pd.DataFrame)
    assert failures.index.tolist() == ["exact_site_a", "exact_site_b"]
    assert "Traceback" not in str(exc.value)


def test_kernel_does_not_mutate_request_owned_frames() -> None:
    protein = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0], dtype=float)
    request = _request(
        matrix=_matrix({"site_a": _site_values(5.0, 6.0, 0.25, protein)}),
        proteins=_proteins({"protein_a": protein}),
        pairs=_pairs({"site_a": ("MAPK14", "protein_a")}),
    )
    matrix_before = request.phosphosite_matrix.copy(deep=True)
    design_before = request.base_design.to_dataframe()
    contrasts_before = request.base_contrasts.to_dataframe()
    pairs_before = request.matched_pairs.copy(deep=True)
    proteins_before = request.resolved_protein_covariates.copy(deep=True)

    result = ProteinCovariateAdjustedDifferentialKernel().run(request)
    exported_coefficients = result.coefficient_dataframe()
    exported_coefficients.iloc[0, 0] = 999.0

    pd.testing.assert_frame_equal(request.phosphosite_matrix, matrix_before)
    pd.testing.assert_frame_equal(request.base_design.to_dataframe(), design_before)
    pd.testing.assert_frame_equal(
        request.base_contrasts.to_dataframe(), contrasts_before
    )
    pd.testing.assert_frame_equal(request.matched_pairs, pairs_before)
    pd.testing.assert_frame_equal(request.resolved_protein_covariates, proteins_before)
    assert result.coefficient_dataframe().iloc[0, 0] != exported_coefficients.iloc[0, 0]


def _site_values(
    beta_a: float,
    beta_b: float,
    gamma: float,
    protein: np.ndarray,
    *,
    seed: np.ndarray = NOISE_SEED,
    noise_scale: float = 0.05,
) -> np.ndarray:
    centered = protein - float(np.mean(protein))
    design = np.column_stack([CONDITION_A, CONDITION_B, centered])
    noise = _orthogonal_noise(design, seed)
    return (
        beta_a * CONDITION_A
        + beta_b * CONDITION_B
        + gamma * centered
        + noise_scale * noise
    )


def _orthogonal_noise(design: np.ndarray, seed: np.ndarray) -> np.ndarray:
    projected = design @ np.linalg.lstsq(design, seed, rcond=None)[0]
    residual = seed - projected
    norm = float(np.linalg.norm(residual))
    if norm <= 0.0:
        raise AssertionError("test fixture noise seed is not independent of design")
    return residual / norm


def _request(
    *,
    matrix: pd.DataFrame,
    proteins: pd.DataFrame,
    pairs: pd.DataFrame,
    design: pd.DataFrame | None = None,
    contrasts: pd.DataFrame | None = None,
    sample_order: tuple[str, ...] = SAMPLES,
    empirical_bayes: EmpiricalBayesConfig | None = None,
    multiple_testing_method: str = "benjamini_hochberg",
) -> ProteinAwareDifferentialComputationRequest:
    return ProteinAwareDifferentialComputationRequest(
        phosphosite_matrix=matrix,
        base_design=_base_design(sample_order) if design is None else design,
        base_contrasts=_contrasts() if contrasts is None else contrasts,
        sample_order=sample_order,
        matched_pairs=pairs,
        resolved_protein_covariates=proteins,
        empirical_bayes=(
            EmpiricalBayesConfig() if empirical_bayes is None else empirical_bayes
        ),
        multiple_testing_method=multiple_testing_method,
    )


def _base_design(samples: tuple[str, ...] = SAMPLES) -> pd.DataFrame:
    if samples != SAMPLES:
        raise AssertionError("custom sample sets must supply an explicit design")
    return pd.DataFrame(
        {"A": CONDITION_A, "B": CONDITION_B},
        index=pd.Index(samples, name="sample_id"),
    )


def _contrasts(*, multiple: bool = False) -> pd.DataFrame:
    data = {"B_vs_A": [-1.0, 1.0]}
    if multiple:
        data["A_level"] = [1.0, 0.0]
    return pd.DataFrame(
        data,
        index=pd.Index(["A", "B"], name="coefficient"),
    )


def _matrix(values_by_site: dict[str, np.ndarray]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            sample: [float(values[position]) for values in values_by_site.values()]
            for position, sample in enumerate(SAMPLES)
        },
        index=pd.Index(tuple(values_by_site), name="site_key"),
    )


def _proteins(
    values_by_protein: dict[str, np.ndarray],
    *,
    samples: tuple[str, ...] = SAMPLES,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            sample: [float(values[position]) for values in values_by_protein.values()]
            for position, sample in enumerate(samples)
        },
        index=pd.Index(tuple(values_by_protein), name="total_protein_row_key"),
    )


def _pairs(mapping: dict[str, tuple[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": list(mapping),
            "protein_identifier": [value[0] for value in mapping.values()],
            "total_protein_row_key": [value[1] for value in mapping.values()],
        }
    )


def _oracle(
    *,
    values: np.ndarray,
    protein: np.ndarray,
    design: pd.DataFrame | None = None,
    contrasts: pd.DataFrame | None = None,
) -> dict[str, np.ndarray | float]:
    base_design = _base_design() if design is None else design
    base_contrasts = _contrasts() if contrasts is None else contrasts
    centered = protein - float(np.mean(protein))
    augmented_design = base_design.copy(deep=True)
    augmented_design.loc[:, PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME] = centered
    augmented_contrasts = pd.concat(
        [
            base_contrasts.copy(deep=True),
            pd.DataFrame(
                np.zeros((1, int(base_contrasts.shape[1])), dtype=float),
                index=pd.Index(
                    [PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME],
                    name=base_contrasts.index.name,
                ),
                columns=base_contrasts.columns.copy(),
            ),
        ],
        axis=0,
    )
    design_values = augmented_design.to_numpy(dtype=float)
    contrast_values = augmented_contrasts.to_numpy(dtype=float)
    coefficients = np.linalg.lstsq(design_values, values, rcond=None)[0]
    fitted = design_values @ coefficients
    residuals = values - fitted
    rank = np.linalg.matrix_rank(design_values)
    residual_dof = float(design_values.shape[0] - rank)
    covariance = np.linalg.inv(design_values.T @ design_values)
    contrast_covariance = contrast_values.T @ covariance @ contrast_values
    return {
        "coefficients": coefficients,
        "residuals": residuals,
        "residual_variance": float(residuals @ residuals / residual_dof),
        "contrast_effects": coefficients @ contrast_values,
        "contrast_scales": np.sqrt(np.diag(contrast_covariance)),
    }
