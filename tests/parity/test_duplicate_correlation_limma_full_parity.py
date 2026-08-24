from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy.science.differential.executor import (
    DuplicateCorrelationDifferentialAnalysisExecutor,
    DuplicateCorrelationDifferentialFit,
)
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisRequest,
    EmpiricalBayesConfig,
)
from tests.support.duplicate_correlation_parity_scopes import (
    FULL_PIPELINE_FIXTURE_IDS,
)

pytestmark = [pytest.mark.parity]

FIXTURE_ROOT = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_duplicate_correlation"
)


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationParityTolerances:
    linear_fit_atol: float = 1.0e-8
    variance_atol: float = 5.0e-9
    degrees_of_freedom_atol: float = 1.0e-12
    standard_error_atol: float = 5.0e-12
    moderated_statistic_atol: float = 1.0e-8
    p_value_atol: float = 5.0e-12
    relative: float = 1.0e-7


TOLERANCES = DuplicateCorrelationParityTolerances()
SERIALIZED_POSITIVE_INFINITY_COLUMNS = frozenset({"df.prior"})


@pytest.mark.parametrize("fixture_id", FULL_PIPELINE_FIXTURE_IDS)
def test_duplicate_correlation_ebayes_priors_match_limma_fixture(
    fixture_id: str,
) -> None:
    fixture = _load_fixture(fixture_id)
    fit = _run_duplicate_correlation_fit(fixture)

    observed = _observed_ebayes_priors(fit, feature_ids=fixture.matrix.index)
    expected = pd.read_csv(
        FIXTURE_ROOT / fixture_id / "ebayes_priors.csv",
        index_col="feature_id",
    )

    _assert_frame_columns_close(
        observed,
        expected,
        absolute_tolerances={
            "sigma": TOLERANCES.variance_atol,
            "residual_variance": TOLERANCES.variance_atol,
            "residual_degrees_of_freedom": TOLERANCES.degrees_of_freedom_atol,
            "s2.prior": TOLERANCES.variance_atol,
            "df.prior": TOLERANCES.degrees_of_freedom_atol,
            "s2.post": TOLERANCES.variance_atol,
            "df.total": TOLERANCES.degrees_of_freedom_atol,
            "average_expression": TOLERANCES.linear_fit_atol,
        },
    )


@pytest.mark.parametrize("fixture_id", FULL_PIPELINE_FIXTURE_IDS)
def test_duplicate_correlation_public_and_derived_statistics_match_limma_fixture(
    fixture_id: str,
) -> None:
    fixture = _load_fixture(fixture_id)
    fit = _run_duplicate_correlation_fit(fixture)
    expected_all = pd.read_csv(FIXTURE_ROOT / fixture_id / "ebayes_statistics.csv")

    assert tuple(fit.gls_fit.contrast_names) == tuple(fixture.contrasts.columns)
    for contrast_name in fixture.contrasts.columns.astype(str).tolist():
        observed = _observed_contrast_statistics(
            fit,
            contrast_name=contrast_name,
            feature_ids=fixture.matrix.index,
        )
        expected = (
            expected_all.loc[expected_all["contrast"].astype(str) == contrast_name]
            .set_index("feature_id")
            .loc[:, observed.columns]
        )
        _assert_frame_columns_close(
            observed,
            expected,
            absolute_tolerances={
                "logFC": TOLERANCES.linear_fit_atol,
                "AveExpr": TOLERANCES.linear_fit_atol,
                "SE": TOLERANCES.standard_error_atol,
                "t": TOLERANCES.moderated_statistic_atol,
                "P.Value": TOLERANCES.p_value_atol,
                "adj.P.Val": TOLERANCES.p_value_atol,
                "stdev.unscaled": TOLERANCES.standard_error_atol,
                "sigma": TOLERANCES.variance_atol,
                "residual_degrees_of_freedom": (TOLERANCES.degrees_of_freedom_atol),
                "s2.prior": TOLERANCES.variance_atol,
                "df.prior": TOLERANCES.degrees_of_freedom_atol,
                "s2.post": TOLERANCES.variance_atol,
                "df.total": TOLERANCES.degrees_of_freedom_atol,
            },
        )


def test_duplicate_correlation_full_fit_is_invariant_to_feature_order() -> None:
    fixture = _load_fixture("fixture_b_three_observation_blocks")
    baseline = _run_duplicate_correlation_fit(fixture)
    permutation = [5, 0, 9, 1, 13, 3, 7, 11, 2, 4, 6, 8, 10, 12]
    permuted_fixture = _FixtureInputs(
        matrix=fixture.matrix.iloc[permutation, :],
        design=fixture.design,
        contrasts=fixture.contrasts,
        blocks=fixture.blocks,
    )

    permuted = _run_duplicate_correlation_fit(permuted_fixture)

    assert permuted.consensus.consensus_correlation == pytest.approx(
        baseline.consensus.consensus_correlation,
        abs=TOLERANCES.p_value_atol,
    )
    _assert_computation_result_matches_after_feature_reindex(
        observed=permuted,
        expected=baseline,
    )


def test_duplicate_correlation_full_fit_is_invariant_to_sample_order_and_block_labels() -> (
    None
):
    fixture = _load_fixture("fixture_b_three_observation_blocks")
    baseline = _run_duplicate_correlation_fit(fixture)
    permutation = [
        5,
        1,
        9,
        0,
        7,
        13,
        3,
        11,
        2,
        4,
        6,
        8,
        10,
        12,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
    ]
    renamed = {
        block: f"renamed_{position}"
        for position, block in enumerate(sorted(set(fixture.blocks)), start=1)
    }
    permuted_fixture = _FixtureInputs(
        matrix=fixture.matrix.iloc[:, permutation],
        design=fixture.design.iloc[permutation, :],
        contrasts=fixture.contrasts,
        blocks=tuple(renamed[fixture.blocks[position]] for position in permutation),
    )

    permuted = _run_duplicate_correlation_fit(permuted_fixture)

    assert permuted.consensus.consensus_correlation == pytest.approx(
        baseline.consensus.consensus_correlation,
        abs=TOLERANCES.p_value_atol,
    )
    _assert_computation_result_matches_after_feature_reindex(
        observed=permuted,
        expected=baseline,
    )


def test_duplicate_correlation_full_fit_accepts_equivalent_contrast_parameterisation() -> (
    None
):
    fixture = _load_fixture("fixture_c_incomplete_unequal_blocks")
    baseline = _run_duplicate_correlation_fit(fixture)
    alternative_design = pd.DataFrame(
        {
            "Intercept": np.ones(fixture.design.shape[0], dtype=float),
            "B_effect": fixture.design["B"].to_numpy(dtype=float),
            "C_effect": fixture.design["C"].to_numpy(dtype=float),
        },
        index=fixture.design.index.copy(),
    )
    alternative_contrasts = pd.DataFrame(
        {
            "B_vs_A": [0.0, 1.0, 0.0],
            "C_vs_A": [0.0, 0.0, 1.0],
        },
        index=pd.Index(("Intercept", "B_effect", "C_effect"), name="coefficient"),
    )
    alternative_fixture = _FixtureInputs(
        matrix=fixture.matrix,
        design=alternative_design,
        contrasts=alternative_contrasts,
        blocks=fixture.blocks,
    )

    alternative = _run_duplicate_correlation_fit(alternative_fixture)

    assert alternative.consensus.consensus_correlation == pytest.approx(
        baseline.consensus.consensus_correlation,
        abs=TOLERANCES.p_value_atol,
    )
    _assert_computation_result_matches_after_feature_reindex(
        observed=alternative,
        expected=baseline,
    )


@dataclass(frozen=True, slots=True)
class _FixtureInputs:
    matrix: pd.DataFrame
    design: pd.DataFrame
    contrasts: pd.DataFrame
    blocks: tuple[str, ...]


def _load_fixture(fixture_id: str) -> _FixtureInputs:
    base = FIXTURE_ROOT / fixture_id
    return _FixtureInputs(
        matrix=pd.read_csv(base / "matrix.csv", index_col="feature_id"),
        design=pd.read_csv(base / "design.csv", index_col="sample_id"),
        contrasts=pd.read_csv(base / "contrasts.csv", index_col="coefficient"),
        blocks=tuple(pd.read_csv(base / "blocks.csv")["block_id"].astype(str).tolist()),
    )


def _run_duplicate_correlation_fit(
    fixture: _FixtureInputs,
) -> DuplicateCorrelationDifferentialFit:
    request = DifferentialAnalysisRequest(
        matrix=fixture.matrix,
        design=DesignMatrix(fixture.design),
        contrasts=ContrastMatrix(fixture.contrasts),
        empirical_bayes=EmpiricalBayesConfig(method="standard"),
    )
    return DuplicateCorrelationDifferentialAnalysisExecutor().run(
        request,
        block_ids=fixture.blocks,
    )


def _observed_ebayes_priors(
    fit: DuplicateCorrelationDifferentialFit,
    *,
    feature_ids: pd.Index,
) -> pd.DataFrame:
    result = fit.computation_result
    residual_dof = float(result.residual_degrees_of_freedom)
    feature_count = int(feature_ids.size)
    prior_dof = result.prior_degrees_of_freedom_series_value.to_numpy(dtype=float)
    return pd.DataFrame(
        {
            "sigma": np.sqrt(result.residual_variance.to_numpy(dtype=float)),
            "residual_variance": result.residual_variance.to_numpy(dtype=float),
            "residual_degrees_of_freedom": np.repeat(residual_dof, feature_count),
            "s2.prior": result.prior_residual_variance.to_numpy(dtype=float),
            "df.prior": prior_dof,
            "s2.post": result.posterior_residual_variance.to_numpy(dtype=float),
            "df.total": _moderated_degrees_of_freedom(
                prior_dof,
                residual_dof=residual_dof,
                feature_count=feature_count,
            ),
            "average_expression": fit.gls_fit.average_expression,
        },
        index=pd.Index(feature_ids.astype(str), name="feature_id"),
    )


def _observed_contrast_statistics(
    fit: DuplicateCorrelationDifferentialFit,
    *,
    contrast_name: str,
    feature_ids: pd.Index,
) -> pd.DataFrame:
    result = fit.computation_result
    contrast_position = fit.gls_fit.contrast_names.index(contrast_name)
    contrast_table = result.table_for(contrast_name)
    residual_dof = float(result.residual_degrees_of_freedom)
    feature_count = int(feature_ids.size)
    prior_dof = result.prior_degrees_of_freedom_series_value.to_numpy(dtype=float)
    posterior_variance = result.posterior_residual_variance.to_numpy(dtype=float)
    contrast_stdev_unscaled = fit.gls_fit.contrast_stdev_unscaled
    assert contrast_stdev_unscaled is not None
    stdev_unscaled = contrast_stdev_unscaled[contrast_position, :]
    return pd.DataFrame(
        {
            "logFC": contrast_table["logFC"].to_numpy(dtype=float),
            "AveExpr": fit.gls_fit.average_expression,
            "SE": np.sqrt(posterior_variance) * stdev_unscaled,
            "t": contrast_table["t"].to_numpy(dtype=float),
            "P.Value": contrast_table["P.Value"].to_numpy(dtype=float),
            "adj.P.Val": contrast_table["adj.P.Val"].to_numpy(dtype=float),
            "stdev.unscaled": stdev_unscaled,
            "sigma": np.sqrt(result.residual_variance.to_numpy(dtype=float)),
            "residual_degrees_of_freedom": np.repeat(residual_dof, feature_count),
            "s2.prior": result.prior_residual_variance.to_numpy(dtype=float),
            "df.prior": prior_dof,
            "s2.post": posterior_variance,
            "df.total": _moderated_degrees_of_freedom(
                prior_dof,
                residual_dof=residual_dof,
                feature_count=feature_count,
            ),
        },
        index=pd.Index(feature_ids.astype(str), name="feature_id"),
    )


def _moderated_degrees_of_freedom(
    prior_dof: np.ndarray,
    *,
    residual_dof: float,
    feature_count: int,
) -> np.ndarray:
    total_residual_dof = residual_dof * float(feature_count)
    moderated = np.full(prior_dof.shape, total_residual_dof, dtype=float)
    finite_prior = np.isfinite(prior_dof)
    moderated[finite_prior] = np.minimum(
        residual_dof + prior_dof[finite_prior],
        total_residual_dof,
    )
    return moderated


def _assert_frame_columns_close(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    absolute_tolerances: dict[str, float],
) -> None:
    assert observed.index.tolist() == expected.index.astype(str).tolist()
    assert (
        observed.columns.astype(str).tolist() == expected.columns.astype(str).tolist()
    )
    for column_name, absolute_tolerance in absolute_tolerances.items():
        observed_values = observed[column_name].to_numpy(dtype=np.float64)
        expected_values = expected[column_name].to_numpy(dtype=np.float64)
        if column_name in SERIALIZED_POSITIVE_INFINITY_COLUMNS:
            serialized_infinite = np.isnan(expected_values) & np.isposinf(
                observed_values
            )
            observed_values = observed_values[~serialized_infinite]
            expected_values = expected_values[~serialized_infinite]
        np.testing.assert_allclose(
            observed_values,
            expected_values,
            rtol=TOLERANCES.relative,
            atol=absolute_tolerance,
            equal_nan=True,
        )


def _assert_computation_result_matches_after_feature_reindex(
    *,
    observed: DuplicateCorrelationDifferentialFit,
    expected: DuplicateCorrelationDifferentialFit,
) -> None:
    observed_result = observed.computation_result
    expected_result = expected.computation_result
    expected_index = expected_result.residual_variance.index
    observed_index = observed_result.residual_variance.index
    assert sorted(observed_index.astype(str)) == sorted(expected_index.astype(str))

    for observed_series, expected_series, absolute_tolerance in (
        (
            observed_result.residual_variance,
            expected_result.residual_variance,
            TOLERANCES.variance_atol,
        ),
        (
            observed_result.posterior_residual_variance,
            expected_result.posterior_residual_variance,
            TOLERANCES.variance_atol,
        ),
        (
            observed_result.prior_residual_variance,
            expected_result.prior_residual_variance,
            TOLERANCES.variance_atol,
        ),
        (
            observed_result.prior_degrees_of_freedom_series_value,
            expected_result.prior_degrees_of_freedom_series_value,
            TOLERANCES.degrees_of_freedom_atol,
        ),
    ):
        np.testing.assert_allclose(
            observed_series.reindex(expected_index).to_numpy(dtype=np.float64),
            expected_series.to_numpy(dtype=np.float64),
            rtol=TOLERANCES.relative,
            atol=absolute_tolerance,
            equal_nan=True,
        )

    for contrast_name in expected.gls_fit.contrast_names:
        observed_table = observed_result.table_for(contrast_name).reindex(
            expected_index
        )
        expected_table = expected_result.table_for(contrast_name)
        _assert_frame_columns_close(
            observed_table.loc[:, ("logFC", "t", "P.Value", "adj.P.Val")],
            expected_table.loc[:, ("logFC", "t", "P.Value", "adj.P.Val")],
            absolute_tolerances={
                "logFC": TOLERANCES.linear_fit_atol,
                "t": TOLERANCES.moderated_statistic_atol,
                "P.Value": TOLERANCES.p_value_atol,
                "adj.P.Val": TOLERANCES.p_value_atol,
            },
        )
