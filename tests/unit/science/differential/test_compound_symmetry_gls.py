from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest

import phospy.science.differential.compound_symmetry_gls as gls_module
from phospy.errors import PhosPyInputError
from phospy.science.differential.compound_symmetry_gls import (
    COMPOUND_SYMMETRY_GLS_STATUS_PARTIAL_RANK_LOSS,
    fit_compound_symmetry_gls,
)
from phospy.science.differential.linear_model import decompose_differential_design
from tests.support.duplicate_correlation_parity_scopes import (
    ESTIMATOR_AND_GLS_ONLY_FIXTURE_IDS,
    FULL_PIPELINE_FIXTURE_IDS,
)

ROOT = Path(__file__).resolve().parents[4]
FIXTURE_ROOT = (
    ROOT
    / "tests"
    / "fixtures"
    / "rewrite_parity"
    / "differential_duplicate_correlation"
)

FIXTURE_IDS = (*FULL_PIPELINE_FIXTURE_IDS, *ESTIMATOR_AND_GLS_ONLY_FIXTURE_IDS)
GLS_LINEAR_ABSOLUTE_TOLERANCE = 5.0e-9
GLS_STDEV_UNSCALED_ABSOLUTE_TOLERANCE = 5.0e-9
GLS_COVARIANCE_ABSOLUTE_TOLERANCE = 5.0e-9
GLS_VARIANCE_ABSOLUTE_TOLERANCE = 5.0e-9
GLS_DOF_ABSOLUTE_TOLERANCE = 1.0e-12
GLS_PARITY_RELATIVE_TOLERANCE = 5.0e-9
STRICT_ABSOLUTE_TOLERANCE = 1.0e-12


@pytest.mark.parametrize("fixture_id", FIXTURE_IDS)
def test_compound_symmetry_gls_matches_limma_upstream_fit_fields(
    fixture_id: str,
) -> None:
    matrix, design, blocks, contrasts = _fixture_inputs(fixture_id)
    summary = _summary(fixture_id)

    fit = fit_compound_symmetry_gls(
        matrix.to_numpy(dtype=np.float64),
        design.to_numpy(dtype=np.float64),
        blocks,
        float(summary["consensus_correlation"]),
        feature_ids=tuple(str(value) for value in matrix.index.tolist()),
        coefficient_names=tuple(str(value) for value in design.columns.tolist()),
        contrasts=contrasts.to_numpy(dtype=np.float64),
        contrast_names=tuple(str(value) for value in contrasts.columns.tolist()),
    )

    _assert_frame_close(
        _coefficient_frame(fit.coefficients, fit),
        _numeric_csv(fixture_id, "fit_coefficients.csv", index_col="feature_id"),
        absolute_tolerance=GLS_LINEAR_ABSOLUTE_TOLERANCE,
    )
    _assert_frame_close(
        _coefficient_frame(fit.stdev_unscaled, fit),
        _numeric_csv(fixture_id, "fit_stdev_unscaled.csv", index_col="feature_id"),
        absolute_tolerance=GLS_STDEV_UNSCALED_ABSOLUTE_TOLERANCE,
    )
    _assert_frame_close(
        _covariance_frame(fit.coefficient_covariance, fit.coefficient_names),
        _numeric_csv(fixture_id, "fit_cov_coefficients.csv", index_col="coefficient"),
        absolute_tolerance=GLS_COVARIANCE_ABSOLUTE_TOLERANCE,
    )
    _assert_frame_columns_close(
        pd.DataFrame(
            {
                "sigma": fit.residual_standard_deviation,
                "residual_variance": fit.residual_variance,
                "residual_degrees_of_freedom": fit.residual_degrees_of_freedom,
                "average_expression": fit.average_expression,
            },
            index=pd.Index(fit.feature_ids, name="feature_id"),
        ),
        _numeric_csv(fixture_id, "fit_sigma_df.csv", index_col="feature_id"),
        absolute_tolerances={
            "sigma": GLS_VARIANCE_ABSOLUTE_TOLERANCE,
            "residual_variance": GLS_VARIANCE_ABSOLUTE_TOLERANCE,
            "residual_degrees_of_freedom": GLS_DOF_ABSOLUTE_TOLERANCE,
            "average_expression": GLS_LINEAR_ABSOLUTE_TOLERANCE,
        },
    )
    assert fit.contrast_coefficients is not None
    assert fit.contrast_stdev_unscaled is not None
    assert fit.contrast_covariance is not None
    _assert_frame_close(
        _contrast_feature_frame(fit.contrast_coefficients, fit),
        _numeric_csv(
            fixture_id,
            "contrast_fit_coefficients.csv",
            index_col="feature_id",
        ),
        absolute_tolerance=GLS_LINEAR_ABSOLUTE_TOLERANCE,
    )
    _assert_frame_close(
        _contrast_feature_frame(fit.contrast_stdev_unscaled, fit),
        _numeric_csv(
            fixture_id,
            "contrast_fit_stdev_unscaled.csv",
            index_col="feature_id",
        ),
        absolute_tolerance=GLS_STDEV_UNSCALED_ABSOLUTE_TOLERANCE,
    )
    _assert_frame_close(
        _covariance_frame(fit.contrast_covariance, fit.contrast_names),
        _numeric_csv(
            fixture_id,
            "contrast_fit_cov_coefficients.csv",
            index_col="contrast",
        ),
        absolute_tolerance=GLS_COVARIANCE_ABSOLUTE_TOLERANCE,
    )


def test_zero_correlation_gls_reduces_to_existing_ols_decomposition() -> None:
    design = np.array(
        [
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [0.0, 1.0, 2.0],
        ],
        dtype=np.float64,
    )
    matrix = np.array(
        [
            [1.0, 1.2, 2.0, 2.2, 2.4],
            [4.0, 4.1, 3.0, 3.3, 3.5],
            [8.0, 8.4, 9.0, 9.2, 9.5],
        ],
        dtype=np.float64,
    )
    blocks = ("b1", "b1", "b2", "b2", "b3")
    contrasts = np.array([[-1.0], [1.0], [0.0]], dtype=np.float64)

    gls = fit_compound_symmetry_gls(
        matrix,
        design,
        blocks,
        0.0,
        contrasts=contrasts,
    )
    ols_decomposition = decompose_differential_design(design)
    ols = ols_decomposition.fit(matrix.T)

    np.testing.assert_allclose(
        gls.coefficients,
        ols.coefficients,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        gls.residual_variance,
        ols.residual_variance,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        gls.residual_sum_of_squares,
        ols.residual_sum_of_squares,
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        gls.stdev_unscaled,
        np.repeat(
            np.sqrt(np.diagonal(ols_decomposition.coefficient_covariance))[
                :,
                np.newaxis,
            ],
            matrix.shape[0],
            axis=1,
        ),
        rtol=1.0e-12,
        atol=1.0e-12,
    )
    assert np.all(gls.residual_degrees_of_freedom == pytest.approx(2.0))
    np.testing.assert_allclose(
        gls.coefficient_covariance,
        ols_decomposition.coefficient_covariance,
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_positive_and_valid_negative_correlations_fit_without_policy_fallback() -> None:
    matrix = np.array(
        [
            [1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
            [3.0, 2.7, 2.3, 2.2, 2.0, 1.8],
        ],
        dtype=np.float64,
    )
    design = np.array(
        [
            [1.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
        ],
        dtype=np.float64,
    )
    blocks = ("p1", "p1", "p2", "p2", "p3", "p3")

    positive = fit_compound_symmetry_gls(matrix, design, blocks, 0.35)
    negative = fit_compound_symmetry_gls(matrix, design, blocks, -0.35)

    assert positive.consensus_correlation == pytest.approx(0.35)
    assert negative.consensus_correlation == pytest.approx(-0.35)
    assert np.isfinite(positive.residual_standard_deviation).all()
    assert np.isfinite(negative.residual_standard_deviation).all()
    assert not np.allclose(positive.stdev_unscaled, negative.stdev_unscaled)


def test_invalid_correlation_is_rejected_before_fitting() -> None:
    matrix = np.array([[1.0, 2.0, 3.0]], dtype=np.float64)
    design = np.ones((3, 1), dtype=np.float64)
    blocks = ("b1", "b1", "b1")

    with pytest.raises(PhosPyInputError, match="positive-definite"):
        fit_compound_symmetry_gls(matrix, design, blocks, -0.5)


def test_cholesky_failure_is_explicit_and_does_not_fall_back_to_ols(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _failing_cholesky(_: npt.ArrayLike) -> npt.NDArray[np.float64]:
        raise np.linalg.LinAlgError("forced test failure")

    monkeypatch.setattr(gls_module.np.linalg, "cholesky", _failing_cholesky)

    with pytest.raises(PhosPyInputError, match="positive-definite|Cholesky"):
        fit_compound_symmetry_gls(
            np.array([[1.0, 2.0, 3.0, 4.0]], dtype=np.float64),
            np.array([[1.0], [1.0], [1.0], [1.0]], dtype=np.float64),
            ("b1", "b1", "b2", "b2"),
            0.1,
        )


def test_sample_order_and_block_label_invariance() -> None:
    matrix, design, blocks, contrasts = _fixture_inputs(
        "fixture_b_three_observation_blocks"
    )
    summary = _summary("fixture_b_three_observation_blocks")
    permutation = np.array(
        [
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
        ],
        dtype=np.int64,
    )
    renamed_by_block = {
        block: f"renamed_{position}"
        for position, block in enumerate(sorted(set(blocks)), start=1)
    }

    baseline = fit_compound_symmetry_gls(
        matrix.to_numpy(dtype=np.float64),
        design.to_numpy(dtype=np.float64),
        blocks,
        float(summary["consensus_correlation"]),
        contrasts=contrasts.to_numpy(dtype=np.float64),
    )
    permuted = fit_compound_symmetry_gls(
        matrix.iloc[:, permutation].to_numpy(dtype=np.float64),
        design.iloc[permutation, :].to_numpy(dtype=np.float64),
        tuple(blocks[int(position)] for position in permutation.tolist()),
        float(summary["consensus_correlation"]),
        contrasts=contrasts.to_numpy(dtype=np.float64),
    )
    renamed = fit_compound_symmetry_gls(
        matrix.to_numpy(dtype=np.float64),
        design.to_numpy(dtype=np.float64),
        tuple(renamed_by_block[block] for block in blocks),
        float(summary["consensus_correlation"]),
        contrasts=contrasts.to_numpy(dtype=np.float64),
    )

    _assert_gls_numeric_fit_close(permuted, baseline)
    _assert_gls_numeric_fit_close(renamed, baseline)


def test_explicit_missingness_mask_omits_values_without_zero_filling() -> None:
    matrix, design, blocks, contrasts = _fixture_inputs(
        "fixture_d_feature_level_failures"
    )
    summary = _summary("fixture_d_feature_level_failures")
    nan_fit = fit_compound_symmetry_gls(
        matrix.to_numpy(dtype=np.float64),
        design.to_numpy(dtype=np.float64),
        blocks,
        float(summary["consensus_correlation"]),
        feature_ids=tuple(str(value) for value in matrix.index.tolist()),
        coefficient_names=tuple(str(value) for value in design.columns.tolist()),
        contrasts=contrasts.to_numpy(dtype=np.float64),
        contrast_names=tuple(str(value) for value in contrasts.columns.tolist()),
    )
    filled_matrix = matrix.fillna(-9999.0)
    explicit_mask = matrix.notna().to_numpy(dtype=bool)

    masked_fit = fit_compound_symmetry_gls(
        filled_matrix.to_numpy(dtype=np.float64),
        design.to_numpy(dtype=np.float64),
        blocks,
        float(summary["consensus_correlation"]),
        feature_ids=tuple(str(value) for value in matrix.index.tolist()),
        coefficient_names=tuple(str(value) for value in design.columns.tolist()),
        observation_mask=explicit_mask,
        contrasts=contrasts.to_numpy(dtype=np.float64),
        contrast_names=tuple(str(value) for value in contrasts.columns.tolist()),
    )

    assert masked_fit.factorization_cache_size == nan_fit.factorization_cache_size
    _assert_gls_numeric_fit_close(masked_fit, nan_fit)


def test_repeated_missingness_masks_reuse_gls_factorizations() -> None:
    matrix, design, blocks, contrasts = _fixture_inputs(
        "fixture_d_feature_level_failures"
    )
    summary = _summary("fixture_d_feature_level_failures")

    fit = fit_compound_symmetry_gls(
        matrix.to_numpy(dtype=np.float64),
        design.to_numpy(dtype=np.float64),
        blocks,
        float(summary["consensus_correlation"]),
        feature_ids=tuple(str(value) for value in matrix.index.tolist()),
        coefficient_names=tuple(str(value) for value in design.columns.tolist()),
        contrasts=contrasts.to_numpy(dtype=np.float64),
        contrast_names=tuple(str(value) for value in contrasts.columns.tolist()),
    )

    row_masks = [tuple(row) for row in np.isfinite(matrix.to_numpy(dtype=float))]
    unique_masks = set(row_masks)
    full_mask = tuple(True for _ in range(matrix.shape[1]))
    seen_masks = {full_mask}
    expected_cache_hits = 0
    for row_mask in row_masks:
        if row_mask in seen_masks:
            expected_cache_hits += 1
        else:
            seen_masks.add(row_mask)

    assert fit.factorization_cache_size == len(unique_masks)
    assert fit.factorization_cache_hit_count == expected_cache_hits
    assert fit.factorization_cache_hit_count > 0


def test_rank_deficient_feature_subset_uses_explicit_non_estimable_outputs() -> None:
    matrix, design, blocks, contrasts = _fixture_inputs(
        "fixture_d_feature_level_failures"
    )
    summary = _summary("fixture_d_feature_level_failures")

    fit = fit_compound_symmetry_gls(
        matrix.to_numpy(dtype=np.float64),
        design.to_numpy(dtype=np.float64),
        blocks,
        float(summary["consensus_correlation"]),
        feature_ids=tuple(str(value) for value in matrix.index.tolist()),
        coefficient_names=tuple(str(value) for value in design.columns.tolist()),
        contrasts=contrasts.to_numpy(dtype=np.float64),
        contrast_names=tuple(str(value) for value in contrasts.columns.tolist()),
    )

    position = fit.feature_ids.index("D_rank_loss_only_A")
    coefficient_position = fit.coefficient_names.index("B")
    assert fit.feature_fit_statuses[position] == (
        COMPOUND_SYMMETRY_GLS_STATUS_PARTIAL_RANK_LOSS
    )
    assert not bool(fit.coefficient_estimability[coefficient_position, position])
    assert math.isnan(float(fit.coefficients[coefficient_position, position]))
    assert fit.contrast_coefficients is not None
    assert math.isnan(float(fit.contrast_coefficients[0, position]))
    assert fit.residual_degrees_of_freedom[position] == pytest.approx(4.0)


def test_block_design_columns_are_rejected_for_gls_correlation_model() -> None:
    with pytest.raises(PhosPyInputError, match="exclude fixed block dummy columns"):
        fit_compound_symmetry_gls(
            np.array([[1.0, 2.0, 3.0, 4.0]], dtype=np.float64),
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [1.0, 0.0, 1.0],
                    [0.0, 1.0, 1.0],
                ],
                dtype=np.float64,
            ),
            ("pair_1", "pair_1", "pair_2", "pair_2"),
            0.1,
            coefficient_names=("A", "B", "block[pair_2]"),
        )


def _fixture_inputs(
    fixture_id: str,
) -> tuple[pd.DataFrame, pd.DataFrame, tuple[str, ...], pd.DataFrame]:
    base = FIXTURE_ROOT / fixture_id
    matrix = pd.read_csv(base / "matrix.csv", index_col=0)
    design = pd.read_csv(base / "design.csv", index_col=0)
    blocks = tuple(
        str(value) for value in pd.read_csv(base / "blocks.csv")["block_id"].tolist()
    )
    contrasts = pd.read_csv(base / "contrasts.csv", index_col=0)
    return matrix, design, blocks, contrasts


def _summary(fixture_id: str) -> dict[str, str]:
    frame = pd.read_csv(FIXTURE_ROOT / fixture_id / "duplicate_correlation_summary.csv")
    return {str(row.field): str(row.value) for row in frame.itertuples(index=False)}


def _numeric_csv(fixture_id: str, filename: str, *, index_col: str) -> pd.DataFrame:
    return pd.read_csv(FIXTURE_ROOT / fixture_id / filename, index_col=index_col)


def _coefficient_frame(
    values: npt.NDArray[np.float64],
    fit: gls_module.CompoundSymmetryGLSFit,
) -> pd.DataFrame:
    return pd.DataFrame(
        values.T,
        index=pd.Index(fit.feature_ids, name="feature_id"),
        columns=pd.Index(fit.coefficient_names),
    )


def _contrast_feature_frame(
    values: npt.NDArray[np.float64],
    fit: gls_module.CompoundSymmetryGLSFit,
) -> pd.DataFrame:
    return pd.DataFrame(
        values.T,
        index=pd.Index(fit.feature_ids, name="feature_id"),
        columns=pd.Index(fit.contrast_names),
    )


def _covariance_frame(
    values: npt.NDArray[np.float64],
    names: tuple[str, ...],
) -> pd.DataFrame:
    return pd.DataFrame(values, index=pd.Index(names), columns=pd.Index(names))


def _assert_frame_close(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    absolute_tolerance: float,
) -> None:
    assert observed.index.tolist() == expected.index.astype(str).tolist()
    assert (
        observed.columns.astype(str).tolist() == expected.columns.astype(str).tolist()
    )
    np.testing.assert_allclose(
        observed.to_numpy(dtype=np.float64),
        expected.to_numpy(dtype=np.float64),
        rtol=GLS_PARITY_RELATIVE_TOLERANCE,
        atol=absolute_tolerance,
        equal_nan=True,
    )


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
        np.testing.assert_allclose(
            observed[column_name].to_numpy(dtype=np.float64),
            expected[column_name].to_numpy(dtype=np.float64),
            rtol=GLS_PARITY_RELATIVE_TOLERANCE,
            atol=absolute_tolerance,
            equal_nan=True,
        )


def _assert_gls_numeric_fit_close(
    observed: gls_module.CompoundSymmetryGLSFit,
    expected: gls_module.CompoundSymmetryGLSFit,
) -> None:
    assert observed.feature_ids == expected.feature_ids
    assert observed.coefficient_names == expected.coefficient_names
    assert observed.contrast_names == expected.contrast_names
    for attribute_name in (
        "coefficients",
        "stdev_unscaled",
        "residual_standard_deviation",
        "residual_variance",
        "residual_degrees_of_freedom",
        "residual_sum_of_squares",
        "average_expression",
        "coefficient_covariance",
        "feature_coefficient_covariances",
    ):
        np.testing.assert_allclose(
            getattr(observed, attribute_name),
            getattr(expected, attribute_name),
            rtol=1.0e-10,
            atol=STRICT_ABSOLUTE_TOLERANCE,
            equal_nan=True,
        )
    if observed.contrast_coefficients is not None:
        assert expected.contrast_coefficients is not None
        np.testing.assert_allclose(
            observed.contrast_coefficients,
            expected.contrast_coefficients,
            rtol=1.0e-10,
            atol=STRICT_ABSOLUTE_TOLERANCE,
            equal_nan=True,
        )
    if observed.contrast_stdev_unscaled is not None:
        assert expected.contrast_stdev_unscaled is not None
        np.testing.assert_allclose(
            observed.contrast_stdev_unscaled,
            expected.contrast_stdev_unscaled,
            rtol=1.0e-10,
            atol=STRICT_ABSOLUTE_TOLERANCE,
            equal_nan=True,
        )
