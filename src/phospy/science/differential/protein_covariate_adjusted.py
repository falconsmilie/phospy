"""Grouped protein-covariate-adjusted differential computation."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import cast

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import stats

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.empirical_bayes import (
    EmpiricalBayesFit,
    fit_empirical_bayes,
)
from phospy.science.differential.linear_model import (
    DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER,
    DifferentialDesignDecomposition,
    DifferentialDesignDecompositionError,
    decompose_differential_design,
)
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.models.empirical_bayes_config import (
    EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH,
)
from phospy.science.differential.models.protein_aware import (
    PROTEIN_AWARE_DIFFERENTIAL_AUGMENTED_DESIGN_FAILURE_DIAGNOSTIC_COLUMNS,
    PROTEIN_AWARE_DIFFERENTIAL_SITE_FAILURE_DIAGNOSTIC_COLUMNS,
    ProteinAwareDifferentialComputationRequest,
    ProteinAwareDifferentialComputationResult,
)
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_FIT_QUANTITIES_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_POSITIVE,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID,
)
from phospy.science.differential.quantification_depth import (
    QUANTIFICATION_DEPTH_LOG2_TREND_COVARIATE_NAME,
    QUANTIFICATION_DEPTH_TREND_TRANSFORMATION,
    log2_quantification_depth_series,
    validate_quantification_depth_series,
)
from phospy.science.statistics.multiple_testing import adjust_p_values

PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME = "protein_covariate"
PROTEIN_AWARE_CENTERING_POLICY = "mean_centered_no_standardization"

_FloatArray = npt.NDArray[np.float64]
_PROTEIN_VARIANCE_RELATIVE_TOLERANCE = float(np.finfo(np.float64).eps ** 0.75)


@dataclass(frozen=True, slots=True)
class _AlignedInputs:
    matrix: pd.DataFrame
    base_design: pd.DataFrame
    base_contrasts: pd.DataFrame
    protein_covariates: pd.DataFrame


@dataclass(frozen=True, slots=True)
class _ProteinCovariateSummary:
    raw_mean: float
    raw_standard_deviation: float
    centered_variance: float
    centered_vector: _FloatArray


@dataclass(frozen=True, slots=True)
class _Failure:
    total_protein_row_key: str
    status: str
    reason: str
    failure_message: str
    sample_count: int
    coefficient_count: int
    rank: int | None
    residual_degrees_of_freedom: float | None
    condition_number: float | None
    max_condition_number: float
    protein_raw_mean: float | None
    protein_raw_standard_deviation: float | None
    protein_centered_variance: float | None


@dataclass(frozen=True, slots=True)
class _ScaledSvdDiagnostics:
    rank: int | None
    residual_degrees_of_freedom: float | None
    condition_number: float | None


@dataclass(frozen=True, slots=True)
class _PostFitValidation:
    failed_positions: frozenset[int]
    failure_rows: tuple[dict[str, object], ...]


class ProteinCovariateAdjustedDifferentialKernel:
    """Run the version-1 protein-covariate-adjusted OLS kernel."""

    def run(
        self,
        request: ProteinAwareDifferentialComputationRequest,
    ) -> ProteinAwareDifferentialComputationResult:
        aligned = _align_inputs(request)
        sample_order = tuple(str(value) for value in request.sample_order)
        base_design = aligned.base_design
        base_contrasts = aligned.base_contrasts
        if PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME in base_design.columns:
            raise PhosPyInputError(
                "protein-aware differential base_design must not already contain "
                f"{PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME!r}"
            )

        site_ids = tuple(str(value) for value in aligned.matrix.index.tolist())
        total_protein_row_keys = tuple(
            str(value)
            for value in request.matched_pairs.loc[:, "total_protein_row_key"].tolist()
        )
        protein_identifiers = tuple(
            str(value)
            for value in request.matched_pairs.loc[:, "protein_identifier"].tolist()
        )
        contrast_names = tuple(str(value) for value in base_contrasts.columns.tolist())
        coefficient_names = (
            *tuple(str(value) for value in base_design.columns.tolist()),
            PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
        )

        residual_variance_by_position: dict[int, float] = {}
        mean_intensity_by_position: dict[int, float] = {}
        protein_coefficient_by_position: dict[int, float] = {}
        coefficients_by_position: dict[int, _FloatArray] = {}
        residuals_by_position: dict[int, _FloatArray] = {}
        contrast_effects_by_position: dict[int, _FloatArray] = {}
        contrast_scale_by_position: dict[int, _FloatArray] = {}
        successful_group_rows_by_key: dict[str, dict[str, object]] = {}
        successful_group_residual_dof_by_key: dict[str, float] = {}
        site_failure_rows: list[dict[str, object]] = []
        group_failure_rows: list[dict[str, object]] = []

        positions_by_row_key = _positions_by_total_protein_row_key(
            total_protein_row_keys
        )
        for row_key, positions in positions_by_row_key.items():
            protein_vector = (
                aligned.protein_covariates.reindex(
                    index=[row_key],
                    columns=list(sample_order),
                )
                .iloc[0]
                .to_numpy(dtype=float)
            )
            summary = _summarize_protein_covariate(protein_vector)
            failure = _protein_covariate_failure(
                row_key=row_key,
                summary=summary,
                sample_count=int(base_design.shape[0]),
                coefficient_count=int(base_design.shape[1]) + 1,
            )
            if failure is None:
                augmented_design = _augmented_design(
                    base_design=base_design,
                    centered_protein=summary.centered_vector,
                )
                augmented_contrasts = _augmented_contrasts(base_contrasts)
                decomposition, failure = _decompose_augmented_design(
                    row_key=row_key,
                    augmented_design=augmented_design,
                    augmented_contrasts=augmented_contrasts,
                    summary=summary,
                )
            else:
                decomposition = None
                augmented_design = None
                augmented_contrasts = None

            if failure is not None:
                group_failure_rows.append(_group_failure_row(failure))
                site_failure_rows.extend(
                    _site_failure_row(
                        site_id=site_ids[position],
                        row_key=row_key,
                        failure=failure,
                    )
                    for position in positions
                )
                continue

            if (
                decomposition is None
                or augmented_design is None
                or augmented_contrasts is None
            ):
                raise PhosPyInputError(
                    "protein-aware differential augmented design state is incomplete"
                )
            _fit_group(
                matrix=aligned.matrix,
                positions=positions,
                decomposition=decomposition,
                augmented_contrasts=augmented_contrasts,
                residual_variance_by_position=residual_variance_by_position,
                mean_intensity_by_position=mean_intensity_by_position,
                protein_coefficient_by_position=protein_coefficient_by_position,
                coefficients_by_position=coefficients_by_position,
                residuals_by_position=residuals_by_position,
                contrast_effects_by_position=contrast_effects_by_position,
                contrast_scale_by_position=contrast_scale_by_position,
            )
            successful_group_rows_by_key[row_key] = _successful_group_row(
                row_key=row_key,
                decomposition=decomposition,
                contrast_names=contrast_names,
            )
            successful_group_residual_dof_by_key[row_key] = float(
                decomposition.residual_degrees_of_freedom
            )

        post_fit_validation = _post_fit_validation(
            site_ids=site_ids,
            total_protein_row_keys=total_protein_row_keys,
            residual_variance_by_position=residual_variance_by_position,
            mean_intensity_by_position=mean_intensity_by_position,
            protein_coefficient_by_position=protein_coefficient_by_position,
            coefficients_by_position=coefficients_by_position,
            residuals_by_position=residuals_by_position,
            contrast_effects_by_position=contrast_effects_by_position,
            contrast_scale_by_position=contrast_scale_by_position,
        )
        site_failure_rows.extend(post_fit_validation.failure_rows)
        tested_positions = tuple(
            position
            for position in range(len(site_ids))
            if position in residual_variance_by_position
            and position not in post_fit_validation.failed_positions
        )
        if not tested_positions:
            failure_report = _failure_report(
                site_failure_rows=site_failure_rows,
                group_failure_rows=group_failure_rows,
            )
            if residual_variance_by_position:
                raise PhosPyInputError(
                    "protein-aware differential computation produced no successfully "
                    "tested sites after post-fit numerical eligibility; "
                    f"reason_counts={_reason_count_text(site_failure_rows)}",
                    diagnostics=failure_report,
                )
            raise PhosPyInputError(
                "protein-aware differential computation produced no successfully "
                "tested sites after protein covariate and augmented-design checks; "
                f"reason_counts={_reason_count_text(site_failure_rows)}",
                diagnostics=failure_report,
            )

        residual_dof = _require_common_residual_dof(
            [
                successful_group_residual_dof_by_key[row_key]
                for row_key in dict.fromkeys(
                    total_protein_row_keys[position] for position in tested_positions
                )
            ]
        )
        successful_group_rows = list(successful_group_rows_by_key.values())
        tested_index = pd.Index(
            [site_ids[position] for position in tested_positions],
            name="site_key",
        )
        residual_variance = np.asarray(
            [residual_variance_by_position[position] for position in tested_positions],
            dtype=np.float64,
        )
        mean_intensity = np.asarray(
            [mean_intensity_by_position[position] for position in tested_positions],
            dtype=np.float64,
        )
        quantification_depth = _quantification_depth_for_tested_sites(
            request=request,
            tested_index=tested_index,
        )
        trend_covariate = _resolve_empirical_bayes_trend_covariate(
            mean_intensity=mean_intensity,
            quantification_depth=quantification_depth,
        )

        try:
            eb_fit = fit_empirical_bayes(
                variances=residual_variance,
                residual_dof=residual_dof,
                method=request.empirical_bayes.method,
                trend=request.empirical_bayes.trend,
                winsor_tail_p=request.empirical_bayes.winsor_tail_p,
                trend_covariate=trend_covariate,
            )
        except ValueError as error:
            raise PhosPyInputError(
                "empirical-Bayes prior estimation failed for protein-aware "
                "differential analysis"
            ) from error

        posterior_variance, moderated_dof = _moderated_variance_and_dof(
            residual_variance=residual_variance,
            residual_dof=residual_dof,
            prior_variance=eb_fit.prior_variance,
            prior_dof=eb_fit.prior_degrees_of_freedom,
            tested_index=tested_index,
        )
        contrast_effects = np.vstack(
            [contrast_effects_by_position[position] for position in tested_positions]
        )
        contrast_scales = np.vstack(
            [contrast_scale_by_position[position] for position in tested_positions]
        )
        contrast_tables = _contrast_tables(
            contrast_effects=contrast_effects,
            contrast_scales=contrast_scales,
            posterior_variance=posterior_variance,
            moderated_dof=moderated_dof,
            contrast_names=contrast_names,
            tested_index=tested_index,
            multiple_testing_method=request.multiple_testing_method,
        )
        residual_variance_series = pd.Series(
            residual_variance.astype(float),
            index=tested_index.copy(),
            name="residual_variance",
        )
        posterior_variance_series = pd.Series(
            posterior_variance.astype(float),
            index=tested_index.copy(),
            name="posterior_residual_variance",
        )
        prior_variance_series = pd.Series(
            eb_fit.prior_variance.astype(float),
            index=tested_index.copy(),
            name="prior_residual_variance",
        )
        prior_dof_series = pd.Series(
            eb_fit.prior_degrees_of_freedom.astype(float),
            index=tested_index.copy(),
            name="prior_degrees_of_freedom",
        )
        protein_coefficient = pd.Series(
            [
                protein_coefficient_by_position[position]
                for position in tested_positions
            ],
            index=tested_index.copy(),
            name="protein_coefficient",
            dtype=float,
        )
        prior_diagnostics = EmpiricalBayesPriorDiagnostics(
            method=request.empirical_bayes.method,
            robust=request.empirical_bayes.method == "robust",
            trend=request.empirical_bayes.trend,
            winsor_tail_p=request.empirical_bayes.winsor_tail_p,
            base_prior_variance=eb_fit.base_prior_variance,
            base_prior_degrees_of_freedom=eb_fit.base_prior_degrees_of_freedom,
            robust_outlier_count=eb_fit.robust_outlier_count,
            robust_outlier_fraction=eb_fit.robust_outlier_fraction,
            winsorized_low_count=eb_fit.winsorized_low_count,
            winsorized_high_count=eb_fit.winsorized_high_count,
            prior_variance=prior_variance_series,
            prior_degrees_of_freedom=prior_dof_series,
            _assume_owned=True,
        )
        trend_diagnostics = _trend_diagnostics(
            request=request,
            eb_fit=eb_fit,
            tested_index=tested_index,
            mean_intensity=mean_intensity,
            enabled=request.empirical_bayes.trend,
            quantification_depth=quantification_depth,
        )

        return ProteinAwareDifferentialComputationResult(
            residual_variance=residual_variance_series,
            posterior_residual_variance=posterior_variance_series,
            prior_residual_variance=prior_variance_series,
            prior_degrees_of_freedom_series_value=prior_dof_series,
            prior_variance=float(np.nanmedian(eb_fit.prior_variance)),
            prior_degrees_of_freedom=float(
                np.nanmedian(eb_fit.prior_degrees_of_freedom)
            ),
            residual_degrees_of_freedom=float(residual_dof),
            empirical_bayes_method=request.empirical_bayes.method,
            empirical_bayes_robust=request.empirical_bayes.method == "robust",
            empirical_bayes_trend=request.empirical_bayes.trend,
            prior_diagnostics=prior_diagnostics,
            mean_variance_trend_diagnostics=trend_diagnostics,
            contrast_tables=contrast_tables,
            protein_coefficient=protein_coefficient,
            site_diagnostics=_site_diagnostics(
                tested_positions=tested_positions,
                site_ids=site_ids,
                protein_identifiers=protein_identifiers,
                total_protein_row_keys=total_protein_row_keys,
                protein_coefficient_by_position=protein_coefficient_by_position,
                mean_intensity_by_position=mean_intensity_by_position,
            ),
            augmented_design_diagnostics=_augmented_design_diagnostics(
                successful_group_rows
            ),
            tested_site_ids=tuple(tested_index.tolist()),
            coefficient_table=_coefficient_table(
                tested_positions=tested_positions,
                tested_index=tested_index,
                coefficient_names=coefficient_names,
                coefficients_by_position=coefficients_by_position,
            ),
            residuals=_residual_table(
                tested_positions=tested_positions,
                tested_index=tested_index,
                sample_order=sample_order,
                residuals_by_position=residuals_by_position,
            ),
            contrast_standard_error_scale=pd.DataFrame(
                contrast_scales,
                index=tested_index.copy(),
                columns=pd.Index(contrast_names, name=base_contrasts.columns.name),
            ),
            site_failure_diagnostics=_site_failure_diagnostics(site_failure_rows),
            augmented_design_failure_diagnostics=(
                _augmented_design_failure_diagnostics(group_failure_rows)
            ),
            method_id=request.method_id,
            _assume_owned=True,
        )


def run_protein_covariate_adjusted_differential(
    request: ProteinAwareDifferentialComputationRequest,
) -> ProteinAwareDifferentialComputationResult:
    """Run the private protein-covariate-adjusted computation kernel."""

    return ProteinCovariateAdjustedDifferentialKernel().run(request)


def _positions_by_total_protein_row_key(
    total_protein_row_keys: tuple[str, ...],
) -> dict[str, tuple[int, ...]]:
    positions_by_row_key: dict[str, list[int]] = {}
    for position, row_key in enumerate(total_protein_row_keys):
        positions_by_row_key.setdefault(row_key, []).append(position)
    return {
        row_key: tuple(positions) for row_key, positions in positions_by_row_key.items()
    }


def _align_inputs(
    request: ProteinAwareDifferentialComputationRequest,
) -> _AlignedInputs:
    base_design_model = request.base_design
    base_contrasts_model = request.base_contrasts
    if not isinstance(base_design_model, DesignMatrix) or not isinstance(
        base_contrasts_model,
        ContrastMatrix,
    ):
        raise PhosPyInputError(
            "protein-aware differential request must contain normalized base design "
            "and contrast matrices"
        )
    sample_order = list(request.sample_order)
    base_design = base_design_model.frame.loc[sample_order, :]
    return _AlignedInputs(
        matrix=request.phosphosite_matrix.loc[:, sample_order],
        base_design=base_design,
        base_contrasts=base_contrasts_model.frame.loc[list(base_design.columns), :],
        protein_covariates=request.resolved_protein_covariates.loc[:, sample_order],
    )


def _summarize_protein_covariate(
    protein_vector: _FloatArray,
) -> _ProteinCovariateSummary:
    vector = np.asarray(protein_vector, dtype=np.float64)
    if not np.isfinite(vector).all():
        return _ProteinCovariateSummary(
            raw_mean=math.nan,
            raw_standard_deviation=math.nan,
            centered_variance=math.nan,
            centered_vector=np.full(vector.shape, np.nan, dtype=np.float64),
        )
    raw_mean = float(np.mean(vector))
    raw_standard_deviation = (
        float(np.std(vector, ddof=1)) if int(vector.size) > 1 else 0.0
    )
    centered = cast(_FloatArray, np.asarray(vector - raw_mean, dtype=np.float64))
    centered_variance = (
        float(np.var(centered, ddof=1)) if int(centered.size) > 1 else 0.0
    )
    return _ProteinCovariateSummary(
        raw_mean=raw_mean,
        raw_standard_deviation=raw_standard_deviation,
        centered_variance=centered_variance,
        centered_vector=centered,
    )


def _protein_covariate_failure(
    *,
    row_key: str,
    summary: _ProteinCovariateSummary,
    sample_count: int,
    coefficient_count: int,
) -> _Failure | None:
    if not np.isfinite(summary.centered_vector).all():
        return _failure(
            row_key=row_key,
            status=DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
            reason=DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
            message="protein covariate contains non-finite values",
            sample_count=sample_count,
            coefficient_count=coefficient_count,
            summary=summary,
        )
    centered_sum_squares = float(summary.centered_vector @ summary.centered_vector)
    scale = max(float(np.linalg.norm(summary.centered_vector)), 1.0)
    tolerance = (
        _PROTEIN_VARIANCE_RELATIVE_TOLERANCE
        * float(max(sample_count, coefficient_count))
        * scale
    )
    if (
        not math.isfinite(centered_sum_squares)
        or centered_sum_squares <= tolerance
        or summary.centered_variance <= 0.0
    ):
        return _failure(
            row_key=row_key,
            status=DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
            reason=DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
            message="protein covariate has zero or numerically unusable variance",
            sample_count=sample_count,
            coefficient_count=coefficient_count,
            summary=summary,
        )
    return None


def _augmented_design(
    *,
    base_design: pd.DataFrame,
    centered_protein: _FloatArray,
) -> pd.DataFrame:
    frame = base_design.copy(deep=True)
    frame.loc[:, PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME] = centered_protein
    return frame


def _augmented_contrasts(base_contrasts: pd.DataFrame) -> pd.DataFrame:
    protein_row = pd.DataFrame(
        np.zeros((1, int(base_contrasts.shape[1])), dtype=float),
        index=pd.Index(
            [PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME],
            name=base_contrasts.index.name,
        ),
        columns=base_contrasts.columns.copy(),
    )
    return pd.concat([base_contrasts.copy(deep=True), protein_row], axis=0)


def _decompose_augmented_design(
    *,
    row_key: str,
    augmented_design: pd.DataFrame,
    augmented_contrasts: pd.DataFrame,
    summary: _ProteinCovariateSummary,
) -> tuple[DifferentialDesignDecomposition | None, _Failure | None]:
    design_values = augmented_design.to_numpy(dtype=float)
    try:
        decomposition = decompose_differential_design(design_values)
    except DifferentialDesignDecompositionError as error:
        diagnostics = _scaled_svd_diagnostics(design_values)
        return None, _failure(
            row_key=row_key,
            status=DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
            reason=_decomposition_failure_reason(str(error)),
            message=str(error),
            sample_count=int(augmented_design.shape[0]),
            coefficient_count=int(augmented_design.shape[1]),
            summary=summary,
            rank=diagnostics.rank,
            residual_degrees_of_freedom=diagnostics.residual_degrees_of_freedom,
            condition_number=diagnostics.condition_number,
            max_condition_number=DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER,
        )

    invalid_contrast_positions = decomposition.invalid_contrast_positions(
        augmented_contrasts.to_numpy(dtype=float)
    )
    if invalid_contrast_positions:
        contrast_names = tuple(
            str(augmented_contrasts.columns[position])
            for position in invalid_contrast_positions
        )
        return None, _failure(
            row_key=row_key,
            status=DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
            reason=DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
            message=(
                "protein-aware augmented contrast is non-estimable or zero-length; "
                f"contrasts={contrast_names}"
            ),
            sample_count=decomposition.sample_count,
            coefficient_count=decomposition.coefficient_count,
            summary=summary,
            rank=decomposition.rank,
            residual_degrees_of_freedom=decomposition.residual_degrees_of_freedom,
            condition_number=decomposition.condition_number,
            max_condition_number=decomposition.max_condition_number,
        )
    return decomposition, None


def _fit_group(
    *,
    matrix: pd.DataFrame,
    positions: tuple[int, ...],
    decomposition: DifferentialDesignDecomposition,
    augmented_contrasts: pd.DataFrame,
    residual_variance_by_position: dict[int, float],
    mean_intensity_by_position: dict[int, float],
    protein_coefficient_by_position: dict[int, float],
    coefficients_by_position: dict[int, _FloatArray],
    residuals_by_position: dict[int, _FloatArray],
    contrast_effects_by_position: dict[int, _FloatArray],
    contrast_scale_by_position: dict[int, _FloatArray],
) -> None:
    response = matrix.iloc[list(positions), :].to_numpy(dtype=float).T
    try:
        linear_fit = decomposition.fit(response)
    except DifferentialDesignDecompositionError as error:
        raise PhosPyInputError(
            "protein-aware differential model fitting failed under an already "
            f"validated augmented decomposition: {error}"
        ) from error

    contrast_values = augmented_contrasts.to_numpy(dtype=float)
    contrast_effects = cast(
        _FloatArray,
        np.asarray(linear_fit.coefficients.T @ contrast_values, dtype=np.float64),
    )
    contrast_scales = decomposition.contrast_scales(contrast_values)
    for local_position, site_position in enumerate(positions):
        residual_variance_by_position[site_position] = float(
            linear_fit.residual_variance[local_position]
        )
        mean_intensity_by_position[site_position] = float(
            np.mean(response[:, local_position])
        )
        protein_coefficient_by_position[site_position] = float(
            linear_fit.coefficients[-1, local_position]
        )
        coefficients_by_position[site_position] = cast(
            _FloatArray,
            np.asarray(linear_fit.coefficients[:, local_position], dtype=np.float64),
        )
        residuals_by_position[site_position] = cast(
            _FloatArray,
            np.asarray(linear_fit.residuals[:, local_position], dtype=np.float64),
        )
        contrast_effects_by_position[site_position] = cast(
            _FloatArray,
            np.asarray(contrast_effects[local_position, :], dtype=np.float64),
        )
        contrast_scale_by_position[site_position] = cast(
            _FloatArray,
            np.asarray(contrast_scales, dtype=np.float64),
        )


def _post_fit_validation(
    *,
    site_ids: tuple[str, ...],
    total_protein_row_keys: tuple[str, ...],
    residual_variance_by_position: dict[int, float],
    mean_intensity_by_position: dict[int, float],
    protein_coefficient_by_position: dict[int, float],
    coefficients_by_position: dict[int, _FloatArray],
    residuals_by_position: dict[int, _FloatArray],
    contrast_effects_by_position: dict[int, _FloatArray],
    contrast_scale_by_position: dict[int, _FloatArray],
) -> _PostFitValidation:
    failed_positions: set[int] = set()
    rows: list[dict[str, object]] = []
    for position in sorted(residual_variance_by_position):
        row = _post_fit_failure_row(
            position=position,
            site_id=site_ids[position],
            row_key=total_protein_row_keys[position],
            residual_variance=residual_variance_by_position[position],
            mean_intensity=mean_intensity_by_position[position],
            protein_coefficient=protein_coefficient_by_position[position],
            coefficients=coefficients_by_position[position],
            residuals=residuals_by_position[position],
            contrast_effects=contrast_effects_by_position[position],
            contrast_scales=contrast_scale_by_position[position],
        )
        if row is None:
            continue
        failed_positions.add(position)
        rows.append(row)
    return _PostFitValidation(
        failed_positions=frozenset(failed_positions),
        failure_rows=tuple(rows),
    )


def _post_fit_failure_row(
    *,
    position: int,
    site_id: str,
    row_key: str,
    residual_variance: float,
    mean_intensity: float,
    protein_coefficient: float,
    coefficients: _FloatArray,
    residuals: _FloatArray,
    contrast_effects: _FloatArray,
    contrast_scales: _FloatArray,
) -> dict[str, object] | None:
    residual_variance_value = float(residual_variance)
    coefficient_non_finite_count = _non_finite_count(coefficients)
    residual_non_finite_count = _non_finite_count(residuals)
    contrast_effect_non_finite_count = _non_finite_count(contrast_effects)
    contrast_scale_non_finite_count = _non_finite_count(contrast_scales)
    contrast_scale_non_positive_count = _non_positive_count(contrast_scales)
    mean_intensity_value = float(mean_intensity)
    protein_coefficient_value = float(protein_coefficient)

    reason: str | None
    message: str
    if not math.isfinite(residual_variance_value):
        reason = DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_FINITE
        message = "protein-aware fitted residual variance is non-finite"
    elif residual_variance_value <= 0.0:
        reason = DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_POSITIVE
        message = "protein-aware fitted residual variance is not greater than zero"
    elif (
        coefficient_non_finite_count
        or residual_non_finite_count
        or contrast_effect_non_finite_count
        or contrast_scale_non_finite_count
        or contrast_scale_non_positive_count
        or not math.isfinite(mean_intensity_value)
        or not math.isfinite(protein_coefficient_value)
    ):
        reason = DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_FIT_QUANTITIES_NON_FINITE
        message = (
            "protein-aware fitted coefficients, residuals, mean intensity, nuisance "
            "coefficient, contrast effects, or contrast scales are non-finite or "
            "inadmissible"
        )
    else:
        return None

    return {
        "site_key": site_id,
        "total_protein_row_key": row_key,
        DIFFERENTIAL_RESULT_STATUS_COLUMN: (
            DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID
        ),
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: reason,
        "failure_message": message,
        "post_fit_validation_stage": "protein_aware_model_fit_numerical_eligibility",
        "raw_fit_position": int(position),
        "residual_variance": _finite_or_none(residual_variance_value),
        "mean_intensity": _finite_or_none(mean_intensity_value),
        "protein_coefficient": _finite_or_none(protein_coefficient_value),
        "coefficient_non_finite_count": coefficient_non_finite_count,
        "residual_non_finite_count": residual_non_finite_count,
        "contrast_effect_non_finite_count": contrast_effect_non_finite_count,
        "contrast_scale_non_finite_count": contrast_scale_non_finite_count,
        "contrast_scale_non_positive_count": contrast_scale_non_positive_count,
    }


def _non_finite_count(values: _FloatArray) -> int:
    numeric = np.asarray(values, dtype=np.float64)
    return int(np.count_nonzero(~np.isfinite(numeric)))


def _non_positive_count(values: _FloatArray) -> int:
    numeric = np.asarray(values, dtype=np.float64)
    return int(np.count_nonzero(np.isfinite(numeric) & (numeric <= 0.0)))


def _finite_or_none(value: float) -> float | None:
    return float(value) if math.isfinite(float(value)) else None


def _moderated_variance_and_dof(
    *,
    residual_variance: _FloatArray,
    residual_dof: float,
    prior_variance: _FloatArray,
    prior_dof: _FloatArray,
    tested_index: pd.Index,
) -> tuple[_FloatArray, _FloatArray]:
    total_residual_dof = residual_dof * float(tested_index.size)
    posterior_variance = np.empty_like(residual_variance, dtype=float)
    finite_prior_dof = np.isfinite(prior_dof)
    if np.any(finite_prior_dof):
        posterior_variance[finite_prior_dof] = (
            prior_dof[finite_prior_dof] * prior_variance[finite_prior_dof]
            + residual_dof * residual_variance[finite_prior_dof]
        ) / (prior_dof[finite_prior_dof] + residual_dof)
    if np.any(~finite_prior_dof):
        posterior_variance[~finite_prior_dof] = prior_variance[~finite_prior_dof]

    moderated_dof = cast(
        _FloatArray,
        np.asarray(
            np.where(
                finite_prior_dof,
                np.minimum(residual_dof + prior_dof, total_residual_dof),
                total_residual_dof,
            ),
            dtype=np.float64,
        ),
    )
    invalid = ~np.isfinite(moderated_dof) | (moderated_dof <= 0.0)
    if np.any(invalid):
        raise PhosPyInputError(
            "protein-aware differential analysis produced invalid moderated "
            "degrees of freedom; degrees of freedom must be finite and > 0.0; "
            f"{_preview_invalid_entries(invalid, moderated_dof, row_index=tested_index)}"
        )
    return cast(_FloatArray, posterior_variance), moderated_dof


def _quantification_depth_for_tested_sites(
    *,
    request: ProteinAwareDifferentialComputationRequest,
    tested_index: pd.Index,
) -> pd.Series | None:
    if (
        request.empirical_bayes.trend_covariate
        != EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH
    ):
        return None
    depth = request.quantification_depth
    if depth is None:
        raise PhosPyInputError(
            "protein_aware_differential_request.quantification_depth must be "
            "provided for quantification-depth empirical-Bayes variance trends"
        )
    try:
        selected = depth.loc[tested_index]
    except KeyError as exc:
        raise PhosPyInputError(
            "protein_aware_differential_request.quantification_depth.index must "
            "contain every protein-aware tested site_key"
        ) from exc
    if not selected.index.equals(tested_index):
        raise PhosPyInputError(
            "protein_aware_differential_request.quantification_depth.index must "
            "align to the protein-aware tested site_key order"
        )
    return validate_quantification_depth_series(
        selected,
        field_name="protein_aware_differential_request.quantification_depth",
        expected_index=tested_index,
    )


def _resolve_empirical_bayes_trend_covariate(
    *,
    mean_intensity: _FloatArray,
    quantification_depth: pd.Series | None,
) -> _FloatArray:
    if quantification_depth is None:
        return np.asarray(mean_intensity, dtype=np.float64)
    return np.asarray(
        log2_quantification_depth_series(
            quantification_depth,
            field_name="protein_aware_differential_request.quantification_depth",
            expected_index=quantification_depth.index,
        ).to_numpy(dtype=float),
        dtype=np.float64,
    )


def _contrast_tables(
    *,
    contrast_effects: _FloatArray,
    contrast_scales: _FloatArray,
    posterior_variance: _FloatArray,
    moderated_dof: _FloatArray,
    contrast_names: tuple[str, ...],
    tested_index: pd.Index,
    multiple_testing_method: str,
) -> dict[str, pd.DataFrame]:
    posterior_sd = np.sqrt(posterior_variance)
    contrast_tables: dict[str, pd.DataFrame] = {}
    for column_idx, contrast_name in enumerate(contrast_names):
        log_fc = contrast_effects[:, column_idx]
        standard_error = posterior_sd * contrast_scales[:, column_idx]
        invalid_standard_error = ~np.isfinite(standard_error) | (standard_error <= 0.0)
        if np.any(invalid_standard_error):
            raise PhosPyInputError(
                "protein-aware differential analysis produced unstable standard "
                f"errors for contrast {contrast_name!r}; standard errors must be "
                "finite and > 0.0; "
                f"{_preview_invalid_entries(invalid_standard_error, standard_error, row_index=tested_index)}"
            )
        moderated_t = log_fc / standard_error
        p_values = np.asarray(
            2.0 * stats.t.sf(np.abs(moderated_t), df=moderated_dof),
            dtype=np.float64,
        )
        invalid_p_values = ~np.isfinite(p_values) | (p_values < 0.0) | (p_values > 1.0)
        if np.any(invalid_p_values):
            raise PhosPyInputError(
                "protein-aware differential analysis produced invalid p-values for "
                f"contrast {contrast_name!r}; P.Value must be finite and within "
                "[0, 1]; "
                f"{_preview_invalid_entries(invalid_p_values, p_values, row_index=tested_index)}"
            )
        adjusted = adjust_p_values(
            p_values,
            method=multiple_testing_method,
        )
        contrast_tables[contrast_name] = pd.DataFrame(
            {
                "logFC": log_fc.astype(float),
                "t": moderated_t.astype(float),
                "P.Value": p_values.astype(float),
                "adj.P.Val": adjusted.astype(float),
            },
            index=tested_index.copy(),
        )
    return contrast_tables


def _trend_diagnostics(
    *,
    request: ProteinAwareDifferentialComputationRequest,
    eb_fit: EmpiricalBayesFit,
    tested_index: pd.Index,
    mean_intensity: _FloatArray,
    enabled: bool,
    quantification_depth: pd.Series | None,
) -> MeanVarianceTrendDiagnostics | None:
    if not enabled:
        return None
    trend_covariate = eb_fit.trend_covariate
    log_residual_variance = eb_fit.log_residual_variance
    fitted_log_prior_variance = eb_fit.fitted_log_prior_variance
    if (
        trend_covariate is None
        or log_residual_variance is None
        or fitted_log_prior_variance is None
    ):
        raise PhosPyInputError(
            "protein-aware differential empirical-Bayes trend diagnostics are "
            "incomplete"
        )
    mean_intensity_series = pd.Series(
        mean_intensity,
        index=tested_index.copy(),
        name="mean_intensity",
    )
    if quantification_depth is not None:
        return MeanVarianceTrendDiagnostics(
            mean_intensity=mean_intensity_series,
            trend_covariate=pd.Series(
                trend_covariate,
                index=tested_index.copy(),
                name=QUANTIFICATION_DEPTH_LOG2_TREND_COVARIATE_NAME,
            ),
            trend_covariate_name=request.empirical_bayes.trend_covariate,
            trend_covariate_transformation=QUANTIFICATION_DEPTH_TREND_TRANSFORMATION,
            quantification_depth=quantification_depth,
            quantification_depth_kind=request.empirical_bayes.quantification_depth_kind,
            log_residual_variance=pd.Series(
                log_residual_variance,
                index=tested_index.copy(),
                name="log_residual_variance",
            ),
            fitted_log_prior_variance=pd.Series(
                fitted_log_prior_variance,
                index=tested_index.copy(),
                name="fitted_log_prior_variance",
            ),
            _assume_owned=True,
        )
    return MeanVarianceTrendDiagnostics(
        mean_intensity=mean_intensity_series,
        log_residual_variance=pd.Series(
            log_residual_variance,
            index=tested_index.copy(),
            name="log_residual_variance",
        ),
        fitted_log_prior_variance=pd.Series(
            fitted_log_prior_variance,
            index=tested_index.copy(),
            name="fitted_log_prior_variance",
        ),
        _assume_owned=True,
    )


def _site_diagnostics(
    *,
    tested_positions: tuple[int, ...],
    site_ids: tuple[str, ...],
    protein_identifiers: tuple[str, ...],
    total_protein_row_keys: tuple[str, ...],
    protein_coefficient_by_position: dict[int, float],
    mean_intensity_by_position: dict[int, float],
) -> pd.DataFrame:
    index = pd.Index(
        [site_ids[position] for position in tested_positions], name="site_key"
    )
    return pd.DataFrame(
        {
            "site_key": [site_ids[position] for position in tested_positions],
            "protein_identifier": [
                protein_identifiers[position] for position in tested_positions
            ],
            "total_protein_row_key": [
                total_protein_row_keys[position] for position in tested_positions
            ],
            "protein_coefficient": [
                protein_coefficient_by_position[position]
                for position in tested_positions
            ],
            "mean_intensity": [
                mean_intensity_by_position[position] for position in tested_positions
            ],
        },
        index=index.copy(),
    )


def _coefficient_table(
    *,
    tested_positions: tuple[int, ...],
    tested_index: pd.Index,
    coefficient_names: tuple[str, ...],
    coefficients_by_position: dict[int, _FloatArray],
) -> pd.DataFrame:
    return pd.DataFrame(
        np.vstack(
            [coefficients_by_position[position] for position in tested_positions]
        ),
        index=tested_index.copy(),
        columns=pd.Index(coefficient_names, name="coefficient"),
    )


def _residual_table(
    *,
    tested_positions: tuple[int, ...],
    tested_index: pd.Index,
    sample_order: tuple[str, ...],
    residuals_by_position: dict[int, _FloatArray],
) -> pd.DataFrame:
    return pd.DataFrame(
        np.vstack([residuals_by_position[position] for position in tested_positions]),
        index=tested_index.copy(),
        columns=pd.Index(sample_order, name="sample_id"),
    )


def _successful_group_row(
    *,
    row_key: str,
    decomposition: DifferentialDesignDecomposition,
    contrast_names: tuple[str, ...],
) -> dict[str, object]:
    return {
        "total_protein_row_key": row_key,
        "sample_count": decomposition.sample_count,
        "coefficient_count": decomposition.coefficient_count,
        "rank": decomposition.rank,
        "residual_degrees_of_freedom": decomposition.residual_degrees_of_freedom,
        "condition_number": decomposition.condition_number,
        "max_condition_number": decomposition.max_condition_number,
        "decomposition_method": decomposition.decomposition_method,
        "solver": decomposition.solver,
        "column_scale_method": decomposition.column_scale_method,
        "rank_tolerance": decomposition.rank_tolerance,
        "rank_tolerance_policy": decomposition.rank_tolerance_policy,
        "singular_values": decomposition.singular_values,
        "contrast_names": contrast_names,
        "protein_covariate_contrast_weights": tuple(0.0 for _ in contrast_names),
        "protein_covariate_policy": PROTEIN_AWARE_CENTERING_POLICY,
    }


def _augmented_design_diagnostics(
    rows: list[dict[str, object]],
) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    frame.index = pd.Index(
        [str(row["total_protein_row_key"]) for row in rows],
        name="total_protein_row_key",
    )
    return frame


def _site_failure_row(
    *,
    site_id: str,
    row_key: str,
    failure: _Failure,
) -> dict[str, object]:
    return {
        "site_key": site_id,
        "total_protein_row_key": row_key,
        DIFFERENTIAL_RESULT_STATUS_COLUMN: failure.status,
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: failure.reason,
        "failure_message": failure.failure_message,
    }


def _group_failure_row(failure: _Failure) -> dict[str, object]:
    return {
        "total_protein_row_key": failure.total_protein_row_key,
        DIFFERENTIAL_RESULT_STATUS_COLUMN: failure.status,
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: failure.reason,
        "sample_count": failure.sample_count,
        "coefficient_count": failure.coefficient_count,
        "rank": failure.rank,
        "residual_degrees_of_freedom": failure.residual_degrees_of_freedom,
        "condition_number": failure.condition_number,
        "max_condition_number": failure.max_condition_number,
        "protein_raw_mean": failure.protein_raw_mean,
        "protein_raw_standard_deviation": failure.protein_raw_standard_deviation,
        "protein_centered_variance": failure.protein_centered_variance,
        "failure_message": failure.failure_message,
    }


def _site_failure_diagnostics(rows: list[dict[str, object]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(
            columns=pd.Index(
                PROTEIN_AWARE_DIFFERENTIAL_SITE_FAILURE_DIAGNOSTIC_COLUMNS,
                dtype=object,
            ),
            index=pd.Index((), name="site_key"),
        )
    frame = pd.DataFrame(rows)
    frame.index = pd.Index(frame["site_key"].astype(str).tolist(), name="site_key")
    return frame


def _augmented_design_failure_diagnostics(
    rows: list[dict[str, object]],
) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(
            columns=pd.Index(
                PROTEIN_AWARE_DIFFERENTIAL_AUGMENTED_DESIGN_FAILURE_DIAGNOSTIC_COLUMNS,
                dtype=object,
            ),
            index=pd.Index((), name="total_protein_row_key"),
        )
    frame = pd.DataFrame(rows)
    frame.index = pd.Index(
        frame["total_protein_row_key"].astype(str).tolist(),
        name="total_protein_row_key",
    )
    return frame


def _failure_report(
    *,
    site_failure_rows: list[dict[str, object]],
    group_failure_rows: list[dict[str, object]],
) -> dict[str, object]:
    return {
        "site_failure_diagnostics": _site_failure_diagnostics(site_failure_rows),
        "augmented_design_failure_diagnostics": (
            _augmented_design_failure_diagnostics(group_failure_rows)
        ),
        "status_counts": _failure_count_mapping(
            rows=site_failure_rows,
            column_name=DIFFERENTIAL_RESULT_STATUS_COLUMN,
        ),
        "reason_counts": _failure_count_mapping(
            rows=site_failure_rows,
            column_name=DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        ),
    }


def _failure(
    *,
    row_key: str,
    status: str,
    reason: str,
    message: str,
    sample_count: int,
    coefficient_count: int,
    summary: _ProteinCovariateSummary,
    rank: int | None = None,
    residual_degrees_of_freedom: float | None = None,
    condition_number: float | None = None,
    max_condition_number: float = DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER,
) -> _Failure:
    return _Failure(
        total_protein_row_key=row_key,
        status=status,
        reason=reason,
        failure_message=message,
        sample_count=int(sample_count),
        coefficient_count=int(coefficient_count),
        rank=rank,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
        condition_number=condition_number,
        max_condition_number=float(max_condition_number),
        protein_raw_mean=(
            summary.raw_mean if math.isfinite(summary.raw_mean) else None
        ),
        protein_raw_standard_deviation=(
            summary.raw_standard_deviation
            if math.isfinite(summary.raw_standard_deviation)
            else None
        ),
        protein_centered_variance=(
            summary.centered_variance
            if math.isfinite(summary.centered_variance)
            else None
        ),
    )


def _decomposition_failure_reason(message: str) -> str:
    if "residual degrees of freedom must be positive" in message:
        return DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF
    if "rank deficient" in message:
        return DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
    if "too ill-conditioned" in message:
        return DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED
    return DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED


def _scaled_svd_diagnostics(design_values: _FloatArray) -> _ScaledSvdDiagnostics:
    sample_count = int(design_values.shape[0])
    coefficient_count = int(design_values.shape[1])
    try:
        column_scales = np.asarray(np.linalg.norm(design_values, axis=0), dtype=float)
        usable = np.isfinite(column_scales) & (column_scales > 0.0)
        if not np.any(usable):
            return _ScaledSvdDiagnostics(
                rank=0,
                residual_degrees_of_freedom=float(sample_count),
                condition_number=math.inf,
            )
        scaled = design_values[:, usable] / column_scales[usable][np.newaxis, :]
        singular_values = np.asarray(
            np.linalg.svd(scaled, compute_uv=False),
            dtype=np.float64,
        )
    except np.linalg.LinAlgError:
        return _ScaledSvdDiagnostics(
            rank=None,
            residual_degrees_of_freedom=None,
            condition_number=None,
        )
    if singular_values.size == 0 or not np.isfinite(singular_values).all():
        return _ScaledSvdDiagnostics(
            rank=None,
            residual_degrees_of_freedom=None,
            condition_number=None,
        )
    largest = float(singular_values[0])
    tolerance = (
        np.finfo(np.float64).eps * float(max(sample_count, coefficient_count)) * largest
    )
    rank = int(np.count_nonzero(singular_values > tolerance))
    smallest = float(singular_values[-1])
    condition_number = math.inf if smallest == 0.0 else float(largest / smallest)
    return _ScaledSvdDiagnostics(
        rank=rank,
        residual_degrees_of_freedom=float(sample_count - rank),
        condition_number=condition_number,
    )


def _require_common_residual_dof(values: list[float]) -> float:
    if not values:
        raise PhosPyInputError(
            "protein-aware differential computation requires at least one "
            "successful augmented design"
        )
    first = float(values[0])
    if not math.isfinite(first) or first <= 0.0:
        raise PhosPyInputError(
            "protein-aware differential successful augmented designs must have "
            "positive residual degrees of freedom"
        )
    disagree = [
        value
        for value in values
        if not math.isclose(float(value), first, rel_tol=0.0, abs_tol=1.0e-12)
    ]
    if disagree:
        raise PhosPyInputError(
            "protein-aware differential successful augmented design groups produced "
            "inconsistent residual degrees of freedom"
        )
    return first


def _reason_count_text(rows: list[dict[str, object]]) -> str:
    counts = Counter(str(row[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]) for row in rows)
    if not counts:
        return "{}"
    return (
        "{"
        + ", ".join(f"{reason}: {counts[reason]}" for reason in sorted(counts))
        + "}"
    )


def _failure_count_mapping(
    *,
    rows: list[dict[str, object]],
    column_name: str,
) -> dict[str, int]:
    counts = Counter(str(row[column_name]) for row in rows if str(row[column_name]))
    return {key: int(counts[key]) for key in sorted(counts)}


def _preview_invalid_entries(
    invalid_mask: npt.NDArray[np.bool_],
    values: _FloatArray,
    *,
    row_index: pd.Index,
) -> str:
    invalid_positions = np.flatnonzero(invalid_mask)
    preview = ", ".join(
        f"({row_index[position]!r}, {values[position]:.6g})"
        for position in invalid_positions[:3]
    )
    suffix = (
        ""
        if invalid_positions.size <= 3
        else f", +{int(invalid_positions.size - 3)} more"
    )
    return (
        f"invalid values: {preview}{suffix}; "
        f"invalid_entry_count={int(invalid_positions.size)}"
    )


__all__ = [
    "PROTEIN_AWARE_CENTERING_POLICY",
    "PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME",
    "ProteinCovariateAdjustedDifferentialKernel",
    "run_protein_covariate_adjusted_differential",
]
