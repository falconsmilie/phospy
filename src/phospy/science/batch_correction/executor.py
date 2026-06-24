"""Deterministic SPS/RUV-style numerical correction executor."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Protocol, TypeAlias

import numpy as np
import pandas as pd

from phospy.contracts.configs.preprocessing import TemporaryImputationMethod
from phospy.errors.input import PhosPyInputError
from phospy.provenance import fingerprint_matrix
from phospy.provenance.models import TableFingerprint

SPS_RUV_STYLE_EXECUTOR_ID = "deterministic_sps_ruv_style_executor_v1"
SPS_RUV_STYLE_METHODS = frozenset(
    {"sps_ruv_style", "control_site_ruv_style", "ruv_iii_style"}
)

JsonScalar: TypeAlias = str | int | float | bool | None


class _ControlRowLike(Protocol):
    site_key: str
    row_position: int
    weight: float | None


class _ObservationMaskLike(Protocol):
    feature_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    originally_missing_cells: tuple[tuple[str, str], ...]

    def to_payload(self) -> dict[str, object]: ...


class _TemporaryImputationPolicyLike(Protocol):
    allowed: bool
    method: object
    method_parameters: tuple[tuple[str, JsonScalar], ...]

    def to_payload(self) -> dict[str, object]: ...


class _ResolvedPlanLike(Protocol):
    method: str
    resolved_design_matrix: pd.DataFrame
    condition_terms_to_preserve: tuple[str, ...]
    eligible_control_site_rows: tuple[_ControlRowLike, ...]
    observation_mask: _ObservationMaskLike
    temporary_imputation_policy: _TemporaryImputationPolicyLike
    n_unwanted_factors: int | None
    stage_order: tuple[str, ...]
    provenance_seed_data: Mapping[str, object]

    def to_payload(self) -> dict[str, object]: ...


@dataclass(frozen=True, slots=True)
class SpsRuvStyleExecutorDiagnostics:
    """Numerical diagnostics emitted by the SPS/RUV-style executor."""

    method: str
    executor_id: str
    status: str
    matrix_shape_before: tuple[int, int]
    matrix_shape_after: tuple[int, int]
    control_site_count: int
    protected_design_terms: tuple[str, ...]
    protected_design_rank: int
    requested_unwanted_factors: int | None
    estimated_unwanted_factors: int
    singular_values: tuple[float, ...]
    originally_missing_cell_count: int
    corrected_observed_cell_count: int
    withheld_cell_count: int
    max_abs_adjustment: float
    mean_abs_adjustment: float
    warnings: tuple[str, ...]

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible diagnostics payload."""

        return {
            "method": self.method,
            "executor_id": self.executor_id,
            "status": self.status,
            "matrix_shape_before": list(self.matrix_shape_before),
            "matrix_shape_after": list(self.matrix_shape_after),
            "control_site_count": self.control_site_count,
            "protected_design_terms": list(self.protected_design_terms),
            "protected_design_rank": self.protected_design_rank,
            "requested_unwanted_factors": self.requested_unwanted_factors,
            "estimated_unwanted_factors": self.estimated_unwanted_factors,
            "singular_values": list(self.singular_values),
            "originally_missing_cell_count": self.originally_missing_cell_count,
            "corrected_observed_cell_count": self.corrected_observed_cell_count,
            "withheld_cell_count": self.withheld_cell_count,
            "max_abs_adjustment": self.max_abs_adjustment,
            "mean_abs_adjustment": self.mean_abs_adjustment,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True, slots=True)
class SpsRuvStyleExecutorResult:
    """Output from deterministic SPS/RUV-style numerical correction."""

    corrected_matrix: pd.DataFrame
    estimated_unwanted_factors: pd.DataFrame
    diagnostics: SpsRuvStyleExecutorDiagnostics
    warnings: tuple[str, ...]
    withheld_rows: tuple[str, ...]
    rejected_rows: tuple[str, ...]
    withheld_cells: tuple[tuple[str, str], ...]
    rejected_cells: tuple[tuple[str, str], ...]
    output_observation_mask: pd.DataFrame
    corrected_cell_status: pd.DataFrame
    provenance_payload: Mapping[str, object]

    @property
    def corrected(self) -> pd.DataFrame:
        """Return the corrected phosphosite matrix."""

        return self.corrected_matrix


@dataclass(frozen=True, slots=True)
class _PreparedMatrix:
    working: pd.DataFrame
    originally_missing: pd.DataFrame
    missing_cells: tuple[tuple[str, str], ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _FactorFit:
    corrected_values: np.ndarray
    estimated_factors: np.ndarray
    singular_values: tuple[float, ...]
    protected_rank: int
    factor_count: int
    adjustment: np.ndarray
    warnings: tuple[str, ...]


class DeterministicSpsRuvStyleExecutor:
    """Estimate unwanted factors from eligible controls and correct a matrix."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        plan: _ResolvedPlanLike,
    ) -> SpsRuvStyleExecutorResult:
        _require_supported_plan(plan)
        _require_matrix(phospho)
        design = _aligned_protected_design(phospho=phospho, plan=plan)
        prepared = _prepare_matrix(phospho=phospho, plan=plan)
        control_positions = _control_positions(phospho=phospho, plan=plan)
        fit = _fit_sps_ruv_style(
            matrix=prepared.working,
            design=design,
            control_positions=control_positions,
            control_weights=_control_weights(plan),
            requested_factors=plan.n_unwanted_factors,
        )

        corrected_matrix = pd.DataFrame(
            fit.corrected_values,
            index=phospho.index.copy(),
            columns=phospho.columns.copy(),
        )
        corrected_matrix = corrected_matrix.mask(prepared.originally_missing, np.nan)
        estimated_factors = pd.DataFrame(
            fit.estimated_factors,
            index=phospho.columns.copy(),
            columns=pd.Index(
                [
                    f"unwanted_factor_{position + 1}"
                    for position in range(fit.factor_count)
                ],
                name="factor",
            ),
        )
        output_observation_mask = ~prepared.originally_missing.copy(deep=True)
        corrected_cell_status = _corrected_cell_status(output_observation_mask)
        withheld_cells = tuple(prepared.missing_cells)
        warnings = (*prepared.warnings, *fit.warnings)
        diagnostics = _diagnostics(
            phospho=phospho,
            corrected=corrected_matrix,
            plan=plan,
            fit=fit,
            originally_missing=prepared.originally_missing,
            withheld_cells=withheld_cells,
            warnings=warnings,
        )
        provenance_payload = _provenance_payload(
            phospho=phospho,
            corrected=corrected_matrix,
            output_observation_mask=output_observation_mask,
            plan=plan,
            diagnostics=diagnostics,
            warnings=warnings,
            estimated_factors=estimated_factors,
        )
        return SpsRuvStyleExecutorResult(
            corrected_matrix=corrected_matrix,
            estimated_unwanted_factors=estimated_factors,
            diagnostics=diagnostics,
            warnings=warnings,
            withheld_rows=(),
            rejected_rows=(),
            withheld_cells=withheld_cells,
            rejected_cells=(),
            output_observation_mask=output_observation_mask,
            corrected_cell_status=corrected_cell_status,
            provenance_payload=provenance_payload,
        )


SpsRuvStyleExecutor = DeterministicSpsRuvStyleExecutor


def _require_supported_plan(plan: _ResolvedPlanLike) -> None:
    method = str(plan.method).strip()
    if method not in SPS_RUV_STYLE_METHODS:
        raise PhosPyInputError(
            "SPS/RUV-style executor requires a resolved SPS/RUV-style plan; "
            f"got method={method!r}"
        )
    if not plan.condition_terms_to_preserve:
        raise PhosPyInputError(
            "SPS/RUV-style executor requires protected condition design terms"
        )
    if not plan.eligible_control_site_rows:
        raise PhosPyInputError(
            "SPS/RUV-style executor requires at least one eligible control row"
        )


def _require_matrix(phospho: pd.DataFrame) -> None:
    if not isinstance(phospho, pd.DataFrame):
        raise PhosPyInputError("SPS/RUV-style executor requires a pandas DataFrame")
    if phospho.shape[0] < 1:
        raise PhosPyInputError(
            "SPS/RUV-style executor requires at least one phosphosite row"
        )
    if phospho.shape[1] < 2:
        raise PhosPyInputError("SPS/RUV-style executor requires at least two samples")
    non_numeric = [
        str(column)
        for column in phospho.columns.tolist()
        if (
            not pd.api.types.is_numeric_dtype(phospho.loc[:, column])
            or pd.api.types.is_bool_dtype(phospho.loc[:, column])
        )
    ]
    if non_numeric:
        raise PhosPyInputError(
            "SPS/RUV-style executor requires numeric phospho columns. "
            "Non-numeric columns: " + ", ".join(non_numeric)
        )


def _aligned_protected_design(
    *,
    phospho: pd.DataFrame,
    plan: _ResolvedPlanLike,
) -> pd.DataFrame:
    design = plan.resolved_design_matrix
    if not isinstance(design, pd.DataFrame):
        raise PhosPyInputError("resolved batch-correction design must be a DataFrame")
    sample_labels = tuple(str(value) for value in phospho.columns.tolist())
    design_index = tuple(str(value) for value in design.index.tolist())
    if design_index != sample_labels:
        raise PhosPyInputError(
            "resolved batch-correction design sample order must match phospho columns"
        )
    missing_terms = [
        term for term in plan.condition_terms_to_preserve if term not in design.columns
    ]
    if missing_terms:
        raise PhosPyInputError(
            "resolved batch-correction design is missing protected terms: "
            + ", ".join(missing_terms)
        )
    protected = design.loc[:, list(plan.condition_terms_to_preserve)].astype("float64")
    values = protected.to_numpy(dtype="float64", copy=True)
    if not np.isfinite(values).all():
        raise PhosPyInputError(
            "resolved batch-correction protected design contains non-finite values"
        )
    if _matrix_rank(values) < int(values.shape[1]):
        raise PhosPyInputError(
            "resolved batch-correction protected design is rank-deficient"
        )
    return protected


def _prepare_matrix(
    *,
    phospho: pd.DataFrame,
    plan: _ResolvedPlanLike,
) -> _PreparedMatrix:
    values = phospho.astype("float64").copy(deep=True)
    feature_ids = tuple(str(value) for value in phospho.index.tolist())
    sample_ids = tuple(str(value) for value in phospho.columns.tolist())
    mask = plan.observation_mask
    if tuple(mask.feature_ids) != feature_ids:
        raise PhosPyInputError("observation mask feature_ids must match phospho.index")
    if tuple(mask.sample_ids) != sample_ids:
        raise PhosPyInputError("observation mask sample_ids must match phospho.columns")

    originally_missing = pd.DataFrame(
        False,
        index=phospho.index.copy(),
        columns=phospho.columns.copy(),
    )
    for feature_id, sample_id in mask.originally_missing_cells:
        originally_missing.loc[feature_id, sample_id] = True

    actual_missing = values.isna()
    if not actual_missing.equals(originally_missing):
        raise PhosPyInputError(
            "SPS/RUV-style executor requires matrix missing cells to match the "
            "resolved observation mask"
        )
    if not bool(actual_missing.to_numpy().any()):
        return _PreparedMatrix(
            working=values,
            originally_missing=originally_missing,
            missing_cells=(),
            warnings=(),
        )

    policy = plan.temporary_imputation_policy
    method = TemporaryImputationMethod.parse(
        policy.method,
        field_name="temporary imputation policy.method",
    )
    if (
        not policy.allowed
        or method is not TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY
    ):
        raise PhosPyInputError(
            "SPS/RUV-style executor currently requires row_median_temporary "
            "temporary imputation for matrices with originally missing cells"
        )
    min_observed_values = _min_observed_values(policy)
    working = _row_median_temporary_impute(
        values,
        missing_mask=actual_missing,
        min_observed_values=min_observed_values,
    )
    missing_cells = tuple(
        (str(feature_id), str(sample_id))
        for feature_id, sample_id in mask.originally_missing_cells
    )
    return _PreparedMatrix(
        working=working,
        originally_missing=originally_missing,
        missing_cells=missing_cells,
        warnings=(
            "originally missing cells were temporarily imputed for numerical "
            "correction and restored to missing in the corrected output",
        ),
    )


def _row_median_temporary_impute(
    matrix: pd.DataFrame,
    *,
    missing_mask: pd.DataFrame,
    min_observed_values: int,
) -> pd.DataFrame:
    imputed = matrix.copy(deep=True)
    for row_id in matrix.index.tolist():
        row = matrix.loc[row_id, :]
        observed = row.dropna()
        if int(observed.shape[0]) < min_observed_values:
            raise PhosPyInputError(
                "SPS/RUV-style executor cannot temporarily impute row "
                f"{str(row_id)!r}; observed values are below min_observed_values"
            )
        row_missing = missing_mask.loc[row_id, :]
        if bool(row_missing.any()):
            imputed.loc[row_id, row_missing] = float(observed.median())
    values = imputed.to_numpy(dtype="float64", copy=True)
    if not np.isfinite(values).all():
        raise PhosPyInputError(
            "SPS/RUV-style executor requires finite values after temporary imputation"
        )
    return imputed


def _min_observed_values(policy: _TemporaryImputationPolicyLike) -> int:
    value = dict(policy.method_parameters).get("min_observed_values", 1)
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(
            "row_median_temporary requires integer min_observed_values"
        )
    return int(value)


def _control_positions(
    *,
    phospho: pd.DataFrame,
    plan: _ResolvedPlanLike,
) -> tuple[int, ...]:
    positions: list[int] = []
    row_count = int(phospho.shape[0])
    for row in plan.eligible_control_site_rows:
        position = int(row.row_position)
        if position < 0 or position >= row_count:
            raise PhosPyInputError(
                "resolved control row_position is outside the phospho matrix"
            )
        expected_site_key = str(row.site_key)
        observed_site_key = str(phospho.index[position])
        if observed_site_key != expected_site_key:
            raise PhosPyInputError(
                "resolved control row_position does not match phospho index "
                f"(expected {expected_site_key!r}, observed {observed_site_key!r})"
            )
        positions.append(position)
    if len(set(positions)) != len(positions):
        raise PhosPyInputError("resolved control rows must not contain duplicates")
    return tuple(positions)


def _control_weights(plan: _ResolvedPlanLike) -> np.ndarray:
    weights: list[float] = []
    for row in plan.eligible_control_site_rows:
        weight = 1.0 if row.weight is None else float(row.weight)
        if not np.isfinite(weight) or weight <= 0.0:
            raise PhosPyInputError(
                "resolved control weights must be positive finite values"
            )
        weights.append(weight)
    return np.asarray(weights, dtype="float64")


def _fit_sps_ruv_style(
    *,
    matrix: pd.DataFrame,
    design: pd.DataFrame,
    control_positions: tuple[int, ...],
    control_weights: np.ndarray,
    requested_factors: int | None,
) -> _FactorFit:
    response = matrix.to_numpy(dtype="float64", copy=True).T
    protected = design.to_numpy(dtype="float64", copy=True)
    protected_rank = _matrix_rank(protected)
    if int(response.shape[0]) <= protected_rank:
        raise PhosPyInputError(
            "SPS/RUV-style executor requires sample degrees of freedom beyond "
            "protected design terms"
        )
    protected_coefficients = np.linalg.pinv(protected) @ response
    protected_fitted = protected @ protected_coefficients
    residual = response - protected_fitted

    control_residual = residual[:, list(control_positions)]
    weighted_control_residual = control_residual * np.sqrt(control_weights)[None, :]
    u_matrix, singular_values, _ = np.linalg.svd(
        weighted_control_residual,
        full_matrices=False,
    )
    numerical_rank = int(
        np.sum(
            singular_values > _svd_tolerance(weighted_control_residual, singular_values)
        )
    )
    requested = 1 if requested_factors is None else int(requested_factors)
    if requested < 1:
        raise PhosPyInputError(
            "SPS/RUV-style executor requires n_unwanted_factors >= 1"
        )
    max_factors = min(
        requested,
        numerical_rank,
        int(response.shape[0]) - protected_rank,
        len(control_positions),
    )
    warnings: tuple[str, ...] = ()
    if max_factors < requested:
        warnings = (
            "requested unwanted factor count exceeded the numerical control "
            "residual rank and was capped",
        )
    if max_factors < 1:
        raise PhosPyInputError(
            "SPS/RUV-style executor could not estimate unwanted factors from "
            "eligible control residuals"
        )

    factors = u_matrix[:, :max_factors] * singular_values[:max_factors][None, :]
    full_design = np.concatenate((protected, factors), axis=1)
    full_rank = _matrix_rank(full_design)
    if full_rank < int(full_design.shape[1]):
        raise PhosPyInputError(
            "SPS/RUV-style executor factor design is rank-deficient after "
            "protecting condition terms"
        )
    coefficients, *_ = np.linalg.lstsq(full_design, response, rcond=None)
    factor_coefficients = coefficients[protected.shape[1] :, :]
    adjustment = factors @ factor_coefficients
    corrected_values = response - adjustment
    return _FactorFit(
        corrected_values=corrected_values.T,
        estimated_factors=factors,
        singular_values=tuple(float(value) for value in singular_values[:max_factors]),
        protected_rank=protected_rank,
        factor_count=max_factors,
        adjustment=adjustment.T,
        warnings=warnings,
    )


def _diagnostics(
    *,
    phospho: pd.DataFrame,
    corrected: pd.DataFrame,
    plan: _ResolvedPlanLike,
    fit: _FactorFit,
    originally_missing: pd.DataFrame,
    withheld_cells: tuple[tuple[str, str], ...],
    warnings: tuple[str, ...],
) -> SpsRuvStyleExecutorDiagnostics:
    adjustment = fit.adjustment
    observed_cell_count = int((~originally_missing).to_numpy().sum())
    return SpsRuvStyleExecutorDiagnostics(
        method=plan.method,
        executor_id=SPS_RUV_STYLE_EXECUTOR_ID,
        status="applied",
        matrix_shape_before=(int(phospho.shape[0]), int(phospho.shape[1])),
        matrix_shape_after=(int(corrected.shape[0]), int(corrected.shape[1])),
        control_site_count=len(plan.eligible_control_site_rows),
        protected_design_terms=tuple(plan.condition_terms_to_preserve),
        protected_design_rank=fit.protected_rank,
        requested_unwanted_factors=plan.n_unwanted_factors,
        estimated_unwanted_factors=fit.factor_count,
        singular_values=fit.singular_values,
        originally_missing_cell_count=int(originally_missing.to_numpy().sum()),
        corrected_observed_cell_count=observed_cell_count,
        withheld_cell_count=len(withheld_cells),
        max_abs_adjustment=float(np.max(np.abs(adjustment))),
        mean_abs_adjustment=float(np.mean(np.abs(adjustment))),
        warnings=warnings,
    )


def _provenance_payload(
    *,
    phospho: pd.DataFrame,
    corrected: pd.DataFrame,
    output_observation_mask: pd.DataFrame,
    plan: _ResolvedPlanLike,
    diagnostics: SpsRuvStyleExecutorDiagnostics,
    warnings: tuple[str, ...],
    estimated_factors: pd.DataFrame,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "executor_id": SPS_RUV_STYLE_EXECUTOR_ID,
        "method": plan.method,
        "algorithm": (
            "control residual SVD factors with protected-design coefficient "
            "preservation"
        ),
        "input_matrix_fingerprint": _fingerprint_payload(
            fingerprint_matrix(phospho, name="batch_correction.sps_ruv.input")
        ),
        "corrected_matrix_fingerprint": _fingerprint_payload(
            fingerprint_matrix(corrected, name="batch_correction.sps_ruv.corrected")
        ),
        "estimated_unwanted_factors_fingerprint": _fingerprint_payload(
            fingerprint_matrix(
                estimated_factors,
                name="batch_correction.sps_ruv.unwanted_factors",
            )
        ),
        "output_observation_mask_fingerprint": _fingerprint_payload(
            fingerprint_matrix(
                output_observation_mask.astype("int8"),
                name="batch_correction.sps_ruv.output_observation_mask",
            )
        ),
        "resolved_plan": plan.to_payload(),
        "provenance_seed_data": dict(plan.provenance_seed_data),
        "diagnostics": diagnostics.to_payload(),
        "warnings": list(warnings),
    }


def _fingerprint_payload(fingerprint: TableFingerprint) -> dict[str, object]:
    return {
        "name": fingerprint.name,
        "rows": int(fingerprint.rows),
        "columns": int(fingerprint.columns),
        "index_name": fingerprint.index_name,
        "column_names": list(fingerprint.column_names),
        "dtypes": list(fingerprint.dtypes),
        "exact_hash_algorithm": fingerprint.exact_hash_algorithm,
        "exact_hash_value": fingerprint.exact_hash_value,
        "tolerance_hash_algorithm": fingerprint.tolerance_hash_algorithm,
        "tolerance_hash_value": fingerprint.tolerance_hash_value,
        "index_structure": fingerprint.index_structure,
        "column_index_structure": fingerprint.column_index_structure,
    }


def _corrected_cell_status(observed_mask: pd.DataFrame) -> pd.DataFrame:
    status = pd.DataFrame(
        "corrected_observed",
        index=observed_mask.index.copy(),
        columns=observed_mask.columns.copy(),
    )
    return status.mask(~observed_mask, "restored_missing")


def _matrix_rank(matrix: np.ndarray) -> int:
    if matrix.size == 0:
        return 0
    return int(np.linalg.matrix_rank(matrix))


def _svd_tolerance(matrix: np.ndarray, singular_values: np.ndarray) -> float:
    if singular_values.size == 0:
        return 0.0
    return float(np.finfo(float).eps * max(matrix.shape) * singular_values[0])


__all__ = [
    "SPS_RUV_STYLE_EXECUTOR_ID",
    "DeterministicSpsRuvStyleExecutor",
    "SpsRuvStyleExecutor",
    "SpsRuvStyleExecutorDiagnostics",
    "SpsRuvStyleExecutorResult",
]
