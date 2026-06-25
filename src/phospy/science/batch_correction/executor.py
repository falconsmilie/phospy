"""Deterministic SPS/RUV-style numerical correction executor."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Protocol, TypeAlias

import numpy as np
import pandas as pd

from phospy.contracts.configs.preprocessing import TemporaryImputationMethod
from phospy.errors.input import PhosPyInputError
from phospy.provenance import fingerprint_matrix
from phospy.provenance.models import TableFingerprint
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)

SPS_RUV_STYLE_EXECUTOR_ID = "deterministic_sps_ruv_style_executor_v1"
SPS_RUV_STYLE_METHODS = frozenset({"sps_ruv_style", "control_site_ruv_style"})
_UNSUPPORTED_RUV_III_STYLE_METHOD_MESSAGE = (
    "ruv_iii_style is not currently supported because replicate-aware RUV-III "
    "numerical semantics are not implemented; use the supported native "
    "SPS/RUV-style method only if applicable and explicitly configured; do not "
    "imply equivalence to RUV-III"
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
    batch_terms: tuple[str, ...]
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
    rejected_control_site_count: int
    rejected_control_sites: tuple[Mapping[str, object], ...]
    rejection_reason_counts: Mapping[str, int]
    design_summary: Mapping[str, object]
    protected_design_terms: tuple[str, ...]
    protected_design_rank: int
    requested_unwanted_factors: int | None
    estimated_unwanted_factors: int
    singular_values: tuple[float, ...]
    batch_associated_variance: Mapping[str, object]
    missingness_imputation_summary: Mapping[str, object]
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
            "eligible_control_site_count": self.control_site_count,
            "rejected_control_site_count": self.rejected_control_site_count,
            "rejected_control_sites": [
                dict(site) for site in self.rejected_control_sites
            ],
            "rejection_reason_counts": dict(self.rejection_reason_counts),
            "design_summary": dict(self.design_summary),
            "protected_design_terms": list(self.protected_design_terms),
            "protected_design_rank": self.protected_design_rank,
            "requested_unwanted_factors": self.requested_unwanted_factors,
            "estimated_unwanted_factors": self.estimated_unwanted_factors,
            "singular_values": list(self.singular_values),
            "batch_associated_variance": dict(self.batch_associated_variance),
            "missingness_imputation_summary": dict(self.missingness_imputation_summary),
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
    corrected_preprocessing_output: CorrectedPreprocessingOutput | None = None

    @property
    def corrected(self) -> pd.DataFrame:
        """Return the corrected phosphosite matrix."""

        return self.corrected_matrix


@dataclass(frozen=True, slots=True)
class _PreparedMatrix:
    working: pd.DataFrame
    originally_missing: pd.DataFrame
    actual_missing: pd.DataFrame
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
        corrected_complete = corrected_matrix.copy(deep=True)
        corrected_matrix = corrected_matrix.mask(prepared.actual_missing, np.nan)
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
            working=prepared.working,
            corrected=corrected_matrix,
            corrected_complete=corrected_complete,
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
        executor_result = SpsRuvStyleExecutorResult(
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
        corrected_preprocessing_output = (
            None
            if bool(corrected_matrix.isna().to_numpy().any())
            else CorrectedPreprocessingOutput.from_sps_ruv_style_result(
                executor_result,
                stage_order=plan.stage_order,
            )
        )
        return SpsRuvStyleExecutorResult(
            corrected_matrix=executor_result.corrected_matrix,
            estimated_unwanted_factors=executor_result.estimated_unwanted_factors,
            diagnostics=executor_result.diagnostics,
            warnings=executor_result.warnings,
            withheld_rows=executor_result.withheld_rows,
            rejected_rows=executor_result.rejected_rows,
            withheld_cells=executor_result.withheld_cells,
            rejected_cells=executor_result.rejected_cells,
            output_observation_mask=executor_result.output_observation_mask,
            corrected_cell_status=executor_result.corrected_cell_status,
            provenance_payload=executor_result.provenance_payload,
            corrected_preprocessing_output=corrected_preprocessing_output,
        )


SpsRuvStyleExecutor = DeterministicSpsRuvStyleExecutor


def _require_supported_plan(plan: _ResolvedPlanLike) -> None:
    method = str(plan.method).strip()
    if method == "ruv_iii_style":
        raise PhosPyInputError(_UNSUPPORTED_RUV_III_STYLE_METHOD_MESSAGE)
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
    untracked_missing = actual_missing & ~originally_missing
    if bool(untracked_missing.to_numpy().any()):
        raise PhosPyInputError(
            "SPS/RUV-style executor found matrix missing cells that are not "
            "tracked by the resolved observation mask"
        )
    missing_cells = tuple(
        (str(feature_id), str(sample_id))
        for feature_id, sample_id in mask.originally_missing_cells
    )
    if not bool(actual_missing.to_numpy().any()):
        return _PreparedMatrix(
            working=values,
            originally_missing=originally_missing,
            actual_missing=actual_missing,
            missing_cells=missing_cells,
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
    return _PreparedMatrix(
        working=working,
        originally_missing=originally_missing,
        actual_missing=actual_missing,
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
    max_estimable_factors = min(
        numerical_rank,
        int(response.shape[0]) - protected_rank,
        len(control_positions),
    )
    if max_estimable_factors < 1:
        raise PhosPyInputError(
            "SPS/RUV-style executor could not estimate unwanted factors from "
            "eligible control residuals"
        )
    if max_estimable_factors < requested:
        raise PhosPyInputError(
            "SPS/RUV-style executor received a non-estimable "
            f"n_unwanted_factors={requested}; only {max_estimable_factors} "
            "factor(s) are supported by the eligible control residual rank and "
            "sample/design degrees of freedom. This should have been rejected "
            "during workflow validation."
        )

    factors = u_matrix[:, :requested] * singular_values[:requested][None, :]
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
        singular_values=tuple(float(value) for value in singular_values[:requested]),
        protected_rank=protected_rank,
        factor_count=requested,
        adjustment=adjustment.T,
        warnings=(),
    )


def _diagnostics(
    *,
    phospho: pd.DataFrame,
    working: pd.DataFrame,
    corrected: pd.DataFrame,
    corrected_complete: pd.DataFrame,
    plan: _ResolvedPlanLike,
    fit: _FactorFit,
    originally_missing: pd.DataFrame,
    withheld_cells: tuple[tuple[str, str], ...],
    warnings: tuple[str, ...],
) -> SpsRuvStyleExecutorDiagnostics:
    adjustment = fit.adjustment
    observed_cell_count = int((~originally_missing).to_numpy().sum())
    rejected_control_sites = _rejected_control_sites(plan)
    full_design = plan.resolved_design_matrix
    return SpsRuvStyleExecutorDiagnostics(
        method=plan.method,
        executor_id=SPS_RUV_STYLE_EXECUTOR_ID,
        status="applied",
        matrix_shape_before=(int(phospho.shape[0]), int(phospho.shape[1])),
        matrix_shape_after=(int(corrected.shape[0]), int(corrected.shape[1])),
        control_site_count=len(plan.eligible_control_site_rows),
        rejected_control_site_count=len(rejected_control_sites),
        rejected_control_sites=rejected_control_sites,
        rejection_reason_counts=_rejection_reason_counts(rejected_control_sites),
        design_summary=_design_summary(
            phospho=phospho,
            plan=plan,
            design=full_design,
        ),
        protected_design_terms=tuple(plan.condition_terms_to_preserve),
        protected_design_rank=fit.protected_rank,
        requested_unwanted_factors=plan.n_unwanted_factors,
        estimated_unwanted_factors=fit.factor_count,
        singular_values=fit.singular_values,
        batch_associated_variance=_batch_associated_variance_summary(
            before=working,
            after=corrected_complete,
            plan=plan,
            design=full_design,
        ),
        missingness_imputation_summary=_missingness_imputation_summary(
            plan=plan,
            originally_missing=originally_missing,
            actual_missing=phospho.isna(),
            withheld_cells=withheld_cells,
        ),
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


def _rejected_control_sites(
    plan: _ResolvedPlanLike,
) -> tuple[Mapping[str, object], ...]:
    mapping = plan.provenance_seed_data.get("control_site_mapping")
    if not isinstance(mapping, Mapping):
        return ()

    rejected: list[Mapping[str, object]] = []
    for scope in ("row_eligibility", "unmapped_annotations"):
        rows = mapping.get(scope)
        if not isinstance(rows, list | tuple):
            continue
        for row in rows:
            if not isinstance(row, Mapping) or not _is_reportable_rejected_control(row):
                continue
            reasons = _control_rejection_reasons(row)
            rejected.append(
                {
                    "site_key": row.get("site_key"),
                    "scope": scope,
                    "control_status": row.get("control_status"),
                    "row_position": row.get("row_position"),
                    "annotation_indices": list(
                        _object_sequence(row.get("annotation_indices"))
                    ),
                    "reasons": list(reasons),
                    "primary_reason": reasons[0],
                    "exclusion_reason": row.get("exclusion_reason"),
                }
            )
    return tuple(rejected)


def _is_reportable_rejected_control(row: Mapping[str, object]) -> bool:
    if row.get("control_status") == "control" and row.get("valid") is True:
        return False
    if _object_sequence(row.get("reasons")):
        return True
    if row.get("control_status") in {"excluded", "invalid", "unknown"}:
        return True
    annotation_count = row.get("annotation_count")
    return isinstance(annotation_count, int) and annotation_count > 0


def _control_rejection_reasons(row: Mapping[str, object]) -> tuple[str, ...]:
    raw_reasons = tuple(str(reason) for reason in _object_sequence(row.get("reasons")))
    if raw_reasons:
        return raw_reasons
    status = str(row.get("control_status") or "").strip()
    if status == "excluded":
        exclusion_reason = row.get("exclusion_reason")
        if exclusion_reason is not None and str(exclusion_reason).strip():
            return (str(exclusion_reason).strip(),)
        return ("excluded_control_site",)
    if status == "non_control":
        return ("not_marked_as_control",)
    if status == "unknown":
        return ("unknown_control_status",)
    if status == "invalid":
        return ("invalid_control_site_annotation",)
    return ("not_eligible_control_site",)


def _rejection_reason_counts(
    rejected_control_sites: tuple[Mapping[str, object], ...],
) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rejected_control_sites:
        reason = str(row.get("primary_reason", "not_eligible_control_site"))
        counts[reason] = counts.get(reason, 0) + 1
    return counts


def _design_summary(
    *,
    phospho: pd.DataFrame,
    plan: _ResolvedPlanLike,
    design: pd.DataFrame,
) -> dict[str, object]:
    batch_by_sample = _string_mapping(plan.provenance_seed_data.get("batch_by_sample"))
    condition_by_sample = _string_mapping(
        plan.provenance_seed_data.get("condition_by_sample")
    )
    batch_levels = _levels_in_order(batch_by_sample.values())
    condition_levels = _levels_in_order(condition_by_sample.values())
    values = design.to_numpy(dtype="float64", copy=True)
    return {
        "sample_count": int(phospho.shape[1]),
        "site_count": int(phospho.shape[0]),
        "design_matrix_shape": [int(design.shape[0]), int(design.shape[1])],
        "design_matrix_rank": _matrix_rank(values),
        "batch_terms": list(plan.batch_terms),
        "condition_terms_to_preserve": list(plan.condition_terms_to_preserve),
        "batch_levels": list(batch_levels),
        "condition_levels": list(condition_levels),
        "number_of_batches": len(batch_levels),
        "number_of_conditions": len(condition_levels),
        "samples_per_batch": _label_counts(batch_by_sample.values()),
        "samples_per_condition": _label_counts(condition_by_sample.values()),
        "batch_condition_sample_counts": _batch_condition_counts(
            batch_by_sample=batch_by_sample,
            condition_by_sample=condition_by_sample,
        ),
    }


def _missingness_imputation_summary(
    *,
    plan: _ResolvedPlanLike,
    originally_missing: pd.DataFrame,
    actual_missing: pd.DataFrame,
    withheld_cells: tuple[tuple[str, str], ...],
) -> dict[str, object]:
    missing_count = int(originally_missing.to_numpy().sum())
    actual_missing_count = int(actual_missing.to_numpy().sum())
    upstream_imputed_count = int(
        (originally_missing & ~actual_missing).to_numpy().sum()
    )
    policy = plan.temporary_imputation_policy.to_payload()
    summary: dict[str, object] = {
        "originally_missing_cell_count": missing_count,
        "withheld_cell_count": len(withheld_cells),
        "restored_missing_cell_count": actual_missing_count,
        "temporary_imputation_applied": actual_missing_count > 0,
        "temporary_imputation_allowed": bool(policy.get("allowed")),
        "temporary_imputation_method": policy.get("method"),
        "temporary_imputation_parameters": dict(
            _mapping_or_empty(policy.get("method_parameters"))
        ),
        "output_policy": (
            "temporarily imputed values are not retained; originally missing cells "
            "are restored to missing in the corrected output"
            if actual_missing_count > 0
            else "no temporary imputation was needed"
        ),
    }
    if upstream_imputed_count > 0:
        summary["upstream_imputed_input_cell_count"] = upstream_imputed_count
        summary["output_policy"] = (
            "upstream-imputed values remain numeric for analysis-ready output, "
            "but their observation mask cells remain false"
        )
    return summary


def _batch_associated_variance_summary(
    *,
    before: pd.DataFrame,
    after: pd.DataFrame,
    plan: _ResolvedPlanLike,
    design: pd.DataFrame,
) -> dict[str, object]:
    batch_terms = tuple(str(term) for term in plan.batch_terms)
    if not batch_terms:
        return {
            "status": "not_available",
            "reason": "no batch design terms were resolved",
        }
    missing_terms = tuple(term for term in batch_terms if term not in design.columns)
    if missing_terms:
        return {
            "status": "not_available",
            "reason": "batch design terms missing from resolved design matrix",
            "missing_batch_terms": list(missing_terms),
        }

    protected = design.loc[:, list(plan.condition_terms_to_preserve)].to_numpy(
        dtype="float64",
        copy=True,
    )
    batch = design.loc[:, list(batch_terms)].to_numpy(dtype="float64", copy=True)
    before_summary = _matrix_batch_variance_summary(
        matrix=before,
        protected=protected,
        batch=batch,
    )
    after_summary = _matrix_batch_variance_summary(
        matrix=after,
        protected=protected,
        batch=batch,
    )
    return {
        "status": "computed",
        "computed_on": (
            "complete numerical correction matrix; originally missing cells use "
            "temporary imputed values for diagnostics only"
        ),
        "batch_terms": list(batch_terms),
        "before": before_summary,
        "after": after_summary,
        "delta_mean_r_squared": _optional_float_delta(
            before_summary.get("mean_r_squared"),
            after_summary.get("mean_r_squared"),
        ),
        "delta_median_r_squared": _optional_float_delta(
            before_summary.get("median_r_squared"),
            after_summary.get("median_r_squared"),
        ),
    }


def _matrix_batch_variance_summary(
    *,
    matrix: pd.DataFrame,
    protected: np.ndarray,
    batch: np.ndarray,
) -> dict[str, object]:
    response = matrix.to_numpy(dtype="float64", copy=True).T
    protected_projection = np.linalg.pinv(protected)
    residual = response - protected @ (protected_projection @ response)
    residualized_batch = batch - protected @ (protected_projection @ batch)
    batch_rank = _matrix_rank(residualized_batch)
    if batch_rank < 1:
        return {
            "status": "not_available",
            "reason": "batch terms have no residual rank after protecting conditions",
            "site_count": int(matrix.shape[0]),
            "finite_site_count": 0,
            "batch_design_rank": batch_rank,
        }

    batch_fit = residualized_batch @ (np.linalg.pinv(residualized_batch) @ residual)
    denominator = np.sum(residual * residual, axis=0)
    numerator = np.sum(batch_fit * batch_fit, axis=0)
    ratios = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan, dtype="float64"),
        where=denominator > _variance_tolerance(denominator),
    )
    finite = ratios[np.isfinite(ratios)]
    if finite.size == 0:
        return {
            "status": "computed",
            "site_count": int(matrix.shape[0]),
            "finite_site_count": 0,
            "batch_design_rank": batch_rank,
            "mean_r_squared": None,
            "median_r_squared": None,
            "max_r_squared": None,
        }
    clipped = np.clip(finite, 0.0, 1.0)
    return {
        "status": "computed",
        "site_count": int(matrix.shape[0]),
        "finite_site_count": int(clipped.size),
        "batch_design_rank": batch_rank,
        "mean_r_squared": float(np.mean(clipped)),
        "median_r_squared": float(np.median(clipped)),
        "max_r_squared": float(np.max(clipped)),
    }


def _variance_tolerance(values: np.ndarray) -> float:
    if values.size == 0:
        return 0.0
    return float(np.finfo(float).eps * max(values.size, 1) * np.nanmax(values))


def _optional_float_delta(before: object, after: object) -> float | None:
    if isinstance(before, int | float) and isinstance(after, int | float):
        return float(after) - float(before)
    return None


def _levels_in_order(labels: Iterable[object]) -> tuple[str, ...]:
    levels: list[str] = []
    seen: set[str] = set()
    for label in tuple(labels):
        text = str(label)
        if text in seen:
            continue
        seen.add(text)
        levels.append(text)
    return tuple(levels)


def _label_counts(labels: Iterable[object]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for label in tuple(labels):
        text = str(label)
        counts[text] = counts.get(text, 0) + 1
    return counts


def _batch_condition_counts(
    *,
    batch_by_sample: Mapping[str, str],
    condition_by_sample: Mapping[str, str],
) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for sample, batch in batch_by_sample.items():
        condition = condition_by_sample.get(sample)
        if condition is None:
            continue
        batch_counts = counts.setdefault(batch, {})
        batch_counts[condition] = batch_counts.get(condition, 0) + 1
    return counts


def _string_mapping(value: object) -> dict[str, str]:
    if not isinstance(value, Mapping):
        return {}
    return {str(key): str(item) for key, item in value.items()}


def _mapping_or_empty(value: object) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        return value
    return {}


def _object_sequence(value: object) -> tuple[object, ...]:
    if isinstance(value, list | tuple):
        return tuple(value)
    return ()


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
