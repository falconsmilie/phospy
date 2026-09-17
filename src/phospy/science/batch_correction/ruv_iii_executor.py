"""Preprocessing executor adapter for replicate-aware RUV-III correction."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Protocol, TypeAlias

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance import fingerprint_matrix
from phospy.provenance.models import TableFingerprint
from phospy.science.batch_correction.ruv_iii import (
    RUV_III_ALGORITHM_ID,
    RuvIIIDiagnostics,
    RuvIIIKernel,
    RuvIIIReplicateStructure,
)
from phospy.science.configs.preprocessing import TemporaryImputationMethod
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)

RUV_III_STYLE_METHOD = "ruv_iii_style"
RUV_III_STYLE_EXECUTOR_ID = "replicate_aware_ruv_iii_style_executor_v1"
RUV_III_STYLE_ALGORITHM_DESCRIPTION = (
    "replicate-aware RUV-III estimates unwanted variation from within-replicate-set "
    "residuals and governed negative-control sites"
)
RUV_III_STYLE_PROTECTED_TERM_ROLE = (
    "protected biological structure is expressed by replicate sets; condition "
    "metadata is validated and recorded but is not injected as the fixed-effect "
    "protected design used by sps_ruv_style"
)
RUV_III_STYLE_BATCH_TERM_ROLE = (
    "batch metadata is validated and recorded for diagnostics; replicate-aware "
    "RUV-III estimates unwanted factors from replicate structure rather than "
    "residualizing batch labels as fixed effects"
)

JsonScalar: TypeAlias = str | int | float | bool | None


class _ControlRowLike(Protocol):
    site_key: str
    row_position: int


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


class _ReplicateStructureLike(Protocol):
    replicate_by_sample: Mapping[str, str] | None

    def to_payload(self) -> dict[str, object]: ...


class _ResolvedPlanLike(Protocol):
    method: str
    condition_terms_to_preserve: tuple[str, ...]
    batch_terms: tuple[str, ...]
    replicate_structure: _ReplicateStructureLike
    eligible_control_site_rows: tuple[_ControlRowLike, ...]
    observation_mask: _ObservationMaskLike
    temporary_imputation_policy: _TemporaryImputationPolicyLike
    n_unwanted_factors: int | None
    stage_order: tuple[str, ...]
    provenance_seed_data: Mapping[str, object]

    def to_payload(self) -> dict[str, object]: ...


@dataclass(frozen=True, slots=True)
class RuvIIIStyleExecutorDiagnostics:
    """Workflow-facing diagnostics for one RUV-III correction."""

    method: str
    executor_id: str
    algorithm_description: str
    term_roles: Mapping[str, str]
    status: str
    matrix_shape_before: tuple[int, int]
    matrix_shape_after: tuple[int, int]
    control_site_count: int
    protected_design_terms: tuple[str, ...]
    requested_unwanted_factors: int
    estimated_unwanted_factors: int
    replicate_definition: Mapping[str, object]
    kernel: RuvIIIDiagnostics
    missingness_imputation_summary: Mapping[str, object]
    max_abs_adjustment: float
    mean_abs_adjustment: float
    warnings: tuple[str, ...]

    def to_payload(self) -> dict[str, object]:
        return {
            "method": self.method,
            "executor_id": self.executor_id,
            "algorithm_id": RUV_III_ALGORITHM_ID,
            "algorithm_description": self.algorithm_description,
            "term_roles": dict(self.term_roles),
            "status": self.status,
            "matrix_shape_before": list(self.matrix_shape_before),
            "matrix_shape_after": list(self.matrix_shape_after),
            "control_site_count": self.control_site_count,
            "eligible_control_site_count": self.control_site_count,
            "protected_design_terms": list(self.protected_design_terms),
            "requested_unwanted_factors": self.requested_unwanted_factors,
            "estimated_unwanted_factors": self.estimated_unwanted_factors,
            "replicate_definition": dict(self.replicate_definition),
            "kernel": self.kernel.to_payload(),
            "missingness_imputation_summary": dict(self.missingness_imputation_summary),
            "max_abs_adjustment": self.max_abs_adjustment,
            "mean_abs_adjustment": self.mean_abs_adjustment,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True, slots=True)
class RuvIIIStyleExecutorResult:
    """RUV-III result adapted to the common correction executor contract."""

    corrected_matrix: pd.DataFrame
    estimated_unwanted_factors: pd.DataFrame
    diagnostics: RuvIIIStyleExecutorDiagnostics
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
        return self.corrected_matrix


@dataclass(frozen=True, slots=True)
class _PreparedMatrix:
    working: pd.DataFrame
    originally_missing: pd.DataFrame
    actual_missing: pd.DataFrame
    missing_cells: tuple[tuple[str, str], ...]
    temporary_completion_applied: bool


class RuvIIIStyleExecutor:
    """Apply the standalone RUV-III kernel through preprocessing contracts."""

    def __init__(self, *, kernel: RuvIIIKernel | None = None) -> None:
        self._kernel = kernel or RuvIIIKernel()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        plan: _ResolvedPlanLike,
    ) -> RuvIIIStyleExecutorResult:
        if str(plan.method).strip() != RUV_III_STYLE_METHOD:
            raise PhosPyInputError(
                "RUV-III style executor requires method='ruv_iii_style'"
            )
        if plan.n_unwanted_factors is None:
            raise PhosPyInputError("RUV-III style executor requires an explicit k")
        replicate_by_sample = plan.replicate_structure.replicate_by_sample
        if replicate_by_sample is None:
            raise PhosPyInputError(
                "RUV-III style executor requires an explicit replicate definition"
            )
        replicate_structure = RuvIIIReplicateStructure.from_assignments(
            sample_order=tuple(str(value) for value in phospho.columns.tolist()),
            replicate_by_sample=replicate_by_sample,
        )
        prepared = _prepare_matrix(phospho=phospho, plan=plan)
        controls = tuple(row.site_key for row in plan.eligible_control_site_rows)
        kernel_result = self._kernel.run(
            phospho=prepared.working,
            control_site_keys=controls,
            replicate_structure=replicate_structure,
            k=plan.n_unwanted_factors,
        )
        corrected_complete = kernel_result.corrected_matrix
        corrected = corrected_complete.mask(prepared.actual_missing, np.nan)
        output_observation_mask = ~prepared.originally_missing.copy(deep=True)
        corrected_cell_status = pd.DataFrame(
            "corrected_observed",
            index=phospho.index.copy(),
            columns=phospho.columns.copy(),
        ).mask(~output_observation_mask, "restored_missing")
        adjustment = prepared.working.to_numpy(dtype="float64", copy=True) - (
            corrected_complete.to_numpy(dtype="float64", copy=True)
        )
        warnings: tuple[str, ...] = ()
        missingness = _missingness_summary(prepared=prepared, plan=plan)
        diagnostics = RuvIIIStyleExecutorDiagnostics(
            method=RUV_III_STYLE_METHOD,
            executor_id=RUV_III_STYLE_EXECUTOR_ID,
            algorithm_description=RUV_III_STYLE_ALGORITHM_DESCRIPTION,
            term_roles={
                "protected_biological_structure": RUV_III_STYLE_PROTECTED_TERM_ROLE,
                "batch_terms": RUV_III_STYLE_BATCH_TERM_ROLE,
                "replicate_metadata": (
                    "replicate assignments are active numerical input to RUV-III"
                ),
            },
            status="applied",
            matrix_shape_before=(int(phospho.shape[0]), int(phospho.shape[1])),
            matrix_shape_after=(int(corrected.shape[0]), int(corrected.shape[1])),
            control_site_count=len(controls),
            protected_design_terms=tuple(plan.condition_terms_to_preserve),
            requested_unwanted_factors=int(plan.n_unwanted_factors),
            estimated_unwanted_factors=kernel_result.diagnostics.effective_k,
            replicate_definition=plan.replicate_structure.to_payload(),
            kernel=kernel_result.diagnostics,
            missingness_imputation_summary=missingness,
            max_abs_adjustment=float(np.max(np.abs(adjustment))),
            mean_abs_adjustment=float(np.mean(np.abs(adjustment))),
            warnings=warnings,
        )
        provenance_payload = _provenance_payload(
            phospho=phospho,
            corrected=corrected,
            output_observation_mask=output_observation_mask,
            estimated_factors=kernel_result.estimated_unwanted_factors,
            plan=plan,
            diagnostics=diagnostics,
        )
        result = RuvIIIStyleExecutorResult(
            corrected_matrix=corrected,
            estimated_unwanted_factors=kernel_result.estimated_unwanted_factors,
            diagnostics=diagnostics,
            warnings=warnings,
            withheld_rows=(),
            rejected_rows=(),
            withheld_cells=prepared.missing_cells,
            rejected_cells=(),
            output_observation_mask=output_observation_mask,
            corrected_cell_status=corrected_cell_status,
            provenance_payload=provenance_payload,
        )
        corrected_output = (
            None
            if bool(corrected.isna().to_numpy().any())
            else CorrectedPreprocessingOutput.from_sps_ruv_style_result(
                result,
                stage_order=plan.stage_order,
            )
        )
        return RuvIIIStyleExecutorResult(
            corrected_matrix=result.corrected_matrix,
            estimated_unwanted_factors=result.estimated_unwanted_factors,
            diagnostics=result.diagnostics,
            warnings=result.warnings,
            withheld_rows=result.withheld_rows,
            rejected_rows=result.rejected_rows,
            withheld_cells=result.withheld_cells,
            rejected_cells=result.rejected_cells,
            output_observation_mask=result.output_observation_mask,
            corrected_cell_status=result.corrected_cell_status,
            provenance_payload=result.provenance_payload,
            corrected_preprocessing_output=corrected_output,
        )


def _prepare_matrix(
    *, phospho: pd.DataFrame, plan: _ResolvedPlanLike
) -> _PreparedMatrix:
    values = phospho.astype("float64").copy(deep=True)
    feature_ids = tuple(str(value) for value in phospho.index.tolist())
    sample_ids = tuple(str(value) for value in phospho.columns.tolist())
    mask = plan.observation_mask
    if tuple(mask.feature_ids) != feature_ids or tuple(mask.sample_ids) != sample_ids:
        raise PhosPyInputError(
            "RUV-III observation mask axes must match the correction matrix"
        )
    originally_missing = pd.DataFrame(
        False, index=phospho.index.copy(), columns=phospho.columns.copy()
    )
    for feature_id, sample_id in mask.originally_missing_cells:
        originally_missing.loc[feature_id, sample_id] = True
    actual_missing = values.isna()
    if bool((actual_missing & ~originally_missing).to_numpy().any()):
        raise PhosPyInputError(
            "RUV-III found missing cells not governed by the observation mask"
        )
    finite_or_missing = np.isfinite(values.to_numpy(dtype="float64", copy=True)) | (
        actual_missing.to_numpy(dtype=bool, copy=True)
    )
    if not bool(finite_or_missing.all()):
        raise PhosPyInputError("RUV-III observed input values must be finite")
    completion_applied = bool(actual_missing.to_numpy().any())
    if completion_applied:
        policy = plan.temporary_imputation_policy
        method = TemporaryImputationMethod.parse(
            policy.method,
            field_name="RUV-III temporary imputation policy.method",
        )
        if (
            not policy.allowed
            or method is not TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY
        ):
            raise PhosPyInputError(
                "RUV-III missing values require explicit row_median_temporary "
                "completion with observation-mask restoration"
            )
        parameters = dict(policy.method_parameters)
        minimum = parameters.get("min_observed_values", 1)
        if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
            raise PhosPyInputError(
                "RUV-III row_median_temporary min_observed_values must be a "
                "positive integer"
            )
        for row_id in values.index.tolist():
            row = values.loc[row_id, :]
            observed = row.dropna()
            if int(observed.shape[0]) < minimum:
                raise PhosPyInputError(
                    f"RUV-III cannot temporarily complete row {str(row_id)!r}; "
                    "too few observed values"
                )
            missing = row.isna()
            if bool(missing.any()):
                values.loc[row_id, missing] = float(observed.median())
    return _PreparedMatrix(
        working=values,
        originally_missing=originally_missing,
        actual_missing=actual_missing,
        missing_cells=tuple(mask.originally_missing_cells),
        temporary_completion_applied=completion_applied,
    )


def _missingness_summary(
    *, prepared: _PreparedMatrix, plan: _ResolvedPlanLike
) -> dict[str, object]:
    policy = plan.temporary_imputation_policy.to_payload()
    actual_count = int(prepared.actual_missing.to_numpy(dtype=bool).sum())
    governed_count = int(prepared.originally_missing.to_numpy(dtype=bool).sum())
    upstream_imputed_count = int(
        (prepared.originally_missing & ~prepared.actual_missing).to_numpy().sum()
    )
    raw_parameters = policy.get("method_parameters")
    parameters = dict(raw_parameters) if isinstance(raw_parameters, Mapping) else {}
    return {
        "strategy": (
            "internal_row_median_temporary_then_restore_governed_missing_positions"
            if prepared.temporary_completion_applied
            else "no_internal_completion_required"
        ),
        "temporary_completion_applied": prepared.temporary_completion_applied,
        "temporary_imputation_method": policy.get("method"),
        "temporary_imputation_parameters": parameters,
        "governed_missing_cell_count": governed_count,
        "restored_missing_cell_count": actual_count,
        "upstream_imputed_input_cell_count": upstream_imputed_count,
        "imputed_values_are_observed_evidence": False,
        "output_policy": (
            "actual missing positions are restored after numerical correction"
            if actual_count
            else "upstream-imputed values remain numeric while their observation "
            "mask remains false"
            if upstream_imputed_count
            else "no governed missing positions"
        ),
    }


def _provenance_payload(
    *,
    phospho: pd.DataFrame,
    corrected: pd.DataFrame,
    output_observation_mask: pd.DataFrame,
    estimated_factors: pd.DataFrame,
    plan: _ResolvedPlanLike,
    diagnostics: RuvIIIStyleExecutorDiagnostics,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "executor_id": RUV_III_STYLE_EXECUTOR_ID,
        "method": RUV_III_STYLE_METHOD,
        "algorithm_id": RUV_III_ALGORITHM_ID,
        "algorithm_description": RUV_III_STYLE_ALGORITHM_DESCRIPTION,
        "resolved_plan": dict(plan.to_payload()),
        "provenance_seed_data": dict(plan.provenance_seed_data),
        "replicate_definition": plan.replicate_structure.to_payload(),
        "k": plan.n_unwanted_factors,
        "missing_data_strategy": dict(diagnostics.missingness_imputation_summary),
        "input_matrix_fingerprint": _fingerprint_payload(
            fingerprint_matrix(phospho, name="batch_correction.ruv_iii.input")
        ),
        "corrected_matrix_fingerprint": _fingerprint_payload(
            fingerprint_matrix(corrected, name="batch_correction.ruv_iii.corrected")
        ),
        "estimated_unwanted_factors_fingerprint": _fingerprint_payload(
            fingerprint_matrix(
                estimated_factors,
                name="batch_correction.ruv_iii.unwanted_factors",
            )
        ),
        "output_observation_mask_fingerprint": _fingerprint_payload(
            fingerprint_matrix(
                output_observation_mask.astype("int8"),
                name="batch_correction.ruv_iii.output_observation_mask",
            )
        ),
        "diagnostics": diagnostics.to_payload(),
        "warnings": [],
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


__all__ = [
    "RUV_III_STYLE_ALGORITHM_DESCRIPTION",
    "RUV_III_STYLE_EXECUTOR_ID",
    "RUV_III_STYLE_METHOD",
    "RuvIIIStyleExecutor",
    "RuvIIIStyleExecutorDiagnostics",
    "RuvIIIStyleExecutorResult",
]
