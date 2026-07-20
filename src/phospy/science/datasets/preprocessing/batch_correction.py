"""Typed batch-correction preprocessing report models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal, TypeAlias

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.provenance.immutability import freeze_json_mapping_with_error_type
from phospy.science.configs.preprocessing.batch_correction import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DatasetBatchCorrectionConfig,
)

BatchCorrectionStatus: TypeAlias = Literal["disabled", "applied", "rejected"]
BatchCorrectionConfoundingCheckStatus: TypeAlias = Literal[
    "not_applicable",
    "not_checked",
    "passed",
    "confounded",
]

BATCH_CORRECTION_STATUS_DISABLED = "disabled"
BATCH_CORRECTION_STATUS_APPLIED = "applied"
BATCH_CORRECTION_STATUS_REJECTED = "rejected"
BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE = "not_applicable"
BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED = "not_checked"
BATCH_CORRECTION_CONFOUNDING_PASSED = "passed"
BATCH_CORRECTION_CONFOUNDING_CONFOUNDED = "confounded"
BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS = (
    "preserve_condition_effects"
)

MatrixShape: TypeAlias = tuple[int, int]


@dataclass(frozen=True, slots=True)
class BatchCorrectionPolicy:
    """Declared batch-correction intent and protected design policy."""

    method: str
    batch_column: str | None
    condition_column: str | None
    design_preservation_policy: str = (
        BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS
    )
    preserve_condition_effects: bool = True
    condition_columns: Sequence[str] | None = ()

    def __post_init__(self) -> None:
        condition_column = _optional_text(self.condition_column)
        condition_columns = _normalize_condition_columns(self.condition_columns)
        if not condition_columns and condition_column is not None:
            condition_columns = (condition_column,)
        if condition_column is None and condition_columns:
            condition_column = condition_columns[0]
        if (
            condition_column is not None
            and condition_columns
            and condition_column != condition_columns[0]
        ):
            raise PhosPyInputError(
                "BatchCorrectionPolicy.condition_column must match the first "
                "condition_columns entry when both are provided"
            )
        object.__setattr__(self, "condition_column", condition_column)
        object.__setattr__(self, "condition_columns", condition_columns)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible policy payload."""

        return {
            "method": self.method,
            "batch_column": self.batch_column,
            "condition_column": self.condition_column,
            "condition_columns": list(self.condition_columns or ()),
            "design_preservation_policy": self.design_preservation_policy,
            "preserve_condition_effects": self.preserve_condition_effects,
        }


@dataclass(frozen=True, slots=True)
class BatchCorrectionDiagnostics:
    """Typed diagnostics for batch-correction provenance."""

    number_of_batches: int | None = None
    batch_levels: tuple[str, ...] = ()
    condition_levels: tuple[str, ...] = ()
    confounding_check_status: BatchCorrectionConfoundingCheckStatus = (
        BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED
    )
    matrix_shape_before: MatrixShape | None = None
    matrix_shape_after: MatrixShape | None = None
    warnings: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible diagnostics payload."""

        return {
            "number_of_batches": self.number_of_batches,
            "batch_levels": list(self.batch_levels),
            "condition_levels": list(self.condition_levels),
            "confounding_check_status": self.confounding_check_status,
            "matrix_shape_before": _shape_to_payload(self.matrix_shape_before),
            "matrix_shape_after": _shape_to_payload(self.matrix_shape_after),
            "warnings": list(self.warnings),
            "limitations": list(self.limitations),
        }


@dataclass(frozen=True, slots=True)
class BatchCorrectionReport:
    """Typed report for batch-correction diagnostics.

    This model records provenance only. It does not validate a design, run
    residualisation, mutate matrices, or imply execution when status is not
    ``"applied"``.
    """

    status: BatchCorrectionStatus
    policy: BatchCorrectionPolicy
    diagnostics: BatchCorrectionDiagnostics = field(
        default_factory=BatchCorrectionDiagnostics
    )

    @property
    def method(self) -> str:
        return self.policy.method

    @property
    def batch_column(self) -> str | None:
        return self.policy.batch_column

    @property
    def condition_column(self) -> str | None:
        return self.policy.condition_column

    @property
    def condition_columns(self) -> tuple[str, ...]:
        return tuple(self.policy.condition_columns or ())

    @property
    def number_of_batches(self) -> int | None:
        return self.diagnostics.number_of_batches

    @property
    def batch_levels(self) -> tuple[str, ...]:
        return self.diagnostics.batch_levels

    @property
    def condition_levels(self) -> tuple[str, ...]:
        return self.diagnostics.condition_levels

    @property
    def design_preservation_policy(self) -> str:
        return self.policy.design_preservation_policy

    @property
    def confounding_check_status(self) -> BatchCorrectionConfoundingCheckStatus:
        return self.diagnostics.confounding_check_status

    @property
    def matrix_shape_before(self) -> MatrixShape | None:
        return self.diagnostics.matrix_shape_before

    @property
    def matrix_shape_after(self) -> MatrixShape | None:
        return self.diagnostics.matrix_shape_after

    @property
    def warnings(self) -> tuple[str, ...]:
        return self.diagnostics.warnings

    @property
    def limitations(self) -> tuple[str, ...]:
        return self.diagnostics.limitations

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible report payload."""

        return {
            "method": self.method,
            "status": self.status,
            "batch_column": self.batch_column,
            "condition_column": self.condition_column,
            "condition_columns": list(self.condition_columns),
            "number_of_batches": self.number_of_batches,
            "batch_levels": list(self.batch_levels),
            "condition_levels": list(self.condition_levels),
            "design_preservation_policy": self.design_preservation_policy,
            "confounding_check_status": self.confounding_check_status,
            "matrix_shape_before": _shape_to_payload(self.matrix_shape_before),
            "matrix_shape_after": _shape_to_payload(self.matrix_shape_after),
            "warnings": list(self.warnings),
            "limitations": list(self.limitations),
            "policy": self.policy.to_payload(),
            "diagnostics": self.diagnostics.to_payload(),
        }


@dataclass(frozen=True, slots=True, init=False)
class BatchCorrectionResult:
    """Matrix correction output with public report and engine diagnostics."""

    report: BatchCorrectionReport
    diagnostics: Mapping[str, object]
    _corrected_matrix: pd.DataFrame = field(init=False, repr=False)

    def __init__(
        self,
        *,
        corrected_matrix: pd.DataFrame,
        report: BatchCorrectionReport,
        diagnostics: Mapping[str, object],
    ) -> None:
        object.__setattr__(
            self,
            "_corrected_matrix",
            own_dataframe(
                corrected_matrix,
                field_name="batch_correction_result.corrected_matrix",
                error_type=PhosPyInputError,
                assume_owned=False,
            ),
        )
        object.__setattr__(self, "report", report)
        object.__setattr__(
            self,
            "diagnostics",
            freeze_json_mapping_with_error_type(
                diagnostics,
                field_name="batch_correction_result.diagnostics",
                error_type=PhosPyInputError,
            ),
        )

    @property
    def corrected_matrix(self) -> pd.DataFrame:
        """Return a corrected phospho-matrix snapshot."""

        return export_dataframe(self._corrected_matrix)

    @property
    def corrected(self) -> pd.DataFrame:
        """Return the corrected phospho matrix."""

        return export_dataframe(self._corrected_matrix)


@dataclass(frozen=True, slots=True)
class _ResidualisationFit:
    corrected_values: np.ndarray
    condition_design_columns: int
    batch_design_columns: int
    full_design_rank: int
    residual_degrees_of_freedom: int
    estimated_batch_contribution: np.ndarray


class LinearResidualizeBatchCorrectionEngine:
    """Apply fixed-effect batch residualisation while preserving condition terms."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        batch_labels: Sequence[object],
        condition_labels: Sequence[object],
        config: DatasetBatchCorrectionConfig,
    ) -> BatchCorrectionResult:
        _require_linear_residualisation_config(config)
        _require_matrix(phospho)
        normalized_batch_labels = _normalize_labels(
            batch_labels,
            expected_length=int(phospho.shape[1]),
            field_name="linear_residualize_batch batch_labels",
        )
        normalized_condition_labels = _normalize_labels(
            condition_labels,
            expected_length=int(phospho.shape[1]),
            field_name="linear_residualize_batch condition_labels",
        )
        batch_levels = _levels_in_order(normalized_batch_labels)
        condition_levels = _levels_in_order(normalized_condition_labels)
        if len(batch_levels) < 2:
            raise PhosPyInputError(
                "linear_residualize_batch requires validation-passed metadata "
                "with at least two batch levels"
            )

        fit = _fit_linear_batch_residualisation(
            phospho=phospho,
            batch_labels=normalized_batch_labels,
            condition_labels=normalized_condition_labels,
        )
        corrected_matrix = pd.DataFrame(
            fit.corrected_values,
            index=phospho.index.copy(),
            columns=phospho.columns.copy(),
        )

        typed_diagnostics = BatchCorrectionDiagnostics(
            number_of_batches=len(batch_levels),
            batch_levels=batch_levels,
            condition_levels=condition_levels,
            confounding_check_status=BATCH_CORRECTION_CONFOUNDING_PASSED,
            matrix_shape_before=_matrix_shape(phospho),
            matrix_shape_after=_matrix_shape(corrected_matrix),
            limitations=("linear residualisation preserves matrix shape",),
        )
        report = BatchCorrectionReport(
            status=BATCH_CORRECTION_STATUS_APPLIED,
            policy=BatchCorrectionPolicy(
                method=config.method,
                batch_column=config.batch_column,
                condition_column=config.condition_column,
                condition_columns=(config.condition_column,),
                preserve_condition_effects=config.preserve_condition_effects,
            ),
            diagnostics=typed_diagnostics,
        )
        diagnostics = _engine_diagnostics(
            phospho=phospho,
            corrected=corrected_matrix,
            batch_levels=batch_levels,
            condition_levels=condition_levels,
            fit=fit,
        )
        return BatchCorrectionResult(
            corrected_matrix=corrected_matrix,
            report=report,
            diagnostics=diagnostics,
        )


BatchCorrectionEngine = LinearResidualizeBatchCorrectionEngine


def _shape_to_payload(shape: MatrixShape | None) -> list[int] | None:
    if shape is None:
        return None
    return [int(shape[0]), int(shape[1])]


def _optional_text(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return None if text == "" else text


def _normalize_condition_columns(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str | bytes | bytearray) or not isinstance(value, Sequence):
        raise PhosPyInputError(
            "BatchCorrectionPolicy.condition_columns must be a sequence of "
            "condition column names"
        )
    columns = tuple(str(column).strip() for column in value)
    if any(column == "" for column in columns):
        raise PhosPyInputError(
            "BatchCorrectionPolicy.condition_columns must contain non-empty "
            "column names"
        )
    if len(set(columns)) != len(columns):
        raise PhosPyInputError(
            "BatchCorrectionPolicy.condition_columns must not contain duplicates"
        )
    return columns


def _require_linear_residualisation_config(
    config: DatasetBatchCorrectionConfig,
) -> None:
    if not isinstance(config, DatasetBatchCorrectionConfig):
        raise PhosPyInputError(
            "linear_residualize_batch requires a validated DatasetBatchCorrectionConfig"
        )
    if config.method != DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH:
        raise PhosPyInputError(
            "linear_residualize_batch engine only executes "
            "method='linear_residualize_batch'"
        )
    if config.preserve_condition_effects is not True:
        raise PhosPyInputError(
            "linear_residualize_batch requires preserve_condition_effects=True"
        )


def _require_matrix(phospho: pd.DataFrame) -> None:
    if not isinstance(phospho, pd.DataFrame):
        raise PhosPyInputError("linear_residualize_batch requires a pandas DataFrame")
    if phospho.shape[0] < 1:
        raise PhosPyInputError(
            "linear_residualize_batch requires at least one phosphosite row"
        )
    if phospho.shape[1] < 2:
        raise PhosPyInputError(
            "linear_residualize_batch requires at least two sample columns"
        )
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
            "linear_residualize_batch requires numeric phospho columns. "
            "Non-numeric columns: " + ", ".join(non_numeric)
        )
    values = phospho.astype("float64").to_numpy(copy=True)
    if not np.isfinite(values).all():
        raise PhosPyInputError(
            "linear_residualize_batch requires finite, non-missing phospho values"
        )


def _normalize_labels(
    labels: Sequence[object],
    *,
    expected_length: int,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(labels, (str, bytes)):
        raise PhosPyInputError(f"{field_name} must be a sequence of labels")

    normalized: list[str] = []
    missing_positions: list[int] = []
    blank_positions: list[int] = []
    for position, value in enumerate(tuple(labels)):
        if _is_missing_value(value):
            missing_positions.append(position)
            continue
        label = str(value).strip()
        if label == "":
            blank_positions.append(position)
            continue
        normalized.append(label)

    if (
        len(normalized) + len(missing_positions) + len(blank_positions)
        != expected_length
    ):
        raise PhosPyInputError(
            f"{field_name} must contain exactly one label per phospho sample "
            f"(expected={expected_length})"
        )
    if missing_positions:
        raise PhosPyInputError(
            f"{field_name} contains missing labels at positions "
            f"{_format_positions(missing_positions)}"
        )
    if blank_positions:
        raise PhosPyInputError(
            f"{field_name} contains blank labels at positions "
            f"{_format_positions(blank_positions)}"
        )
    return tuple(normalized)


def _levels_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    levels: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        levels.append(label)
    return tuple(levels)


def _fit_linear_batch_residualisation(
    *,
    phospho: pd.DataFrame,
    batch_labels: Sequence[str],
    condition_labels: Sequence[str],
) -> _ResidualisationFit:
    condition_design = _treatment_coded_design(
        condition_labels,
        include_intercept=True,
    )
    batch_design = _treatment_coded_design(
        batch_labels,
        include_intercept=False,
    )
    full_design = np.concatenate((condition_design, batch_design), axis=1)
    full_rank = _matrix_rank(full_design)
    if full_rank < int(full_design.shape[1]) or len(batch_labels) <= full_rank:
        raise PhosPyInputError(
            "linear_residualize_batch requires a validation-passed full-rank "
            "condition-plus-batch design with residual degrees of freedom; run "
            "BatchCorrectionAdequacyValidator before correction"
        )

    response = phospho.astype("float64").to_numpy(copy=True).T
    coefficients, *_ = np.linalg.lstsq(full_design, response, rcond=None)
    condition_column_count = int(condition_design.shape[1])
    batch_coefficients = coefficients[condition_column_count:, :]
    estimated_batch_contribution = batch_design @ batch_coefficients
    corrected_values = response - estimated_batch_contribution
    return _ResidualisationFit(
        corrected_values=corrected_values.T,
        condition_design_columns=condition_column_count,
        batch_design_columns=int(batch_design.shape[1]),
        full_design_rank=full_rank,
        residual_degrees_of_freedom=int(len(batch_labels) - full_rank),
        estimated_batch_contribution=estimated_batch_contribution,
    )


def _treatment_coded_design(
    labels: Sequence[str],
    *,
    include_intercept: bool,
) -> np.ndarray:
    levels = _levels_in_order(labels)
    column_count = (1 if include_intercept else 0) + max(len(levels) - 1, 0)
    if column_count == 0:
        return np.empty((len(labels), 0), dtype="float64")

    rows: list[list[float]] = []
    for label in labels:
        row: list[float] = []
        if include_intercept:
            row.append(1.0)
        row.extend(1.0 if label == level else 0.0 for level in levels[1:])
        rows.append(row)
    return np.asarray(rows, dtype="float64")


def _matrix_rank(matrix: np.ndarray) -> int:
    if matrix.size == 0:
        return 0
    return int(np.linalg.matrix_rank(matrix))


def _engine_diagnostics(
    *,
    phospho: pd.DataFrame,
    corrected: pd.DataFrame,
    batch_levels: tuple[str, ...],
    condition_levels: tuple[str, ...],
    fit: _ResidualisationFit,
) -> dict[str, object]:
    estimated_batch_contribution = fit.estimated_batch_contribution
    return {
        "method": DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
        "status": BATCH_CORRECTION_STATUS_APPLIED,
        "number_of_sites": int(phospho.shape[0]),
        "number_of_samples": int(phospho.shape[1]),
        "number_of_batches": len(batch_levels),
        "batch_levels": list(batch_levels),
        "condition_levels": list(condition_levels),
        "condition_design_columns": fit.condition_design_columns,
        "batch_design_columns": fit.batch_design_columns,
        "full_design_rank": fit.full_design_rank,
        "residual_degrees_of_freedom": fit.residual_degrees_of_freedom,
        "matrix_shape_before": list(_matrix_shape(phospho)),
        "matrix_shape_after": list(_matrix_shape(corrected)),
        "max_abs_estimated_batch_contribution": float(
            np.max(np.abs(estimated_batch_contribution))
        ),
        "mean_abs_estimated_batch_contribution": float(
            np.mean(np.abs(estimated_batch_contribution))
        ),
    }


def _matrix_shape(matrix: pd.DataFrame) -> MatrixShape:
    return (int(matrix.shape[0]), int(matrix.shape[1]))


def _is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _format_positions(positions: Sequence[int]) -> str:
    preview = ", ".join(str(position) for position in tuple(positions)[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


__all__ = [
    "BATCH_CORRECTION_CONFOUNDING_CONFOUNDED",
    "BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE",
    "BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED",
    "BATCH_CORRECTION_CONFOUNDING_PASSED",
    "BATCH_CORRECTION_STATUS_APPLIED",
    "BATCH_CORRECTION_STATUS_DISABLED",
    "BATCH_CORRECTION_STATUS_REJECTED",
    "BatchCorrectionEngine",
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "BatchCorrectionResult",
    "LinearResidualizeBatchCorrectionEngine",
]
