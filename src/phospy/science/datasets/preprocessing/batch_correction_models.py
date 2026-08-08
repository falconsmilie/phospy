"""Batch-correction policy, diagnostics, report, and result models."""

from __future__ import annotations

__phospy_contracts_facade_role__ = "science_owned_public_model"

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal, TypeAlias

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.provenance.immutability import freeze_json_mapping_with_error_type

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


__all__ = [
    "BATCH_CORRECTION_CONFOUNDING_CONFOUNDED",
    "BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE",
    "BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED",
    "BATCH_CORRECTION_CONFOUNDING_PASSED",
    "BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS",
    "BATCH_CORRECTION_STATUS_APPLIED",
    "BATCH_CORRECTION_STATUS_DISABLED",
    "BATCH_CORRECTION_STATUS_REJECTED",
    "BatchCorrectionConfoundingCheckStatus",
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "BatchCorrectionResult",
    "BatchCorrectionStatus",
    "MatrixShape",
]
