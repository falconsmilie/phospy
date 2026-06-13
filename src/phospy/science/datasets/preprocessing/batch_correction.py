"""Typed batch-correction preprocessing report models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, TypeAlias

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

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible policy payload."""

        return {
            "method": self.method,
            "batch_column": self.batch_column,
            "condition_column": self.condition_column,
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


def _shape_to_payload(shape: MatrixShape | None) -> list[int] | None:
    if shape is None:
        return None
    return [int(shape[0]), int(shape[1])]


__all__ = [
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
]
