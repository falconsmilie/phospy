"""Internal duplicate-correlation scientific contracts.

These models define the planned duplicate-correlation semantics before the
numerical estimator is implemented. They are intentionally not re-exported from
public request/config modules.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import TableFingerprint
from phospy.provenance.serialization.tables import table_fingerprint_to_payload

DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN = (
    "feature_reml_fisher_atanh_trimmed_mean"
)
DUPLICATE_CORRELATION_TRIM_FRACTION = 0.15


class DuplicateCorrelationFeatureStatus(StrEnum):
    """Feature-level status for REML correlation estimation."""

    ESTIMATED = "estimated"
    NOT_ESTIMABLE = "not_estimable"
    NON_CONVERGED = "non_converged"
    NON_FINITE = "non_finite"
    INVALID_COVARIANCE = "invalid_covariance"


class DuplicateCorrelationFailureReason(StrEnum):
    """Workflow-level duplicate-correlation estimator failure reasons."""

    NO_REPEATED_BLOCKS = "no_repeated_blocks"
    MISSING_OR_EMPTY_BLOCK_IDENTITIES = "missing_or_empty_block_identities"
    INSUFFICIENT_OBSERVATIONS_FOR_DESIGN_RANK = (
        "insufficient_observations_for_design_rank"
    )
    RANK_DEFICIENT_FIXED_EFFECTS_DESIGN = "rank_deficient_fixed_effects_design"
    NO_FEATURE_WITH_ESTIMABLE_CORRELATION = "no_feature_with_estimable_correlation"
    NUMERICAL_NON_CONVERGENCE = "numerical_non_convergence"
    INVALID_OR_NON_POSITIVE_DEFINITE_COVARIANCE = (
        "invalid_or_non_positive_definite_covariance"
    )
    UNSUPPORTED_OBSERVATION_WEIGHTING = "unsupported_observation_weighting"
    ALL_ELIGIBLE_FEATURE_ESTIMATES_NON_FINITE = (
        "all_eligible_feature_estimates_non_finite"
    )


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationBlockStructureSummary:
    """Summary of block metadata used as correlation groups."""

    block_id_field_name: str
    sample_count: int
    block_count: int
    repeated_block_count: int
    singleton_block_count: int
    correlated_pair_count: int
    block_levels: tuple[str, ...]

    def __post_init__(self) -> None:
        block_id_field_name = _require_non_empty_string(
            self.block_id_field_name,
            field_name="duplicate_correlation.block_structure.block_id_field_name",
        )
        sample_count = _require_positive_int(
            self.sample_count,
            field_name="duplicate_correlation.block_structure.sample_count",
        )
        block_count = _require_positive_int(
            self.block_count,
            field_name="duplicate_correlation.block_structure.block_count",
        )
        repeated_block_count = _require_non_negative_int(
            self.repeated_block_count,
            field_name="duplicate_correlation.block_structure.repeated_block_count",
        )
        singleton_block_count = _require_non_negative_int(
            self.singleton_block_count,
            field_name="duplicate_correlation.block_structure.singleton_block_count",
        )
        correlated_pair_count = _require_non_negative_int(
            self.correlated_pair_count,
            field_name="duplicate_correlation.block_structure.correlated_pair_count",
        )
        block_levels = _require_non_empty_string_tuple(
            self.block_levels,
            field_name="duplicate_correlation.block_structure.block_levels",
        )
        if len(set(block_levels)) != len(block_levels):
            raise PhosPyInputError(
                "duplicate_correlation.block_structure.block_levels must be unique"
            )
        if block_count != len(block_levels):
            raise PhosPyInputError(
                "duplicate_correlation.block_structure.block_count must match "
                "block_levels length"
            )
        if repeated_block_count + singleton_block_count != block_count:
            raise PhosPyInputError(
                "duplicate_correlation.block_structure repeated plus singleton "
                "blocks must equal block_count"
            )
        if block_count > sample_count:
            raise PhosPyInputError(
                "duplicate_correlation.block_structure.block_count cannot exceed "
                "sample_count"
            )
        minimum_samples = (2 * repeated_block_count) + singleton_block_count
        if sample_count < minimum_samples:
            raise PhosPyInputError(
                "duplicate_correlation.block_structure.sample_count is too small "
                "for the repeated and singleton block counts"
            )
        if repeated_block_count == 0 and correlated_pair_count != 0:
            raise PhosPyInputError(
                "duplicate_correlation.block_structure.correlated_pair_count must "
                "be 0 when no repeated blocks are present"
            )
        if repeated_block_count > 0 and correlated_pair_count < repeated_block_count:
            raise PhosPyInputError(
                "duplicate_correlation.block_structure.correlated_pair_count must "
                "be at least repeated_block_count"
            )
        object.__setattr__(self, "block_id_field_name", block_id_field_name)
        object.__setattr__(self, "sample_count", sample_count)
        object.__setattr__(self, "block_count", block_count)
        object.__setattr__(self, "repeated_block_count", repeated_block_count)
        object.__setattr__(self, "singleton_block_count", singleton_block_count)
        object.__setattr__(self, "correlated_pair_count", correlated_pair_count)
        object.__setattr__(self, "block_levels", block_levels)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible block-structure summary."""

        return {
            "block_id_field_name": self.block_id_field_name,
            "sample_count": self.sample_count,
            "block_count": self.block_count,
            "repeated_block_count": self.repeated_block_count,
            "singleton_block_count": self.singleton_block_count,
            "correlated_pair_count": self.correlated_pair_count,
            "block_levels": list(self.block_levels),
        }


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationFeatureEstimate:
    """Feature-level REML duplicate-correlation estimate."""

    feature_id: str
    status: DuplicateCorrelationFeatureStatus
    correlation: float | None = None
    failure_reason: DuplicateCorrelationFailureReason | None = None

    def __post_init__(self) -> None:
        feature_id = _require_non_empty_string(
            self.feature_id,
            field_name="duplicate_correlation.feature_estimate.feature_id",
        )
        status = _require_feature_status(
            self.status,
            field_name="duplicate_correlation.feature_estimate.status",
        )
        failure_reason = _optional_failure_reason(
            self.failure_reason,
            field_name="duplicate_correlation.feature_estimate.failure_reason",
        )
        if status is DuplicateCorrelationFeatureStatus.ESTIMATED:
            if failure_reason is not None:
                raise PhosPyInputError(
                    "successful duplicate-correlation feature estimates must not "
                    "carry a failure reason"
                )
            correlation = _require_correlation(
                self.correlation,
                field_name="duplicate_correlation.feature_estimate.correlation",
            )
        else:
            if self.correlation is not None:
                raise PhosPyInputError(
                    "failed duplicate-correlation feature estimates must not carry "
                    "a correlation value"
                )
            if failure_reason is None:
                raise PhosPyInputError(
                    "failed duplicate-correlation feature estimates must carry a "
                    "failure reason"
                )
            correlation = None
        object.__setattr__(self, "feature_id", feature_id)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "correlation", correlation)
        object.__setattr__(self, "failure_reason", failure_reason)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible internal feature-estimate payload."""

        return {
            "feature_id": self.feature_id,
            "status": self.status.value,
            "correlation": self.correlation,
            "failure_reason": (
                None if self.failure_reason is None else self.failure_reason.value
            ),
        }


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationConsensusSummary:
    """Feature-free consensus-correlation summary for workflow provenance."""

    method: str
    trim_fraction: float
    success: bool
    consensus_correlation: float | None
    eligible_feature_count: int
    estimated_feature_count: int
    failed_feature_count: int
    non_finite_feature_count: int
    failure_reason: DuplicateCorrelationFailureReason | None = None

    def __post_init__(self) -> None:
        (
            method,
            trim_fraction,
            eligible,
            estimated,
            failed,
            non_finite,
            failure_reason,
            consensus_correlation,
        ) = _validate_consensus_state(
            method=self.method,
            trim_fraction=self.trim_fraction,
            success=self.success,
            consensus_correlation=self.consensus_correlation,
            eligible_feature_count=self.eligible_feature_count,
            estimated_feature_count=self.estimated_feature_count,
            failed_feature_count=self.failed_feature_count,
            non_finite_feature_count=self.non_finite_feature_count,
            failure_reason=self.failure_reason,
            field_prefix="duplicate_correlation.consensus",
        )
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "trim_fraction", trim_fraction)
        object.__setattr__(self, "eligible_feature_count", eligible)
        object.__setattr__(self, "estimated_feature_count", estimated)
        object.__setattr__(self, "failed_feature_count", failed)
        object.__setattr__(self, "non_finite_feature_count", non_finite)
        object.__setattr__(self, "failure_reason", failure_reason)
        object.__setattr__(self, "consensus_correlation", consensus_correlation)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible feature-free consensus payload."""

        return {
            "method": self.method,
            "trim_fraction": self.trim_fraction,
            "success": self.success,
            "consensus_correlation": self.consensus_correlation,
            "eligible_feature_count": self.eligible_feature_count,
            "estimated_feature_count": self.estimated_feature_count,
            "failed_feature_count": self.failed_feature_count,
            "non_finite_feature_count": self.non_finite_feature_count,
            "failure_reason": (
                None if self.failure_reason is None else self.failure_reason.value
            ),
        }


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationConsensusResult:
    """Internal estimator result retaining feature-wise estimates."""

    method: str
    trim_fraction: float
    success: bool
    consensus_correlation: float | None
    eligible_feature_count: int
    estimated_feature_count: int
    failed_feature_count: int
    non_finite_feature_count: int
    failure_reason: DuplicateCorrelationFailureReason | None = None
    feature_estimates: tuple[DuplicateCorrelationFeatureEstimate, ...] = ()

    def __post_init__(self) -> None:
        (
            method,
            trim_fraction,
            eligible,
            estimated,
            failed,
            non_finite,
            failure_reason,
            consensus_correlation,
        ) = _validate_consensus_state(
            method=self.method,
            trim_fraction=self.trim_fraction,
            success=self.success,
            consensus_correlation=self.consensus_correlation,
            eligible_feature_count=self.eligible_feature_count,
            estimated_feature_count=self.estimated_feature_count,
            failed_feature_count=self.failed_feature_count,
            non_finite_feature_count=self.non_finite_feature_count,
            failure_reason=self.failure_reason,
            field_prefix="duplicate_correlation.consensus",
        )
        estimates = _require_feature_estimate_tuple(
            self.feature_estimates,
            field_name="duplicate_correlation.consensus.feature_estimates",
        )
        if estimates:
            _validate_retained_feature_estimate_counts(
                estimates=estimates,
                eligible_feature_count=eligible,
                estimated_feature_count=estimated,
                failed_feature_count=failed,
                non_finite_feature_count=non_finite,
            )
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "trim_fraction", trim_fraction)
        object.__setattr__(self, "eligible_feature_count", eligible)
        object.__setattr__(self, "estimated_feature_count", estimated)
        object.__setattr__(self, "failed_feature_count", failed)
        object.__setattr__(self, "non_finite_feature_count", non_finite)
        object.__setattr__(self, "failure_reason", failure_reason)
        object.__setattr__(self, "consensus_correlation", consensus_correlation)
        object.__setattr__(self, "feature_estimates", estimates)

    def to_summary(self) -> DuplicateCorrelationConsensusSummary:
        """Return the feature-free summary suitable for workflow provenance."""

        return DuplicateCorrelationConsensusSummary(
            method=self.method,
            trim_fraction=self.trim_fraction,
            success=self.success,
            consensus_correlation=self.consensus_correlation,
            eligible_feature_count=self.eligible_feature_count,
            estimated_feature_count=self.estimated_feature_count,
            failed_feature_count=self.failed_feature_count,
            non_finite_feature_count=self.non_finite_feature_count,
            failure_reason=self.failure_reason,
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible internal estimator payload."""

        payload = self.to_summary().to_payload()
        payload["feature_estimates"] = [
            estimate.to_payload() for estimate in self.feature_estimates
        ]
        return payload


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationWorkflowProvenance:
    """Public-workflow-sized duplicate-correlation provenance summary."""

    model: str
    matrix_authority: str
    authoritative_matrix_fingerprint: TableFingerprint
    design_authority: str
    block_authority: str
    estimator_authority: str
    gls_authority: str
    failure_authority: str
    block_structure: DuplicateCorrelationBlockStructureSummary
    consensus: DuplicateCorrelationConsensusSummary
    imputed_values_participated: bool
    imputed_feature_count: int = 0
    imputed_cell_count: int = 0

    def __post_init__(self) -> None:
        for field_name in (
            "model",
            "matrix_authority",
            "design_authority",
            "block_authority",
            "estimator_authority",
            "gls_authority",
            "failure_authority",
        ):
            object.__setattr__(
                self,
                field_name,
                _require_non_empty_string(
                    getattr(self, field_name),
                    field_name=f"duplicate_correlation.workflow.{field_name}",
                ),
            )
        if self.model != "duplicate_correlation":
            raise PhosPyInputError(
                "duplicate_correlation.workflow.model must be 'duplicate_correlation'"
            )
        if not isinstance(
            cast(object, self.authoritative_matrix_fingerprint),
            TableFingerprint,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.authoritative_matrix_fingerprint "
                "must be a TableFingerprint"
            )
        if not isinstance(
            cast(object, self.block_structure),
            DuplicateCorrelationBlockStructureSummary,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.block_structure must be a "
                "DuplicateCorrelationBlockStructureSummary"
            )
        if not isinstance(
            cast(object, self.consensus),
            DuplicateCorrelationConsensusSummary,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.consensus must be a "
                "DuplicateCorrelationConsensusSummary"
            )
        if not isinstance(cast(object, self.imputed_values_participated), bool):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.imputed_values_participated must be "
                "bool"
            )
        imputed_feature_count = _require_non_negative_int(
            self.imputed_feature_count,
            field_name="duplicate_correlation.workflow.imputed_feature_count",
        )
        imputed_cell_count = _require_non_negative_int(
            self.imputed_cell_count,
            field_name="duplicate_correlation.workflow.imputed_cell_count",
        )
        if not self.imputed_values_participated and (
            imputed_feature_count or imputed_cell_count
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow imputed counts require "
                "imputed_values_participated=True"
            )
        object.__setattr__(self, "imputed_feature_count", imputed_feature_count)
        object.__setattr__(self, "imputed_cell_count", imputed_cell_count)

    def to_payload(self) -> dict[str, object]:
        """Return JSON-compatible workflow provenance without feature estimates."""

        return {
            "model": self.model,
            "matrix_authority": self.matrix_authority,
            "authoritative_matrix_fingerprint": table_fingerprint_to_payload(
                self.authoritative_matrix_fingerprint
            ),
            "design_authority": self.design_authority,
            "block_authority": self.block_authority,
            "estimator_authority": self.estimator_authority,
            "gls_authority": self.gls_authority,
            "failure_authority": self.failure_authority,
            "block_structure": self.block_structure.to_payload(),
            "consensus": self.consensus.to_payload(),
            "imputed_values_participated": self.imputed_values_participated,
            "imputed_feature_count": self.imputed_feature_count,
            "imputed_cell_count": self.imputed_cell_count,
        }


def _validate_consensus_state(
    *,
    method: object,
    trim_fraction: object,
    success: object,
    consensus_correlation: object,
    eligible_feature_count: object,
    estimated_feature_count: object,
    failed_feature_count: object,
    non_finite_feature_count: object,
    failure_reason: object,
    field_prefix: str,
) -> tuple[
    str,
    float,
    int,
    int,
    int,
    int,
    DuplicateCorrelationFailureReason | None,
    float | None,
]:
    method_value = _require_canonical_method(
        method,
        field_name=f"{field_prefix}.method",
    )
    trim_fraction_value = _require_fixed_trim_fraction(
        trim_fraction,
        field_name=f"{field_prefix}.trim_fraction",
    )
    if not isinstance(success, bool):
        raise PhosPyInputError(f"{field_prefix}.success must be bool")
    eligible = _require_non_negative_int(
        eligible_feature_count,
        field_name=f"{field_prefix}.eligible_feature_count",
    )
    estimated = _require_non_negative_int(
        estimated_feature_count,
        field_name=f"{field_prefix}.estimated_feature_count",
    )
    failed = _require_non_negative_int(
        failed_feature_count,
        field_name=f"{field_prefix}.failed_feature_count",
    )
    non_finite = _require_non_negative_int(
        non_finite_feature_count,
        field_name=f"{field_prefix}.non_finite_feature_count",
    )
    if estimated + failed != eligible:
        raise PhosPyInputError(
            f"{field_prefix} estimated plus failed feature counts must equal "
            "eligible_feature_count"
        )
    if non_finite > failed:
        raise PhosPyInputError(
            f"{field_prefix}.non_finite_feature_count cannot exceed "
            "failed_feature_count"
        )
    failure_reason_value = _optional_failure_reason(
        failure_reason,
        field_name=f"{field_prefix}.failure_reason",
    )
    if success:
        if failure_reason_value is not None:
            raise PhosPyInputError(
                "successful duplicate-correlation consensus must not carry a "
                "failure reason"
            )
        consensus = _require_correlation(
            consensus_correlation,
            field_name=f"{field_prefix}.consensus_correlation",
        )
        if estimated < 1:
            raise PhosPyInputError(
                "successful duplicate-correlation consensus requires at least "
                "one estimated feature"
            )
    else:
        if consensus_correlation is not None:
            raise PhosPyInputError(
                "failed duplicate-correlation consensus must not carry a "
                "consensus correlation"
            )
        if failure_reason_value is None:
            raise PhosPyInputError(
                "failed duplicate-correlation consensus requires a failure reason"
            )
        consensus = None
    return (
        method_value,
        trim_fraction_value,
        eligible,
        estimated,
        failed,
        non_finite,
        failure_reason_value,
        consensus,
    )


def _validate_retained_feature_estimate_counts(
    *,
    estimates: tuple[DuplicateCorrelationFeatureEstimate, ...],
    eligible_feature_count: int,
    estimated_feature_count: int,
    failed_feature_count: int,
    non_finite_feature_count: int,
) -> None:
    if len(estimates) != eligible_feature_count:
        raise PhosPyInputError(
            "duplicate_correlation.consensus.feature_estimates length must match "
            "eligible_feature_count"
        )
    estimated = sum(
        estimate.status is DuplicateCorrelationFeatureStatus.ESTIMATED
        for estimate in estimates
    )
    failed = len(estimates) - estimated
    non_finite = sum(
        estimate.status is DuplicateCorrelationFeatureStatus.NON_FINITE
        for estimate in estimates
    )
    if estimated != estimated_feature_count:
        raise PhosPyInputError(
            "duplicate_correlation.consensus.estimated_feature_count must match "
            "retained feature estimate statuses"
        )
    if failed != failed_feature_count:
        raise PhosPyInputError(
            "duplicate_correlation.consensus.failed_feature_count must match "
            "retained feature estimate statuses"
        )
    if non_finite != non_finite_feature_count:
        raise PhosPyInputError(
            "duplicate_correlation.consensus.non_finite_feature_count must match "
            "retained feature estimate statuses"
        )


def _require_feature_estimate_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[DuplicateCorrelationFeatureEstimate, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    estimates: list[DuplicateCorrelationFeatureEstimate] = []
    for estimate in cast(Sequence[object], values):
        if not isinstance(estimate, DuplicateCorrelationFeatureEstimate):
            raise PhosPyInputError(
                f"{field_name} must contain DuplicateCorrelationFeatureEstimate values"
            )
        estimates.append(estimate)
    return tuple(estimates)


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value.strip()


def _require_non_empty_string_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    coerced = tuple(
        _require_non_empty_string(value, field_name=f"{field_name}[]")
        for value in cast(Sequence[object], values)
    )
    if not coerced:
        raise PhosPyInputError(f"{field_name} must be non-empty")
    return coerced


def _require_positive_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise PhosPyInputError(f"{field_name} must be a positive integer")
    return int(value)


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise PhosPyInputError(f"{field_name} must be a non-negative integer")
    return int(value)


def _require_canonical_method(value: object, *, field_name: str) -> str:
    method = _require_non_empty_string(value, field_name=field_name)
    if method != DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN:
        raise PhosPyInputError(
            f"{field_name} must be "
            f"{DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN!r}"
        )
    return method


def _require_fixed_trim_fraction(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric")
    fraction = float(value)
    if not math.isfinite(fraction):
        raise PhosPyInputError(f"{field_name} must be finite")
    if not math.isclose(
        fraction,
        DUPLICATE_CORRELATION_TRIM_FRACTION,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise PhosPyInputError(
            f"{field_name} must be the fixed duplicate-correlation trim value "
            f"{DUPLICATE_CORRELATION_TRIM_FRACTION}"
        )
    return fraction


def _require_correlation(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric")
    correlation = float(value)
    if not math.isfinite(correlation) or not -1.0 < correlation < 1.0:
        raise PhosPyInputError(f"{field_name} must be finite and in (-1.0, 1.0)")
    return correlation


def _require_feature_status(
    value: object,
    *,
    field_name: str,
) -> DuplicateCorrelationFeatureStatus:
    try:
        return DuplicateCorrelationFeatureStatus(value)
    except ValueError as error:
        raise PhosPyInputError(
            f"{field_name} must be a supported duplicate-correlation feature status"
        ) from error


def _optional_failure_reason(
    value: object,
    *,
    field_name: str,
) -> DuplicateCorrelationFailureReason | None:
    if value is None:
        return None
    try:
        return DuplicateCorrelationFailureReason(value)
    except ValueError as error:
        raise PhosPyInputError(
            f"{field_name} must be a supported duplicate-correlation failure reason"
        ) from error


__all__ = [
    "DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN",
    "DUPLICATE_CORRELATION_TRIM_FRACTION",
    "DuplicateCorrelationBlockStructureSummary",
    "DuplicateCorrelationConsensusResult",
    "DuplicateCorrelationConsensusSummary",
    "DuplicateCorrelationFailureReason",
    "DuplicateCorrelationFeatureEstimate",
    "DuplicateCorrelationFeatureStatus",
    "DuplicateCorrelationWorkflowProvenance",
]
