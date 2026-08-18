"""Internal duplicate-correlation scientific contracts.

These models define duplicate-correlation estimator, GLS, and workflow
provenance semantics. They are intentionally not re-exported from public
request/config modules.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import TableFingerprint
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)

DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN = (
    "feature_reml_fisher_atanh_trimmed_mean"
)
DUPLICATE_CORRELATION_TRIM_FRACTION = 0.15
DUPLICATE_CORRELATION_MINIMUM_RESIDUAL_DEGREES_OF_FREEDOM = 2.0
DUPLICATE_CORRELATION_WORKFLOW_PROVENANCE_VERSION = "duplicate_correlation_reml_gls_v1"
DUPLICATE_CORRELATION_ESTIMATOR_POLICY_VERSION = (
    "duplicate_correlation_feature_reml_fisher_atanh_trimmed_mean_v1"
)
DUPLICATE_CORRELATION_BLOCK_TREATMENT_CONSENSUS_CORRELATION = "consensus_correlation"
DUPLICATE_CORRELATION_COVARIANCE_STRUCTURE_COMPOUND_SYMMETRY = "compound_symmetry"
DUPLICATE_CORRELATION_ESTIMATOR_FEATURE_WISE_REML = "feature-wise REML"
DUPLICATE_CORRELATION_GLS_FIT_STATUS_FIT = "fit"


class DuplicateCorrelationFeatureStatus(StrEnum):
    """Feature-level status for REML correlation estimation."""

    ESTIMATED = "estimated"
    BOUNDARY_CONVERGED = "boundary_converged"
    INSUFFICIENT_FINITE_OBSERVATIONS = "insufficient_finite_observations"
    LOST_FIXED_EFFECT_ESTIMABILITY = "lost_fixed_effect_estimability"
    INSUFFICIENT_RESIDUAL_DEGREES_OF_FREEDOM = (
        "insufficient_residual_degrees_of_freedom"
    )
    NO_REPEATED_OBSERVATIONS = "no_repeated_observations"
    ZERO_OR_UNUSABLE_RESIDUAL_VARIATION = "zero_or_unusable_residual_variation"
    OPTIMISATION_FAILED = "optimisation_failed"
    NON_FINITE_OBJECTIVE_OR_ESTIMATE = "non_finite_objective_or_estimate"
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
    INSUFFICIENT_FINITE_OBSERVATIONS = "insufficient_finite_observations"
    LOST_FIXED_EFFECT_ESTIMABILITY = "lost_fixed_effect_estimability"
    INSUFFICIENT_RESIDUAL_DEGREES_OF_FREEDOM = (
        "insufficient_residual_degrees_of_freedom"
    )
    NO_REPEATED_OBSERVATIONS_AFTER_SUBSETTING = (
        "no_repeated_observations_after_subsetting"
    )
    ZERO_OR_UNUSABLE_RESIDUAL_VARIATION = "zero_or_unusable_residual_variation"
    OPTIMISATION_FAILED = "optimisation_failed"
    NON_FINITE_OBJECTIVE_OR_ESTIMATE = "non_finite_objective_or_estimate"


class DuplicateCorrelationBoundary(StrEnum):
    """Numerical boundary reached by a feature-level optimisation."""

    LOWER = "lower"
    UPPER = "upper"


_SUCCESS_FEATURE_STATUSES = frozenset(
    {
        DuplicateCorrelationFeatureStatus.ESTIMATED,
        DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED,
    }
)

_NON_FINITE_FEATURE_STATUSES = frozenset(
    {
        DuplicateCorrelationFeatureStatus.NON_FINITE,
        DuplicateCorrelationFeatureStatus.NON_FINITE_OBJECTIVE_OR_ESTIMATE,
    }
)


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationReasonCount:
    """Count of feature-level failures for one typed reason."""

    reason: DuplicateCorrelationFailureReason
    count: int

    def __post_init__(self) -> None:
        reason = _optional_failure_reason(
            self.reason,
            field_name="duplicate_correlation.reason_count.reason",
        )
        if reason is None:
            raise PhosPyInputError(
                "duplicate_correlation.reason_count.reason must not be None"
            )
        object.__setattr__(self, "reason", reason)
        object.__setattr__(
            self,
            "count",
            _require_non_negative_int(
                self.count,
                field_name="duplicate_correlation.reason_count.count",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible reason-count payload."""

        return {"reason": self.reason.value, "count": self.count}

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> DuplicateCorrelationReasonCount:
        """Reconstruct a reason count from a JSON-compatible payload."""

        mapping = _require_mapping_payload(
            payload,
            field_name="duplicate_correlation.reason_count",
        )
        return cls(
            reason=cast(
                DuplicateCorrelationFailureReason,
                mapping.get("reason"),
            ),
            count=cast(int, mapping.get("count")),
        )


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationConvergenceSummary:
    """Feature-level optimisation convergence counts."""

    converged_feature_count: int
    boundary_feature_count: int
    failed_optimisation_feature_count: int
    non_finite_objective_or_estimate_feature_count: int

    def __post_init__(self) -> None:
        for field_name in (
            "converged_feature_count",
            "boundary_feature_count",
            "failed_optimisation_feature_count",
            "non_finite_objective_or_estimate_feature_count",
        ):
            object.__setattr__(
                self,
                field_name,
                _require_non_negative_int(
                    getattr(self, field_name),
                    field_name=f"duplicate_correlation.convergence.{field_name}",
                ),
            )
        if self.boundary_feature_count > self.converged_feature_count:
            raise PhosPyInputError(
                "duplicate_correlation.convergence.boundary_feature_count cannot "
                "exceed converged_feature_count"
            )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible convergence summary."""

        return {
            "converged_feature_count": self.converged_feature_count,
            "boundary_feature_count": self.boundary_feature_count,
            "failed_optimisation_feature_count": (
                self.failed_optimisation_feature_count
            ),
            "non_finite_objective_or_estimate_feature_count": (
                self.non_finite_objective_or_estimate_feature_count
            ),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> DuplicateCorrelationConvergenceSummary:
        """Reconstruct a convergence summary from a JSON-compatible payload."""

        mapping = _require_mapping_payload(
            payload,
            field_name="duplicate_correlation.convergence",
        )
        return cls(
            converged_feature_count=cast(
                int,
                mapping.get("converged_feature_count"),
            ),
            boundary_feature_count=cast(int, mapping.get("boundary_feature_count")),
            failed_optimisation_feature_count=cast(
                int,
                mapping.get("failed_optimisation_feature_count"),
            ),
            non_finite_objective_or_estimate_feature_count=cast(
                int,
                mapping.get("non_finite_objective_or_estimate_feature_count"),
            ),
        )


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationBoundarySummary:
    """Workflow-level positive-definite interval and boundary counts."""

    lower_correlation_bound: float
    upper_correlation_bound: float
    lower_boundary_feature_count: int
    upper_boundary_feature_count: int
    positive_definite_tolerance: float
    fisher_boundary_tolerance: float

    def __post_init__(self) -> None:
        lower = _require_finite_float(
            self.lower_correlation_bound,
            field_name="duplicate_correlation.boundary.lower_correlation_bound",
        )
        upper = _require_finite_float(
            self.upper_correlation_bound,
            field_name="duplicate_correlation.boundary.upper_correlation_bound",
        )
        if not -1.0 < lower < upper < 1.0:
            raise PhosPyInputError(
                "duplicate_correlation.boundary correlation bounds must satisfy "
                "-1.0 < lower < upper < 1.0"
            )
        lower_count = _require_non_negative_int(
            self.lower_boundary_feature_count,
            field_name=("duplicate_correlation.boundary.lower_boundary_feature_count"),
        )
        upper_count = _require_non_negative_int(
            self.upper_boundary_feature_count,
            field_name=("duplicate_correlation.boundary.upper_boundary_feature_count"),
        )
        pd_tolerance = _require_positive_float(
            self.positive_definite_tolerance,
            field_name=("duplicate_correlation.boundary.positive_definite_tolerance"),
        )
        fisher_tolerance = _require_positive_float(
            self.fisher_boundary_tolerance,
            field_name="duplicate_correlation.boundary.fisher_boundary_tolerance",
        )
        object.__setattr__(self, "lower_correlation_bound", lower)
        object.__setattr__(self, "upper_correlation_bound", upper)
        object.__setattr__(self, "lower_boundary_feature_count", lower_count)
        object.__setattr__(self, "upper_boundary_feature_count", upper_count)
        object.__setattr__(self, "positive_definite_tolerance", pd_tolerance)
        object.__setattr__(self, "fisher_boundary_tolerance", fisher_tolerance)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible boundary summary."""

        return {
            "lower_correlation_bound": self.lower_correlation_bound,
            "upper_correlation_bound": self.upper_correlation_bound,
            "lower_boundary_feature_count": self.lower_boundary_feature_count,
            "upper_boundary_feature_count": self.upper_boundary_feature_count,
            "positive_definite_tolerance": self.positive_definite_tolerance,
            "fisher_boundary_tolerance": self.fisher_boundary_tolerance,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> DuplicateCorrelationBoundarySummary:
        """Reconstruct a boundary summary from a JSON-compatible payload."""

        mapping = _require_mapping_payload(
            payload,
            field_name="duplicate_correlation.boundary",
        )
        return cls(
            lower_correlation_bound=cast(
                float,
                mapping.get("lower_correlation_bound"),
            ),
            upper_correlation_bound=cast(
                float,
                mapping.get("upper_correlation_bound"),
            ),
            lower_boundary_feature_count=cast(
                int,
                mapping.get("lower_boundary_feature_count"),
            ),
            upper_boundary_feature_count=cast(
                int,
                mapping.get("upper_boundary_feature_count"),
            ),
            positive_definite_tolerance=cast(
                float,
                mapping.get("positive_definite_tolerance"),
            ),
            fisher_boundary_tolerance=cast(
                float,
                mapping.get("fisher_boundary_tolerance"),
            ),
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
    minimum_block_size: int | None = None
    maximum_block_size: int | None = None

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
        minimum_block_size = _optional_positive_int(
            self.minimum_block_size,
            field_name="duplicate_correlation.block_structure.minimum_block_size",
        )
        maximum_block_size = _optional_positive_int(
            self.maximum_block_size,
            field_name="duplicate_correlation.block_structure.maximum_block_size",
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
        if minimum_block_size is not None and maximum_block_size is not None:
            if minimum_block_size > maximum_block_size:
                raise PhosPyInputError(
                    "duplicate_correlation.block_structure minimum_block_size "
                    "cannot exceed maximum_block_size"
                )
            if maximum_block_size > sample_count:
                raise PhosPyInputError(
                    "duplicate_correlation.block_structure.maximum_block_size "
                    "cannot exceed sample_count"
                )
            if singleton_block_count > 0 and minimum_block_size != 1:
                raise PhosPyInputError(
                    "duplicate_correlation.block_structure.minimum_block_size "
                    "must be 1 when singleton blocks are present"
                )
            if singleton_block_count == 0 and minimum_block_size < 2:
                raise PhosPyInputError(
                    "duplicate_correlation.block_structure.minimum_block_size "
                    "must be at least 2 when no singleton blocks are present"
                )
            if repeated_block_count > 0 and maximum_block_size < 2:
                raise PhosPyInputError(
                    "duplicate_correlation.block_structure.maximum_block_size "
                    "must be at least 2 when repeated blocks are present"
                )
            if repeated_block_count == 0 and maximum_block_size != 1:
                raise PhosPyInputError(
                    "duplicate_correlation.block_structure.maximum_block_size "
                    "must be 1 when no repeated blocks are present"
                )
        object.__setattr__(self, "block_id_field_name", block_id_field_name)
        object.__setattr__(self, "sample_count", sample_count)
        object.__setattr__(self, "block_count", block_count)
        object.__setattr__(self, "repeated_block_count", repeated_block_count)
        object.__setattr__(self, "singleton_block_count", singleton_block_count)
        object.__setattr__(self, "correlated_pair_count", correlated_pair_count)
        object.__setattr__(self, "block_levels", block_levels)
        object.__setattr__(self, "minimum_block_size", minimum_block_size)
        object.__setattr__(self, "maximum_block_size", maximum_block_size)

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
            "minimum_block_size": self.minimum_block_size,
            "maximum_block_size": self.maximum_block_size,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> DuplicateCorrelationBlockStructureSummary:
        """Reconstruct block structure from a JSON-compatible payload."""

        mapping = _require_mapping_payload(
            payload,
            field_name="duplicate_correlation.block_structure",
        )
        return cls(
            block_id_field_name=cast(str, mapping.get("block_id_field_name")),
            sample_count=cast(int, mapping.get("sample_count")),
            block_count=cast(int, mapping.get("block_count")),
            repeated_block_count=cast(int, mapping.get("repeated_block_count")),
            singleton_block_count=cast(int, mapping.get("singleton_block_count")),
            correlated_pair_count=cast(int, mapping.get("correlated_pair_count")),
            block_levels=_string_tuple_from_payload(
                mapping.get("block_levels"),
                field_name="duplicate_correlation.block_structure.block_levels",
            ),
            minimum_block_size=cast(int | None, mapping.get("minimum_block_size")),
            maximum_block_size=cast(int | None, mapping.get("maximum_block_size")),
        )


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationFeatureEstimate:
    """Feature-level REML duplicate-correlation estimate."""

    feature_id: str
    status: DuplicateCorrelationFeatureStatus
    correlation: float | None = None
    failure_reason: DuplicateCorrelationFailureReason | None = None
    observed_value_count: int | None = None
    design_rank: int | None = None
    residual_degrees_of_freedom: float | None = None
    objective_value: float | None = None
    lower_correlation_bound: float | None = None
    upper_correlation_bound: float | None = None
    boundary: DuplicateCorrelationBoundary | None = None

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
        if status in _SUCCESS_FEATURE_STATUSES:
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
        observed_value_count = _optional_non_negative_int(
            self.observed_value_count,
            field_name="duplicate_correlation.feature_estimate.observed_value_count",
        )
        design_rank = _optional_non_negative_int(
            self.design_rank,
            field_name="duplicate_correlation.feature_estimate.design_rank",
        )
        residual_degrees_of_freedom = _optional_non_negative_float(
            self.residual_degrees_of_freedom,
            field_name=(
                "duplicate_correlation.feature_estimate.residual_degrees_of_freedom"
            ),
        )
        objective_value = _optional_finite_float(
            self.objective_value,
            field_name="duplicate_correlation.feature_estimate.objective_value",
        )
        lower_correlation_bound = _optional_correlation_bound(
            self.lower_correlation_bound,
            field_name=(
                "duplicate_correlation.feature_estimate.lower_correlation_bound"
            ),
        )
        upper_correlation_bound = _optional_correlation_bound(
            self.upper_correlation_bound,
            field_name=(
                "duplicate_correlation.feature_estimate.upper_correlation_bound"
            ),
        )
        if (
            lower_correlation_bound is not None
            and upper_correlation_bound is not None
            and lower_correlation_bound >= upper_correlation_bound
        ):
            raise PhosPyInputError(
                "duplicate_correlation.feature_estimate correlation bounds must "
                "satisfy lower < upper"
            )
        boundary = _optional_boundary(
            self.boundary,
            field_name="duplicate_correlation.feature_estimate.boundary",
        )
        if status is DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED:
            if boundary is None:
                raise PhosPyInputError(
                    "boundary-converged duplicate-correlation feature estimates "
                    "must identify the reached boundary"
                )
            if objective_value is None:
                raise PhosPyInputError(
                    "boundary-converged duplicate-correlation feature estimates "
                    "must carry the finite REML objective value"
                )
        if (
            status is DuplicateCorrelationFeatureStatus.ESTIMATED
            and boundary is not None
        ):
            raise PhosPyInputError(
                "interior duplicate-correlation feature estimates must not carry "
                "a boundary marker"
            )
        object.__setattr__(self, "feature_id", feature_id)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "correlation", correlation)
        object.__setattr__(self, "failure_reason", failure_reason)
        object.__setattr__(self, "observed_value_count", observed_value_count)
        object.__setattr__(self, "design_rank", design_rank)
        object.__setattr__(
            self,
            "residual_degrees_of_freedom",
            residual_degrees_of_freedom,
        )
        object.__setattr__(self, "objective_value", objective_value)
        object.__setattr__(
            self,
            "lower_correlation_bound",
            lower_correlation_bound,
        )
        object.__setattr__(
            self,
            "upper_correlation_bound",
            upper_correlation_bound,
        )
        object.__setattr__(self, "boundary", boundary)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible internal feature-estimate payload."""

        return {
            "feature_id": self.feature_id,
            "status": self.status.value,
            "correlation": self.correlation,
            "failure_reason": (
                None if self.failure_reason is None else self.failure_reason.value
            ),
            "observed_value_count": self.observed_value_count,
            "design_rank": self.design_rank,
            "residual_degrees_of_freedom": self.residual_degrees_of_freedom,
            "objective_value": self.objective_value,
            "lower_correlation_bound": self.lower_correlation_bound,
            "upper_correlation_bound": self.upper_correlation_bound,
            "boundary": None if self.boundary is None else self.boundary.value,
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

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> DuplicateCorrelationConsensusSummary:
        """Reconstruct a feature-free consensus summary from payload."""

        mapping = _require_mapping_payload(
            payload,
            field_name="duplicate_correlation.consensus",
        )
        return cls(
            method=cast(str, mapping.get("method")),
            trim_fraction=cast(float, mapping.get("trim_fraction")),
            success=_require_bool(
                mapping.get("success"),
                field_name="duplicate_correlation.consensus.success",
            ),
            consensus_correlation=cast(
                float | None,
                mapping.get("consensus_correlation"),
            ),
            eligible_feature_count=cast(
                int,
                mapping.get("eligible_feature_count"),
            ),
            estimated_feature_count=cast(
                int,
                mapping.get("estimated_feature_count"),
            ),
            failed_feature_count=cast(int, mapping.get("failed_feature_count")),
            non_finite_feature_count=cast(
                int,
                mapping.get("non_finite_feature_count"),
            ),
            failure_reason=cast(
                DuplicateCorrelationFailureReason | None,
                mapping.get("failure_reason"),
            ),
        )


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
    attempted_feature_count: int | None = None
    trimmed_feature_count_each_tail: int = 0
    retained_feature_count_after_trimming: int | None = None
    failure_reason_counts: tuple[DuplicateCorrelationReasonCount, ...] = ()
    convergence_summary: DuplicateCorrelationConvergenceSummary | None = None
    boundary_summary: DuplicateCorrelationBoundarySummary | None = None
    block_structure: DuplicateCorrelationBlockStructureSummary | None = None
    design_rank: int | None = None
    sample_count: int | None = None

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
        attempted = _optional_non_negative_int(
            self.attempted_feature_count,
            field_name="duplicate_correlation.consensus.attempted_feature_count",
        )
        if attempted is None:
            attempted = eligible
        if attempted < eligible:
            raise PhosPyInputError(
                "duplicate_correlation.consensus.attempted_feature_count cannot "
                "be smaller than eligible_feature_count"
            )
        trimmed = _require_non_negative_int(
            self.trimmed_feature_count_each_tail,
            field_name=(
                "duplicate_correlation.consensus.trimmed_feature_count_each_tail"
            ),
        )
        retained = _optional_non_negative_int(
            self.retained_feature_count_after_trimming,
            field_name=(
                "duplicate_correlation.consensus.retained_feature_count_after_trimming"
            ),
        )
        if retained is None:
            retained = max(estimated - (2 * trimmed), 0)
        if retained + (2 * trimmed) > estimated:
            raise PhosPyInputError(
                "duplicate_correlation.consensus trimmed plus retained feature "
                "counts cannot exceed estimated_feature_count"
            )
        failure_reason_counts = _require_reason_count_tuple(
            self.failure_reason_counts,
            field_name="duplicate_correlation.consensus.failure_reason_counts",
        )
        convergence_summary = _optional_convergence_summary(
            self.convergence_summary,
            field_name="duplicate_correlation.consensus.convergence_summary",
        )
        boundary_summary = _optional_boundary_summary(
            self.boundary_summary,
            field_name="duplicate_correlation.consensus.boundary_summary",
        )
        block_structure = _optional_block_structure(
            self.block_structure,
            field_name="duplicate_correlation.consensus.block_structure",
        )
        design_rank = _optional_non_negative_int(
            self.design_rank,
            field_name="duplicate_correlation.consensus.design_rank",
        )
        sample_count = _optional_non_negative_int(
            self.sample_count,
            field_name="duplicate_correlation.consensus.sample_count",
        )
        if estimates:
            _validate_retained_feature_estimate_counts(
                estimates=estimates,
                eligible_feature_count=eligible,
                estimated_feature_count=estimated,
                failed_feature_count=failed,
                non_finite_feature_count=non_finite,
            )
            if failure_reason_counts:
                _validate_failure_reason_counts(
                    estimates=estimates,
                    failure_reason_counts=failure_reason_counts,
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
        object.__setattr__(self, "attempted_feature_count", attempted)
        object.__setattr__(self, "trimmed_feature_count_each_tail", trimmed)
        object.__setattr__(
            self,
            "retained_feature_count_after_trimming",
            retained,
        )
        object.__setattr__(self, "failure_reason_counts", failure_reason_counts)
        object.__setattr__(self, "convergence_summary", convergence_summary)
        object.__setattr__(self, "boundary_summary", boundary_summary)
        object.__setattr__(self, "block_structure", block_structure)
        object.__setattr__(self, "design_rank", design_rank)
        object.__setattr__(self, "sample_count", sample_count)

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
        payload["attempted_feature_count"] = self.attempted_feature_count
        payload["trimmed_feature_count_each_tail"] = (
            self.trimmed_feature_count_each_tail
        )
        payload["retained_feature_count_after_trimming"] = (
            self.retained_feature_count_after_trimming
        )
        payload["failure_reason_counts"] = [
            reason_count.to_payload() for reason_count in self.failure_reason_counts
        ]
        payload["convergence_summary"] = (
            None
            if self.convergence_summary is None
            else self.convergence_summary.to_payload()
        )
        payload["boundary_summary"] = (
            None
            if self.boundary_summary is None
            else self.boundary_summary.to_payload()
        )
        payload["block_structure"] = (
            None if self.block_structure is None else self.block_structure.to_payload()
        )
        payload["design_rank"] = self.design_rank
        payload["sample_count"] = self.sample_count
        payload["feature_estimates"] = [
            estimate.to_payload() for estimate in self.feature_estimates
        ]
        return payload


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationWorkflowProvenance:
    """Public-workflow-sized duplicate-correlation provenance summary."""

    model: str
    provenance_version: str
    requested_paired_design_policy: str
    normalised_paired_design_policy: str
    block_treatment: str
    covariance_structure: str
    estimator: str
    estimator_policy_version: str
    trim_fraction: float
    matrix_authority: str
    analysis_matrix_fingerprint: TableFingerprint
    authoritative_matrix_fingerprint: TableFingerprint
    design_authority: str
    design_fingerprint: TableFingerprint
    block_authority: str
    block_assignment_fingerprint: TableFingerprint
    estimator_authority: str
    gls_authority: str
    failure_authority: str
    block_structure: DuplicateCorrelationBlockStructureSummary
    consensus: DuplicateCorrelationConsensusSummary
    attempted_feature_count: int
    trimmed_feature_count_each_tail: int
    retained_feature_count_after_trimming: int
    failure_reason_counts: tuple[DuplicateCorrelationReasonCount, ...]
    convergence_summary: DuplicateCorrelationConvergenceSummary
    boundary_summary: DuplicateCorrelationBoundarySummary
    sample_count: int
    block_count: int
    repeated_block_count: int
    singleton_block_count: int
    minimum_block_size: int
    maximum_block_size: int
    design_rank: int
    gls_fit_status: str
    imputed_values_participated: bool
    imputed_feature_count: int = 0
    imputed_cell_count: int = 0

    def __post_init__(self) -> None:
        for field_name in (
            "model",
            "provenance_version",
            "requested_paired_design_policy",
            "normalised_paired_design_policy",
            "block_treatment",
            "covariance_structure",
            "estimator",
            "estimator_policy_version",
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
        if self.provenance_version != DUPLICATE_CORRELATION_WORKFLOW_PROVENANCE_VERSION:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.provenance_version must be "
                f"{DUPLICATE_CORRELATION_WORKFLOW_PROVENANCE_VERSION!r}"
            )
        if self.normalised_paired_design_policy != "duplicate_correlation":
            raise PhosPyInputError(
                "duplicate_correlation.workflow.normalised_paired_design_policy "
                "must be 'duplicate_correlation'"
            )
        if self.requested_paired_design_policy != "duplicate_correlation":
            raise PhosPyInputError(
                "duplicate_correlation.workflow.requested_paired_design_policy "
                "must be 'duplicate_correlation'"
            )
        if (
            self.block_treatment
            != DUPLICATE_CORRELATION_BLOCK_TREATMENT_CONSENSUS_CORRELATION
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.block_treatment must be "
                f"{DUPLICATE_CORRELATION_BLOCK_TREATMENT_CONSENSUS_CORRELATION!r}"
            )
        if (
            self.covariance_structure
            != DUPLICATE_CORRELATION_COVARIANCE_STRUCTURE_COMPOUND_SYMMETRY
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.covariance_structure must be "
                f"{DUPLICATE_CORRELATION_COVARIANCE_STRUCTURE_COMPOUND_SYMMETRY!r}"
            )
        if self.estimator != DUPLICATE_CORRELATION_ESTIMATOR_FEATURE_WISE_REML:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.estimator must be "
                f"{DUPLICATE_CORRELATION_ESTIMATOR_FEATURE_WISE_REML!r}"
            )
        if (
            self.estimator_policy_version
            != DUPLICATE_CORRELATION_ESTIMATOR_POLICY_VERSION
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.estimator_policy_version must be "
                f"{DUPLICATE_CORRELATION_ESTIMATOR_POLICY_VERSION!r}"
            )
        trim_fraction = _require_fixed_trim_fraction(
            self.trim_fraction,
            field_name="duplicate_correlation.workflow.trim_fraction",
        )
        if not isinstance(
            cast(object, self.analysis_matrix_fingerprint),
            TableFingerprint,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.analysis_matrix_fingerprint must "
                "be a TableFingerprint"
            )
        if not isinstance(
            cast(object, self.authoritative_matrix_fingerprint),
            TableFingerprint,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.authoritative_matrix_fingerprint "
                "must be a TableFingerprint"
            )
        if self.analysis_matrix_fingerprint != self.authoritative_matrix_fingerprint:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.analysis_matrix_fingerprint must "
                "match authoritative_matrix_fingerprint"
            )
        if not isinstance(cast(object, self.design_fingerprint), TableFingerprint):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.design_fingerprint must be a "
                "TableFingerprint"
            )
        if not isinstance(
            cast(object, self.block_assignment_fingerprint),
            TableFingerprint,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.block_assignment_fingerprint must "
                "be a TableFingerprint"
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
        if not self.consensus.success or self.consensus.consensus_correlation is None:
            raise PhosPyInputError(
                "duplicate_correlation.workflow requires a successful consensus "
                "correlation"
            )
        if not math.isclose(
            trim_fraction,
            self.consensus.trim_fraction,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.trim_fraction must match consensus"
            )
        attempted_feature_count = _require_non_negative_int(
            self.attempted_feature_count,
            field_name="duplicate_correlation.workflow.attempted_feature_count",
        )
        if attempted_feature_count < self.consensus.eligible_feature_count:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.attempted_feature_count cannot be "
                "smaller than consensus.eligible_feature_count"
            )
        trimmed = _require_non_negative_int(
            self.trimmed_feature_count_each_tail,
            field_name=(
                "duplicate_correlation.workflow.trimmed_feature_count_each_tail"
            ),
        )
        retained = _require_non_negative_int(
            self.retained_feature_count_after_trimming,
            field_name=(
                "duplicate_correlation.workflow.retained_feature_count_after_trimming"
            ),
        )
        if retained + (2 * trimmed) > self.consensus.estimated_feature_count:
            raise PhosPyInputError(
                "duplicate_correlation.workflow retained plus trimmed feature "
                "counts cannot exceed estimated_feature_count"
            )
        failure_reason_counts = _require_reason_count_tuple(
            self.failure_reason_counts,
            field_name="duplicate_correlation.workflow.failure_reason_counts",
        )
        if (
            sum(reason_count.count for reason_count in failure_reason_counts)
            != self.consensus.failed_feature_count
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.failure_reason_counts must sum to "
                "consensus.failed_feature_count"
            )
        if not isinstance(
            cast(object, self.convergence_summary),
            DuplicateCorrelationConvergenceSummary,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.convergence_summary must be a "
                "DuplicateCorrelationConvergenceSummary"
            )
        if not isinstance(
            cast(object, self.boundary_summary),
            DuplicateCorrelationBoundarySummary,
        ):
            raise PhosPyInputError(
                "duplicate_correlation.workflow.boundary_summary must be a "
                "DuplicateCorrelationBoundarySummary"
            )
        sample_count = _require_positive_int(
            self.sample_count,
            field_name="duplicate_correlation.workflow.sample_count",
        )
        block_count = _require_positive_int(
            self.block_count,
            field_name="duplicate_correlation.workflow.block_count",
        )
        repeated_block_count = _require_non_negative_int(
            self.repeated_block_count,
            field_name="duplicate_correlation.workflow.repeated_block_count",
        )
        singleton_block_count = _require_non_negative_int(
            self.singleton_block_count,
            field_name="duplicate_correlation.workflow.singleton_block_count",
        )
        minimum_block_size = _require_positive_int(
            self.minimum_block_size,
            field_name="duplicate_correlation.workflow.minimum_block_size",
        )
        maximum_block_size = _require_positive_int(
            self.maximum_block_size,
            field_name="duplicate_correlation.workflow.maximum_block_size",
        )
        design_rank = _require_positive_int(
            self.design_rank,
            field_name="duplicate_correlation.workflow.design_rank",
        )
        if sample_count != self.block_structure.sample_count:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.sample_count must match "
                "block_structure.sample_count"
            )
        if block_count != self.block_structure.block_count:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.block_count must match "
                "block_structure.block_count"
            )
        if repeated_block_count != self.block_structure.repeated_block_count:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.repeated_block_count must match "
                "block_structure.repeated_block_count"
            )
        if singleton_block_count != self.block_structure.singleton_block_count:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.singleton_block_count must match "
                "block_structure.singleton_block_count"
            )
        if minimum_block_size != self.block_structure.minimum_block_size:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.minimum_block_size must match "
                "block_structure.minimum_block_size"
            )
        if maximum_block_size != self.block_structure.maximum_block_size:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.maximum_block_size must match "
                "block_structure.maximum_block_size"
            )
        if design_rank >= sample_count:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.design_rank must leave positive "
                "residual information"
            )
        gls_fit_status = _require_non_empty_string(
            self.gls_fit_status,
            field_name="duplicate_correlation.workflow.gls_fit_status",
        )
        if gls_fit_status != DUPLICATE_CORRELATION_GLS_FIT_STATUS_FIT:
            raise PhosPyInputError(
                "duplicate_correlation.workflow.gls_fit_status must be "
                f"{DUPLICATE_CORRELATION_GLS_FIT_STATUS_FIT!r}"
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
        object.__setattr__(self, "trim_fraction", trim_fraction)
        object.__setattr__(self, "attempted_feature_count", attempted_feature_count)
        object.__setattr__(self, "trimmed_feature_count_each_tail", trimmed)
        object.__setattr__(
            self,
            "retained_feature_count_after_trimming",
            retained,
        )
        object.__setattr__(self, "failure_reason_counts", failure_reason_counts)
        object.__setattr__(self, "sample_count", sample_count)
        object.__setattr__(self, "block_count", block_count)
        object.__setattr__(self, "repeated_block_count", repeated_block_count)
        object.__setattr__(self, "singleton_block_count", singleton_block_count)
        object.__setattr__(self, "minimum_block_size", minimum_block_size)
        object.__setattr__(self, "maximum_block_size", maximum_block_size)
        object.__setattr__(self, "design_rank", design_rank)
        object.__setattr__(self, "gls_fit_status", gls_fit_status)
        object.__setattr__(self, "imputed_feature_count", imputed_feature_count)
        object.__setattr__(self, "imputed_cell_count", imputed_cell_count)

    def to_payload(self) -> dict[str, object]:
        """Return JSON-compatible workflow provenance without feature estimates."""

        return {
            "model": self.model,
            "provenance_version": self.provenance_version,
            "requested_paired_design_policy": self.requested_paired_design_policy,
            "normalised_paired_design_policy": self.normalised_paired_design_policy,
            "block_treatment": self.block_treatment,
            "covariance_structure": self.covariance_structure,
            "estimator": self.estimator,
            "estimator_policy_version": self.estimator_policy_version,
            "trim_fraction": self.trim_fraction,
            "matrix_authority": self.matrix_authority,
            "analysis_matrix_fingerprint": table_fingerprint_to_payload(
                self.analysis_matrix_fingerprint
            ),
            "authoritative_matrix_fingerprint": table_fingerprint_to_payload(
                self.authoritative_matrix_fingerprint
            ),
            "design_authority": self.design_authority,
            "design_fingerprint": table_fingerprint_to_payload(self.design_fingerprint),
            "block_authority": self.block_authority,
            "block_assignment_fingerprint": table_fingerprint_to_payload(
                self.block_assignment_fingerprint
            ),
            "estimator_authority": self.estimator_authority,
            "gls_authority": self.gls_authority,
            "failure_authority": self.failure_authority,
            "block_structure": self.block_structure.to_payload(),
            "consensus": self.consensus.to_payload(),
            "attempted_feature_count": self.attempted_feature_count,
            "trimmed_feature_count_each_tail": (self.trimmed_feature_count_each_tail),
            "retained_feature_count_after_trimming": (
                self.retained_feature_count_after_trimming
            ),
            "failure_reason_counts": [
                reason_count.to_payload() for reason_count in self.failure_reason_counts
            ],
            "convergence_summary": self.convergence_summary.to_payload(),
            "boundary_summary": self.boundary_summary.to_payload(),
            "sample_count": self.sample_count,
            "block_count": self.block_count,
            "repeated_block_count": self.repeated_block_count,
            "singleton_block_count": self.singleton_block_count,
            "minimum_block_size": self.minimum_block_size,
            "maximum_block_size": self.maximum_block_size,
            "design_rank": self.design_rank,
            "gls_fit_status": self.gls_fit_status,
            "imputed_values_participated": self.imputed_values_participated,
            "imputed_feature_count": self.imputed_feature_count,
            "imputed_cell_count": self.imputed_cell_count,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> DuplicateCorrelationWorkflowProvenance:
        """Reconstruct duplicate-correlation workflow provenance from payload."""

        mapping = _require_mapping_payload(
            payload,
            field_name="duplicate_correlation.workflow",
        )
        analysis_matrix_fingerprint = table_fingerprint_from_payload(
            _require_mapping_payload(
                mapping.get("analysis_matrix_fingerprint"),
                field_name=(
                    "duplicate_correlation.workflow.analysis_matrix_fingerprint"
                ),
            )
        )
        return cls(
            model=cast(str, mapping.get("model")),
            provenance_version=cast(str, mapping.get("provenance_version")),
            requested_paired_design_policy=cast(
                str,
                mapping.get("requested_paired_design_policy"),
            ),
            normalised_paired_design_policy=cast(
                str,
                mapping.get("normalised_paired_design_policy"),
            ),
            block_treatment=cast(str, mapping.get("block_treatment")),
            covariance_structure=cast(str, mapping.get("covariance_structure")),
            estimator=cast(str, mapping.get("estimator")),
            estimator_policy_version=cast(
                str,
                mapping.get("estimator_policy_version"),
            ),
            trim_fraction=cast(float, mapping.get("trim_fraction")),
            matrix_authority=cast(str, mapping.get("matrix_authority")),
            analysis_matrix_fingerprint=analysis_matrix_fingerprint,
            authoritative_matrix_fingerprint=table_fingerprint_from_payload(
                _require_mapping_payload(
                    mapping.get("authoritative_matrix_fingerprint"),
                    field_name=(
                        "duplicate_correlation.workflow."
                        "authoritative_matrix_fingerprint"
                    ),
                )
            ),
            design_authority=cast(str, mapping.get("design_authority")),
            design_fingerprint=table_fingerprint_from_payload(
                _require_mapping_payload(
                    mapping.get("design_fingerprint"),
                    field_name="duplicate_correlation.workflow.design_fingerprint",
                )
            ),
            block_authority=cast(str, mapping.get("block_authority")),
            block_assignment_fingerprint=table_fingerprint_from_payload(
                _require_mapping_payload(
                    mapping.get("block_assignment_fingerprint"),
                    field_name=(
                        "duplicate_correlation.workflow.block_assignment_fingerprint"
                    ),
                )
            ),
            estimator_authority=cast(str, mapping.get("estimator_authority")),
            gls_authority=cast(str, mapping.get("gls_authority")),
            failure_authority=cast(str, mapping.get("failure_authority")),
            block_structure=DuplicateCorrelationBlockStructureSummary.from_payload(
                _require_mapping_payload(
                    mapping.get("block_structure"),
                    field_name="duplicate_correlation.workflow.block_structure",
                )
            ),
            consensus=DuplicateCorrelationConsensusSummary.from_payload(
                _require_mapping_payload(
                    mapping.get("consensus"),
                    field_name="duplicate_correlation.workflow.consensus",
                )
            ),
            attempted_feature_count=cast(
                int,
                mapping.get("attempted_feature_count"),
            ),
            trimmed_feature_count_each_tail=cast(
                int,
                mapping.get("trimmed_feature_count_each_tail"),
            ),
            retained_feature_count_after_trimming=cast(
                int,
                mapping.get("retained_feature_count_after_trimming"),
            ),
            failure_reason_counts=_reason_counts_from_payload(
                mapping.get("failure_reason_counts"),
                field_name="duplicate_correlation.workflow.failure_reason_counts",
            ),
            convergence_summary=DuplicateCorrelationConvergenceSummary.from_payload(
                _require_mapping_payload(
                    mapping.get("convergence_summary"),
                    field_name="duplicate_correlation.workflow.convergence_summary",
                )
            ),
            boundary_summary=DuplicateCorrelationBoundarySummary.from_payload(
                _require_mapping_payload(
                    mapping.get("boundary_summary"),
                    field_name="duplicate_correlation.workflow.boundary_summary",
                )
            ),
            sample_count=cast(int, mapping.get("sample_count")),
            block_count=cast(int, mapping.get("block_count")),
            repeated_block_count=cast(int, mapping.get("repeated_block_count")),
            singleton_block_count=cast(int, mapping.get("singleton_block_count")),
            minimum_block_size=cast(int, mapping.get("minimum_block_size")),
            maximum_block_size=cast(int, mapping.get("maximum_block_size")),
            design_rank=cast(int, mapping.get("design_rank")),
            gls_fit_status=cast(str, mapping.get("gls_fit_status")),
            imputed_values_participated=_require_bool(
                mapping.get("imputed_values_participated"),
                field_name=(
                    "duplicate_correlation.workflow.imputed_values_participated"
                ),
            ),
            imputed_feature_count=cast(int, mapping.get("imputed_feature_count", 0)),
            imputed_cell_count=cast(int, mapping.get("imputed_cell_count", 0)),
        )


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
        estimate.status in _SUCCESS_FEATURE_STATUSES for estimate in estimates
    )
    failed = len(estimates) - estimated
    non_finite = sum(
        estimate.status in _NON_FINITE_FEATURE_STATUSES for estimate in estimates
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


def _validate_failure_reason_counts(
    *,
    estimates: tuple[DuplicateCorrelationFeatureEstimate, ...],
    failure_reason_counts: tuple[DuplicateCorrelationReasonCount, ...],
) -> None:
    observed: dict[DuplicateCorrelationFailureReason, int] = {}
    for estimate in estimates:
        if estimate.failure_reason is None:
            continue
        observed[estimate.failure_reason] = observed.get(estimate.failure_reason, 0) + 1
    expected = {item.reason: item.count for item in failure_reason_counts}
    if expected != observed:
        raise PhosPyInputError(
            "duplicate_correlation.consensus.failure_reason_counts must match "
            "retained feature estimate failure reasons"
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


def _require_reason_count_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[DuplicateCorrelationReasonCount, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    reason_counts: list[DuplicateCorrelationReasonCount] = []
    seen_reasons: set[DuplicateCorrelationFailureReason] = set()
    for value in cast(Sequence[object], values):
        if not isinstance(value, DuplicateCorrelationReasonCount):
            raise PhosPyInputError(
                f"{field_name} must contain DuplicateCorrelationReasonCount values"
            )
        if value.reason in seen_reasons:
            raise PhosPyInputError(f"{field_name} cannot contain duplicate reasons")
        seen_reasons.add(value.reason)
        reason_counts.append(value)
    return tuple(reason_counts)


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


def _optional_non_negative_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return _require_non_negative_int(value, field_name=field_name)


def _optional_positive_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return _require_positive_int(value, field_name=field_name)


def _require_finite_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric")
    coerced = float(value)
    if not math.isfinite(coerced):
        raise PhosPyInputError(f"{field_name} must be finite")
    return coerced


def _require_positive_float(value: object, *, field_name: str) -> float:
    coerced = _require_finite_float(value, field_name=field_name)
    if coerced <= 0.0:
        raise PhosPyInputError(f"{field_name} must be > 0.0")
    return coerced


def _optional_non_negative_float(value: object, *, field_name: str) -> float | None:
    if value is None:
        return None
    coerced = _require_finite_float(value, field_name=field_name)
    if coerced < 0.0:
        raise PhosPyInputError(f"{field_name} must be non-negative")
    return coerced


def _optional_finite_float(value: object, *, field_name: str) -> float | None:
    if value is None:
        return None
    return _require_finite_float(value, field_name=field_name)


def _optional_correlation_bound(value: object, *, field_name: str) -> float | None:
    if value is None:
        return None
    bound = _require_finite_float(value, field_name=field_name)
    if not -1.0 < bound < 1.0:
        raise PhosPyInputError(f"{field_name} must be in (-1.0, 1.0)")
    return bound


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


def _optional_boundary(
    value: object,
    *,
    field_name: str,
) -> DuplicateCorrelationBoundary | None:
    if value is None:
        return None
    try:
        return DuplicateCorrelationBoundary(value)
    except ValueError as error:
        raise PhosPyInputError(
            f"{field_name} must be a supported duplicate-correlation boundary"
        ) from error


def _optional_convergence_summary(
    value: object,
    *,
    field_name: str,
) -> DuplicateCorrelationConvergenceSummary | None:
    if value is None:
        return None
    if not isinstance(value, DuplicateCorrelationConvergenceSummary):
        raise PhosPyInputError(
            f"{field_name} must be a DuplicateCorrelationConvergenceSummary"
        )
    return value


def _optional_boundary_summary(
    value: object,
    *,
    field_name: str,
) -> DuplicateCorrelationBoundarySummary | None:
    if value is None:
        return None
    if not isinstance(value, DuplicateCorrelationBoundarySummary):
        raise PhosPyInputError(
            f"{field_name} must be a DuplicateCorrelationBoundarySummary"
        )
    return value


def _optional_block_structure(
    value: object,
    *,
    field_name: str,
) -> DuplicateCorrelationBlockStructureSummary | None:
    if value is None:
        return None
    if not isinstance(value, DuplicateCorrelationBlockStructureSummary):
        raise PhosPyInputError(
            f"{field_name} must be a DuplicateCorrelationBlockStructureSummary"
        )
    return value


def duplicate_correlation_workflow_provenance_from_payload(
    payload: Mapping[str, object],
) -> DuplicateCorrelationWorkflowProvenance:
    """Reconstruct duplicate-correlation workflow provenance from payload."""

    return DuplicateCorrelationWorkflowProvenance.from_payload(payload)


def _require_mapping_payload(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    return cast(Mapping[str, object], value)


def _string_tuple_from_payload(
    value: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    return tuple(
        _require_non_empty_string(item, field_name=f"{field_name}[]")
        for item in cast(Sequence[object], value)
    )


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be bool")
    return value


def _reason_counts_from_payload(
    value: object,
    *,
    field_name: str,
) -> tuple[DuplicateCorrelationReasonCount, ...]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    return tuple(
        DuplicateCorrelationReasonCount.from_payload(
            _require_mapping_payload(item, field_name=f"{field_name}[]")
        )
        for item in cast(Sequence[object], value)
    )


__all__ = [
    "DUPLICATE_CORRELATION_BLOCK_TREATMENT_CONSENSUS_CORRELATION",
    "DUPLICATE_CORRELATION_COVARIANCE_STRUCTURE_COMPOUND_SYMMETRY",
    "DUPLICATE_CORRELATION_ESTIMATOR_FEATURE_WISE_REML",
    "DUPLICATE_CORRELATION_ESTIMATOR_POLICY_VERSION",
    "DUPLICATE_CORRELATION_GLS_FIT_STATUS_FIT",
    "DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN",
    "DUPLICATE_CORRELATION_MINIMUM_RESIDUAL_DEGREES_OF_FREEDOM",
    "DUPLICATE_CORRELATION_TRIM_FRACTION",
    "DUPLICATE_CORRELATION_WORKFLOW_PROVENANCE_VERSION",
    "DuplicateCorrelationBoundary",
    "DuplicateCorrelationBoundarySummary",
    "DuplicateCorrelationBlockStructureSummary",
    "DuplicateCorrelationConvergenceSummary",
    "DuplicateCorrelationConsensusResult",
    "DuplicateCorrelationConsensusSummary",
    "DuplicateCorrelationFailureReason",
    "DuplicateCorrelationFeatureEstimate",
    "DuplicateCorrelationFeatureStatus",
    "DuplicateCorrelationReasonCount",
    "DuplicateCorrelationWorkflowProvenance",
    "duplicate_correlation_workflow_provenance_from_payload",
]
