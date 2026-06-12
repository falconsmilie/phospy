"""Public differential workflow configuration models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, cast

from phospy.errors.validation import WorkflowValidationError
from phospy.science.differential.models import EmpiricalBayesConfig
from phospy.science.differential.policy_models import TechnicalReplicatePolicy

PairedDesignPolicy = Literal["reject", "fixed_block"]
MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG = "benjamini_hochberg"
SUPPORTED_MULTIPLE_TESTING_METHODS: tuple[str, ...] = (
    MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG,
)
PAIRED_DESIGN_POLICY_REJECT: PairedDesignPolicy = "reject"
PAIRED_DESIGN_POLICY_FIXED_BLOCK: PairedDesignPolicy = "fixed_block"
SUPPORTED_PAIRED_DESIGN_POLICIES: tuple[PairedDesignPolicy, ...] = (
    PAIRED_DESIGN_POLICY_REJECT,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
)


@dataclass(frozen=True, slots=True)
class MultipleTestingConfig:
    """Public multiple-testing policy for differential analysis."""

    method: str = MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG

    def __post_init__(self) -> None:
        if self.method not in SUPPORTED_MULTIPLE_TESTING_METHODS:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_MULTIPLE_TESTING_METHODS
            )
            raise WorkflowValidationError(
                f"differential.multiple_testing.method must be one of: {supported}"
            )


@dataclass(frozen=True, slots=True)
class DifferentialAnalysisConfig:
    """Public configuration for differential analysis.

    Fixed-effect covariate declarations live on ``ExperimentalDesign``.
    ``paired_design_policy`` records user intent for explicit paired or blocked
    designs only. It does not infer ``block_id`` values and does not enable
    mixed-effects modelling. ``"fixed_block"`` requires complete block metadata
    and validates a fixed-effect block design matrix before execution.
    """

    technical_replicate_policy: TechnicalReplicatePolicy = (
        TechnicalReplicatePolicy.REJECT
    )
    paired_design_policy: PairedDesignPolicy = PAIRED_DESIGN_POLICY_REJECT
    allow_design_subset: bool = False
    minimum_condition_replicates: int = 2
    empirical_bayes: EmpiricalBayesConfig = field(default_factory=EmpiricalBayesConfig)
    multiple_testing: MultipleTestingConfig = field(
        default_factory=MultipleTestingConfig
    )

    def __post_init__(self) -> None:
        if self.paired_design_policy not in SUPPORTED_PAIRED_DESIGN_POLICIES:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_PAIRED_DESIGN_POLICIES
            )
            raise WorkflowValidationError(
                f"differential.paired_design_policy must be one of: {supported}"
            )
        object.__setattr__(
            self,
            "paired_design_policy",
            cast(PairedDesignPolicy, self.paired_design_policy),
        )


__all__ = [
    "DifferentialAnalysisConfig",
    "MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG",
    "MultipleTestingConfig",
    "PAIRED_DESIGN_POLICY_FIXED_BLOCK",
    "PAIRED_DESIGN_POLICY_REJECT",
    "PairedDesignPolicy",
    "SUPPORTED_PAIRED_DESIGN_POLICIES",
    "SUPPORTED_MULTIPLE_TESTING_METHODS",
]
