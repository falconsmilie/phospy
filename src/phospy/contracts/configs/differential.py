"""Public differential workflow configuration models."""

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.errors.validation import WorkflowValidationError
from phospy.science.differential.models import EmpiricalBayesConfig
from phospy.science.differential.policy_models import TechnicalReplicatePolicy

MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG = "benjamini_hochberg"
SUPPORTED_MULTIPLE_TESTING_METHODS: tuple[str, ...] = (
    MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG,
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

    Fixed-effect covariate declarations live on ``ExperimentalDesign``. This
    config does not enable adjusted-model execution in the current release.
    """

    technical_replicate_policy: TechnicalReplicatePolicy = (
        TechnicalReplicatePolicy.REJECT
    )
    allow_design_subset: bool = False
    minimum_condition_replicates: int = 2
    empirical_bayes: EmpiricalBayesConfig = field(default_factory=EmpiricalBayesConfig)
    multiple_testing: MultipleTestingConfig = field(
        default_factory=MultipleTestingConfig
    )


__all__ = [
    "DifferentialAnalysisConfig",
    "MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG",
    "MultipleTestingConfig",
    "SUPPORTED_MULTIPLE_TESTING_METHODS",
]
