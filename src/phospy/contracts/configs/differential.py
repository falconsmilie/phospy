"""Public differential workflow configuration models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import cast

from phospy.errors.validation import ContractValidationError
from phospy.science.configs.differential import (
    DIFFERENTIAL_EXPLORATORY_MINIMUM_CONDITION_REPLICATES,
    DIFFERENTIAL_PRODUCTION_MINIMUM_CONDITION_REPLICATES,
    DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE,
    DIFFERENTIAL_RELIABILITY_PROFILE_PRODUCTION,
    IMPUTED_VALUE_POLICY_REJECT,
    IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
    MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG,
    MULTIPLE_TESTING_METHOD_BENJAMINI_YEKUTIELI,
    MULTIPLE_TESTING_METHOD_BONFERRONI,
    MULTIPLE_TESTING_METHOD_HOLM,
    MULTIPLE_TESTING_METHOD_NONE,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
    SUPPORTED_DIFFERENTIAL_IMPUTED_VALUE_POLICIES,
    SUPPORTED_DIFFERENTIAL_RELIABILITY_PROFILES,
    SUPPORTED_MULTIPLE_TESTING_METHODS,
    SUPPORTED_PAIRED_DESIGN_POLICIES,
    DifferentialImputedValuePolicy,
    DifferentialReliabilityProfile,
    MultipleTestingMethod,
    PairedDesignPolicy,
)
from phospy.science.differential.models import (
    EmpiricalBayesConfig,
)
from phospy.science.differential.policy_models import TechnicalReplicatePolicy


@dataclass(frozen=True, slots=True)
class MultipleTestingConfig:
    """Public multiple-testing policy for differential analysis."""

    method: MultipleTestingMethod = MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG

    def __post_init__(self) -> None:
        if self.method not in SUPPORTED_MULTIPLE_TESTING_METHODS:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_MULTIPLE_TESTING_METHODS
            )
            raise ContractValidationError(
                f"differential.multiple_testing.method must be one of: {supported}"
            )
        object.__setattr__(
            self,
            "method",
            cast(MultipleTestingMethod, self.method),
        )


@dataclass(frozen=True, slots=True)
class DifferentialAnalysisConfig:
    """Public configuration for differential analysis.

    Fixed-effect covariate declarations live on ``ExperimentalDesign``.
    ``paired_design_policy`` records user intent for explicit paired or blocked
    designs only. It does not infer ``block_id`` values and does not enable
    mixed-effects modelling. ``"fixed_block"`` requires complete block metadata
    and validates a fixed-effect block design matrix before execution.
    ``imputed_value_policy`` defaults to ``"reject"``. Non-default policies are
    explicit opt-ins and require dataset-owned imputation observation metadata.
    ``reliability_profile`` separates production-supported inference from the
    explicitly exploratory single-biological-replicate lane. Lowering
    ``minimum_condition_replicates`` is not a supported production override.
    ``allow_suspicious_declared_input_scale`` is an explicit scientific override
    for declared log2 input-scale provenance that recorded diagnostic warnings.
    """

    reliability_profile: DifferentialReliabilityProfile = (
        DIFFERENTIAL_RELIABILITY_PROFILE_PRODUCTION
    )
    technical_replicate_policy: TechnicalReplicatePolicy = (
        TechnicalReplicatePolicy.REJECT
    )
    paired_design_policy: PairedDesignPolicy = PAIRED_DESIGN_POLICY_REJECT
    imputed_value_policy: DifferentialImputedValuePolicy = IMPUTED_VALUE_POLICY_REJECT
    imputed_value_max_fraction: float = 0.0
    allow_design_subset: bool = False
    allow_suspicious_declared_input_scale: bool = False
    minimum_condition_replicates: int = 2
    empirical_bayes: EmpiricalBayesConfig = field(default_factory=EmpiricalBayesConfig)
    multiple_testing: MultipleTestingConfig = field(
        default_factory=MultipleTestingConfig
    )

    def __post_init__(self) -> None:
        if self.reliability_profile not in SUPPORTED_DIFFERENTIAL_RELIABILITY_PROFILES:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_DIFFERENTIAL_RELIABILITY_PROFILES
            )
            raise ContractValidationError(
                f"differential.reliability_profile must be one of: {supported}"
            )
        if self.paired_design_policy not in SUPPORTED_PAIRED_DESIGN_POLICIES:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_PAIRED_DESIGN_POLICIES
            )
            raise ContractValidationError(
                f"differential.paired_design_policy must be one of: {supported}"
            )
        if (
            self.imputed_value_policy
            not in SUPPORTED_DIFFERENTIAL_IMPUTED_VALUE_POLICIES
        ):
            supported = ", ".join(
                repr(value) for value in SUPPORTED_DIFFERENTIAL_IMPUTED_VALUE_POLICIES
            )
            raise ContractValidationError(
                f"differential.imputed_value_policy must be one of: {supported}"
            )
        if isinstance(cast(object, self.imputed_value_max_fraction), bool) or not (
            isinstance(cast(object, self.imputed_value_max_fraction), int | float)
        ):
            raise ContractValidationError(
                "differential.imputed_value_max_fraction must be a numeric value"
            )
        imputed_value_max_fraction = float(self.imputed_value_max_fraction)
        if not 0.0 <= imputed_value_max_fraction <= 1.0:
            raise ContractValidationError(
                "differential.imputed_value_max_fraction must be in [0.0, 1.0]"
            )
        object.__setattr__(
            self,
            "reliability_profile",
            cast(DifferentialReliabilityProfile, self.reliability_profile),
        )
        object.__setattr__(
            self,
            "paired_design_policy",
            cast(PairedDesignPolicy, self.paired_design_policy),
        )
        object.__setattr__(
            self,
            "imputed_value_policy",
            cast(DifferentialImputedValuePolicy, self.imputed_value_policy),
        )
        object.__setattr__(
            self,
            "imputed_value_max_fraction",
            imputed_value_max_fraction,
        )


__all__ = [
    "DIFFERENTIAL_EXPLORATORY_MINIMUM_CONDITION_REPLICATES",
    "DIFFERENTIAL_PRODUCTION_MINIMUM_CONDITION_REPLICATES",
    "DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE",
    "DIFFERENTIAL_RELIABILITY_PROFILE_PRODUCTION",
    "DifferentialImputedValuePolicy",
    "DifferentialAnalysisConfig",
    "DifferentialReliabilityProfile",
    "IMPUTED_VALUE_POLICY_REJECT",
    "IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES",
    "MULTIPLE_TESTING_METHOD_BENJAMINI_YEKUTIELI",
    "MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG",
    "MULTIPLE_TESTING_METHOD_BONFERRONI",
    "MULTIPLE_TESTING_METHOD_HOLM",
    "MULTIPLE_TESTING_METHOD_NONE",
    "MultipleTestingMethod",
    "MultipleTestingConfig",
    "PAIRED_DESIGN_POLICY_FIXED_BLOCK",
    "PAIRED_DESIGN_POLICY_REJECT",
    "PairedDesignPolicy",
    "SUPPORTED_DIFFERENTIAL_IMPUTED_VALUE_POLICIES",
    "SUPPORTED_DIFFERENTIAL_RELIABILITY_PROFILES",
    "SUPPORTED_PAIRED_DESIGN_POLICIES",
    "SUPPORTED_MULTIPLE_TESTING_METHODS",
]
