"""Public compatibility wrapper for internal contract ownership."""

from phospy.contracts.configs.differential import (
    MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
    SUPPORTED_MULTIPLE_TESTING_METHODS,
    SUPPORTED_PAIRED_DESIGN_POLICIES,
    DifferentialAnalysisConfig,
    MultipleTestingConfig,
    PairedDesignPolicy,
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
