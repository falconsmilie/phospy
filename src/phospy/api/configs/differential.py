"""Public compatibility wrapper for internal contract ownership."""

from phospy.contracts.configs.differential import (
    MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG,
    SUPPORTED_MULTIPLE_TESTING_METHODS,
    DifferentialAnalysisConfig,
    MultipleTestingConfig,
)

__all__ = [
    "DifferentialAnalysisConfig",
    "MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG",
    "MultipleTestingConfig",
    "SUPPORTED_MULTIPLE_TESTING_METHODS",
]
