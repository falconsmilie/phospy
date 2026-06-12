"""Experimental-design domain exports."""

from phospy.science.design.matrix_builder import (
    DesignMatrixBuilder,
    DesignMatrixBuildResult,
)
from phospy.science.design.models import (
    FIXED_EFFECT_COVARIATE_KIND_BATCH,
    FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
    FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
    SUPPORTED_FIXED_EFFECT_COVARIATE_KINDS,
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    ExperimentalDesign,
    FixedEffectCovariate,
    FixedEffectCovariateKind,
    SampleDesignRecord,
)

__all__ = [
    "BatchCovariate",
    "CategoricalCovariate",
    "Contrast",
    "ContinuousCovariate",
    "DesignMatrixBuilder",
    "DesignMatrixBuildResult",
    "ExperimentalDesign",
    "FIXED_EFFECT_COVARIATE_KIND_BATCH",
    "FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL",
    "FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS",
    "FixedEffectCovariate",
    "FixedEffectCovariateKind",
    "SUPPORTED_FIXED_EFFECT_COVARIATE_KINDS",
    "SampleDesignRecord",
]
