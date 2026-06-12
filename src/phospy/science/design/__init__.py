"""Experimental-design domain exports."""

from phospy.science.design.matrix_builder import (
    DesignMatrixBuilder,
    DesignMatrixBuildResult,
    describe_fixed_effect_design,
)
from phospy.science.design.models import (
    FIXED_EFFECT_COVARIATE_KIND_BATCH,
    FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
    FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
    SUPPORTED_FIXED_EFFECT_COVARIATE_KINDS,
    SUPPORTED_PAIRED_DESIGN_POLICIES,
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    ExperimentalDesign,
    FixedEffectCovariate,
    FixedEffectCovariateKind,
    PairedDesignPolicy,
    SampleDesignRecord,
)

__all__ = [
    "BatchCovariate",
    "CategoricalCovariate",
    "Contrast",
    "ContinuousCovariate",
    "DesignMatrixBuilder",
    "DesignMatrixBuildResult",
    "describe_fixed_effect_design",
    "ExperimentalDesign",
    "FIXED_EFFECT_COVARIATE_KIND_BATCH",
    "FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL",
    "FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS",
    "FixedEffectCovariate",
    "FixedEffectCovariateKind",
    "PAIRED_DESIGN_POLICY_FIXED_BLOCK",
    "PAIRED_DESIGN_POLICY_REJECT",
    "PairedDesignPolicy",
    "SUPPORTED_FIXED_EFFECT_COVARIATE_KINDS",
    "SUPPORTED_PAIRED_DESIGN_POLICIES",
    "SampleDesignRecord",
]
