"""Differential-analysis domain exports."""

from phospy.science.differential.aggregation import (
    PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE,
    PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META,
    PEPTIDE_TO_SITE_STRATEGY_INVERSE_VARIANCE_WEIGHTED,
    PEPTIDE_TO_SITE_STRATEGY_RANDOM_EFFECT_META,
    PEPTIDE_TO_SITE_STRATEGY_STOUFFER_Z,
    PeptideToSiteAggregationConfig,
    PeptideToSiteAggregationResult,
    PeptideToSiteAggregator,
)
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
    EmpiricalBayesConfig,
)
from phospy.science.differential.policy_models import TechnicalReplicatePolicy

__all__ = [
    "ContrastMatrix",
    "DesignMatrix",
    "DifferentialAnalysisRequest",
    "DifferentialAnalysisResult",
    "EmpiricalBayesConfig",
    "TechnicalReplicatePolicy",
    "PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE",
    "PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META",
    "PEPTIDE_TO_SITE_STRATEGY_INVERSE_VARIANCE_WEIGHTED",
    "PEPTIDE_TO_SITE_STRATEGY_RANDOM_EFFECT_META",
    "PEPTIDE_TO_SITE_STRATEGY_STOUFFER_Z",
    "PeptideToSiteAggregationConfig",
    "PeptideToSiteAggregationResult",
    "PeptideToSiteAggregator",
]
