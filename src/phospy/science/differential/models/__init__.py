"""Public models for differential analysis."""

from __future__ import annotations

from phospy.science.differential.models.design import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisRequest,
)
from phospy.science.differential.models.diagnostics import (
    DifferentialModelDiagnostics,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.models.empirical_bayes_config import (
    EMPIRICAL_BAYES_METHOD_ROBUST,
    EMPIRICAL_BAYES_METHOD_STANDARD,
    SUPPORTED_EMPIRICAL_BAYES_METHODS,
    EmpiricalBayesConfig,
)
from phospy.science.differential.models.fit import DifferentialComputationResult
from phospy.science.differential.models.provenance import (
    DifferentialContrastDefinition,
    DifferentialDesignMatrixSummary,
    DifferentialEmpiricalBayesProvenance,
    DifferentialFixedEffectCovariateProvenance,
    DifferentialMissingValuePolicyProvenance,
    DifferentialPolicyProvenance,
    DifferentialReplicatePolicyProvenance,
    DifferentialStatisticalTestingProvenance,
    DifferentialTechnicalReplicateGroup,
    DifferentialUnsupportedDesignPolicyProvenance,
)
from phospy.science.differential.models.results import DifferentialAnalysisResult
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_IMPUTATION_RESULT_COLUMNS,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
    DIFFERENTIAL_RESULT_WITHHELD_STATUSES,
)

__all__ = [
    "ContrastMatrix",
    "DesignMatrix",
    "DifferentialAnalysisRequest",
    "DifferentialAnalysisResult",
    "DifferentialComputationResult",
    "DifferentialContrastDefinition",
    "DifferentialDesignMatrixSummary",
    "DifferentialEmpiricalBayesProvenance",
    "DifferentialModelDiagnostics",
    "DifferentialFixedEffectCovariateProvenance",
    "DifferentialMissingValuePolicyProvenance",
    "DifferentialPolicyProvenance",
    "DifferentialReplicatePolicyProvenance",
    "DifferentialStatisticalTestingProvenance",
    "DifferentialTechnicalReplicateGroup",
    "DifferentialUnsupportedDesignPolicyProvenance",
    "DIFFERENTIAL_IMPUTATION_RESULT_COLUMNS",
    "DIFFERENTIAL_RESULT_STATUS_COLUMN",
    "DIFFERENTIAL_RESULT_STATUS_TESTED",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED",
    "DIFFERENTIAL_RESULT_WITHHELD_STATUSES",
    "EMPIRICAL_BAYES_METHOD_ROBUST",
    "EMPIRICAL_BAYES_METHOD_STANDARD",
    "EmpiricalBayesConfig",
    "EmpiricalBayesPriorDiagnostics",
    "MeanVarianceTrendDiagnostics",
    "SUPPORTED_EMPIRICAL_BAYES_METHODS",
]
