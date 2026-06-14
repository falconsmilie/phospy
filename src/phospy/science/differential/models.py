"""Public models for differential analysis."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import (
    export_dataframe,
    export_series,
    own_dataframe,
    own_series,
)
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_string_index,
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.identity_contracts import (
    RESULT_IDENTITY_COLUMNS,
    RESULT_TABLE_IDENTITY_CONTRACT,
    enforce_phosphosite_identity_contract,
    enforce_required_identity_text_columns,
    enforce_result_identity_metadata_coherence,
)

if TYPE_CHECKING:
    from phospy.science.datasets.models import DatasetPreprocessingReport

EMPIRICAL_BAYES_METHOD_STANDARD = "standard"
EMPIRICAL_BAYES_METHOD_ROBUST = "robust"
SUPPORTED_EMPIRICAL_BAYES_METHODS: tuple[str, ...] = (
    EMPIRICAL_BAYES_METHOD_STANDARD,
    EMPIRICAL_BAYES_METHOD_ROBUST,
)
_RESULT_STATISTIC_COLUMNS: tuple[str, ...] = ("logFC", "t", "P.Value", "adj.P.Val")
_PUBLIC_RESULT_IDENTITY_COLUMNS: tuple[str, ...] = RESULT_IDENTITY_COLUMNS
DIFFERENTIAL_RESULT_STATUS_COLUMN = "result_status"
DIFFERENTIAL_RESULT_STATUS_TESTED = "tested"
DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION = "withheld_high_imputation"
DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED = (
    "withheld_insufficient_observed_samples"
)
DIFFERENTIAL_RESULT_WITHHELD_STATUSES: tuple[str, ...] = (
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
)
DIFFERENTIAL_IMPUTATION_RESULT_COLUMNS: tuple[str, ...] = (
    "imputed_cell_count",
    "observed_cell_count",
    "imputed_fraction",
    "imputation_policy",
    "imputation_fraction_threshold",
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
)


@dataclass(frozen=True, slots=True)
class EmpiricalBayesConfig:
    """Empirical-Bayes configuration for moderated statistics."""

    method: str = EMPIRICAL_BAYES_METHOD_STANDARD
    trend: bool = False
    winsor_tail_p: tuple[float, float] = (0.05, 0.1)

    def __post_init__(self) -> None:
        if self.method not in SUPPORTED_EMPIRICAL_BAYES_METHODS:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_EMPIRICAL_BAYES_METHODS
            )
            raise PhosPyInputError(
                f"empirical_bayes.method must be one of: {supported}"
            )
        if not isinstance(cast(object, self.trend), bool):
            raise PhosPyInputError("empirical_bayes.trend must be a bool")
        winsor_tail_p = self.winsor_tail_p
        if (
            not isinstance(cast(object, winsor_tail_p), tuple)
            or len(winsor_tail_p) != 2
            or not all(
                isinstance(cast(object, value), int | float) for value in winsor_tail_p
            )
        ):
            raise PhosPyInputError(
                "empirical_bayes.winsor_tail_p must be a tuple of two numeric values"
            )
        left_tail_p = float(winsor_tail_p[0])
        right_tail_p = float(winsor_tail_p[1])
        if not (0.0 <= left_tail_p < 1.0 and 0.0 <= right_tail_p < 1.0):
            raise PhosPyInputError(
                "empirical_bayes.winsor_tail_p values must each be in [0.0, 1.0)"
            )
        if left_tail_p + right_tail_p >= 1.0:
            raise PhosPyInputError(
                "empirical_bayes.winsor_tail_p values must sum to less than 1.0"
            )
        object.__setattr__(self, "winsor_tail_p", (left_tail_p, right_tail_p))


@dataclass(frozen=True, slots=True, init=False)
class EmpiricalBayesPriorDiagnostics:
    """Diagnostics for prior-variance and prior-df estimation."""

    method: str
    robust: bool
    trend: bool
    winsor_tail_p: tuple[float, float]
    base_prior_variance: float
    base_prior_degrees_of_freedom: float
    robust_outlier_count: int
    robust_outlier_fraction: float
    winsorized_low_count: int
    winsorized_high_count: int
    prior_variance: pd.Series
    prior_degrees_of_freedom: pd.Series

    def __init__(
        self,
        *,
        method: str,
        robust: bool,
        trend: bool,
        winsor_tail_p: tuple[float, float],
        base_prior_variance: float,
        base_prior_degrees_of_freedom: float,
        robust_outlier_count: int,
        robust_outlier_fraction: float,
        winsorized_low_count: int,
        winsorized_high_count: int,
        prior_variance: pd.Series,
        prior_degrees_of_freedom: pd.Series,
        _assume_owned: bool = False,
    ) -> None:
        prior_variance = own_series(
            prior_variance,
            field_name="differential_result.prior_diagnostics.prior_variance",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        prior_degrees_of_freedom = own_series(
            prior_degrees_of_freedom,
            field_name=(
                "differential_result.prior_diagnostics.prior_degrees_of_freedom"
            ),
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        if not prior_variance.index.equals(prior_degrees_of_freedom.index):
            raise PhosPyInputError(
                "differential_result prior diagnostics index mismatch between "
                "prior_variance and prior_degrees_of_freedom"
            )
        object.__setattr__(self, "method", str(method))
        object.__setattr__(self, "robust", bool(robust))
        object.__setattr__(self, "trend", bool(trend))
        object.__setattr__(
            self,
            "winsor_tail_p",
            (float(winsor_tail_p[0]), float(winsor_tail_p[1])),
        )
        object.__setattr__(self, "base_prior_variance", float(base_prior_variance))
        object.__setattr__(
            self,
            "base_prior_degrees_of_freedom",
            float(base_prior_degrees_of_freedom),
        )
        object.__setattr__(self, "robust_outlier_count", int(robust_outlier_count))
        object.__setattr__(
            self,
            "robust_outlier_fraction",
            float(robust_outlier_fraction),
        )
        object.__setattr__(self, "winsorized_low_count", int(winsorized_low_count))
        object.__setattr__(self, "winsorized_high_count", int(winsorized_high_count))
        object.__setattr__(self, "prior_variance", prior_variance)
        object.__setattr__(self, "prior_degrees_of_freedom", prior_degrees_of_freedom)

    def prior_variance_series(self) -> pd.Series:
        return export_series(self.prior_variance)

    def prior_degrees_of_freedom_series(self) -> pd.Series:
        return export_series(self.prior_degrees_of_freedom)


@dataclass(frozen=True, slots=True, init=False)
class MeanVarianceTrendDiagnostics:
    """Diagnostics payload for mean-intensity vs variance trend fitting."""

    mean_intensity: pd.Series
    log_residual_variance: pd.Series
    fitted_log_prior_variance: pd.Series

    def __init__(
        self,
        *,
        mean_intensity: pd.Series,
        log_residual_variance: pd.Series,
        fitted_log_prior_variance: pd.Series,
        _assume_owned: bool = False,
    ) -> None:
        mean_intensity = own_series(
            mean_intensity,
            field_name="differential_result.mean_variance_trend.mean_intensity",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        log_residual_variance = own_series(
            log_residual_variance,
            field_name=(
                "differential_result.mean_variance_trend.log_residual_variance"
            ),
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        fitted_log_prior_variance = own_series(
            fitted_log_prior_variance,
            field_name=(
                "differential_result.mean_variance_trend.fitted_log_prior_variance"
            ),
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        if not mean_intensity.index.equals(log_residual_variance.index):
            raise PhosPyInputError(
                "mean-variance trend diagnostics index mismatch for mean_intensity and "
                "log_residual_variance"
            )
        if not mean_intensity.index.equals(fitted_log_prior_variance.index):
            raise PhosPyInputError(
                "mean-variance trend diagnostics index mismatch for mean_intensity and "
                "fitted_log_prior_variance"
            )
        object.__setattr__(self, "mean_intensity", mean_intensity)
        object.__setattr__(self, "log_residual_variance", log_residual_variance)
        object.__setattr__(
            self,
            "fitted_log_prior_variance",
            fitted_log_prior_variance,
        )

    def mean_intensity_series(self) -> pd.Series:
        return export_series(self.mean_intensity)

    def log_residual_variance_series(self) -> pd.Series:
        return export_series(self.log_residual_variance)

    def fitted_log_prior_variance_series(self) -> pd.Series:
        return export_series(self.fitted_log_prior_variance)


@dataclass(frozen=True, slots=True)
class DesignMatrix:
    """Validated design matrix with samples on rows and coefficients on columns."""

    frame: pd.DataFrame

    def __post_init__(self) -> None:
        frame = own_dataframe(
            self.frame,
            field_name="differential.design",
            error_type=PhosPyInputError,
        )
        _validate_numeric_matrix(
            frame,
            field_name="differential.design",
        )
        object.__setattr__(self, "frame", frame)

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self.frame)


@dataclass(frozen=True, slots=True)
class ContrastMatrix:
    """Validated contrast matrix with design coefficients on rows."""

    frame: pd.DataFrame

    def __post_init__(self) -> None:
        frame = own_dataframe(
            self.frame,
            field_name="differential.contrasts",
            error_type=PhosPyInputError,
        )
        _validate_numeric_matrix(
            frame,
            field_name="differential.contrasts",
        )
        object.__setattr__(self, "frame", frame)

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self.frame)


@dataclass(frozen=True, slots=True)
class DifferentialAnalysisRequest:
    """Request payload for limma-style moderated differential analysis."""

    matrix: pd.DataFrame
    design: DesignMatrix | pd.DataFrame
    contrasts: ContrastMatrix | pd.DataFrame
    empirical_bayes: EmpiricalBayesConfig = field(default_factory=EmpiricalBayesConfig)

    def __post_init__(self) -> None:
        matrix = own_dataframe(
            self.matrix,
            field_name="differential.matrix",
            error_type=PhosPyInputError,
        )
        _validate_numeric_matrix(
            matrix,
            field_name="differential.matrix",
        )
        design = self.design
        if isinstance(design, pd.DataFrame):
            design = DesignMatrix(design)
        if not isinstance(cast(object, design), DesignMatrix):
            raise PhosPyInputError(
                "differential.design must be a DesignMatrix or pandas DataFrame"
            )
        contrasts = self.contrasts
        if isinstance(contrasts, pd.DataFrame):
            contrasts = ContrastMatrix(contrasts)
        if not isinstance(cast(object, contrasts), ContrastMatrix):
            raise PhosPyInputError(
                "differential.contrasts must be a ContrastMatrix or pandas DataFrame"
            )
        if not isinstance(cast(object, self.empirical_bayes), EmpiricalBayesConfig):
            raise PhosPyInputError(
                "differential.empirical_bayes must be an EmpiricalBayesConfig"
            )
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "design", design)
        object.__setattr__(self, "contrasts", contrasts)


@dataclass(frozen=True, slots=True)
class DifferentialFixedEffectCovariateProvenance:
    """Resolved fixed-effect covariate columns included in the fitted design."""

    name: str
    kind: str
    columns: tuple[str, ...]
    levels: tuple[str, ...] = ()
    reference_level: str | None = None
    unused_levels: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name:
            raise PhosPyInputError(
                "differential_policy_provenance.design.covariates[].name must be "
                "non-empty"
            )
        if not self.kind:
            raise PhosPyInputError(
                "differential_policy_provenance.design.covariates[].kind must be "
                "non-empty"
            )
        if not self.columns:
            raise PhosPyInputError(
                "differential_policy_provenance.design.covariates[].columns must be "
                "non-empty"
            )
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(self, "kind", str(self.kind))
        object.__setattr__(self, "columns", tuple(str(value) for value in self.columns))
        object.__setattr__(self, "levels", tuple(str(value) for value in self.levels))
        object.__setattr__(
            self,
            "reference_level",
            None if self.reference_level is None else str(self.reference_level),
        )
        object.__setattr__(
            self,
            "unused_levels",
            tuple(str(value) for value in self.unused_levels),
        )


@dataclass(frozen=True, slots=True)
class DifferentialDesignMatrixSummary:
    """Structured summary of the resolved differential design matrix."""

    formula: str
    sample_labels: tuple[str, ...]
    coefficient_labels: tuple[str, ...]
    sample_count: int
    coefficient_count: int
    rank: int
    residual_degrees_of_freedom: float
    description: str = ""
    condition_columns: tuple[str, ...] = ()
    covariates: tuple[DifferentialFixedEffectCovariateProvenance, ...] = ()
    paired_design_policy: str = "reject"
    block_id_field_name: str = "block_id"
    block_count: int = 0
    block_levels: tuple[str, ...] = ()
    block_levels_included: tuple[str, ...] = ()
    block_reference_level: str | None = None
    block_columns: tuple[tuple[str, str], ...] = ()
    block_column_names: tuple[str, ...] = ()
    condition_coverage_rule: str = ""
    limitations: tuple[str, ...] = ()
    rank_validation_status: str = "not_recorded"
    estimability_validation_status: str = "not_recorded"

    def __post_init__(self) -> None:
        if not self.formula:
            raise PhosPyInputError(
                "differential_policy_provenance.design.formula must be non-empty"
            )
        if not self.sample_labels:
            raise PhosPyInputError(
                "differential_policy_provenance.design.sample_labels must be non-empty"
            )
        if not self.coefficient_labels:
            raise PhosPyInputError(
                "differential_policy_provenance.design.coefficient_labels must be "
                "non-empty"
            )
        if self.sample_count < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.design.sample_count must be >= 1"
            )
        if self.coefficient_count < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.design.coefficient_count must be >= 1"
            )
        if self.rank < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.design.rank must be >= 1"
            )
        if self.rank > self.coefficient_count:
            raise PhosPyInputError(
                "differential_policy_provenance.design.rank cannot exceed "
                "coefficient_count"
            )
        if self.residual_degrees_of_freedom <= 0.0:
            raise PhosPyInputError(
                "differential_policy_provenance.design.residual_degrees_of_freedom "
                "must be > 0.0"
            )
        if self.block_count < 0:
            raise PhosPyInputError(
                "differential_policy_provenance.design.block_count must be >= 0"
            )
        if not self.rank_validation_status:
            raise PhosPyInputError(
                "differential_policy_provenance.design.rank_validation_status must "
                "be non-empty"
            )
        if not self.estimability_validation_status:
            raise PhosPyInputError(
                "differential_policy_provenance.design."
                "estimability_validation_status must be non-empty"
            )
        covariates = tuple(self.covariates)
        for covariate in covariates:
            if not isinstance(
                cast(object, covariate),
                DifferentialFixedEffectCovariateProvenance,
            ):
                raise PhosPyInputError(
                    "differential_policy_provenance.design.covariates must contain "
                    "DifferentialFixedEffectCovariateProvenance values"
                )
        object.__setattr__(self, "formula", str(self.formula))
        object.__setattr__(self, "description", str(self.description))
        object.__setattr__(
            self,
            "sample_labels",
            tuple(str(value) for value in self.sample_labels),
        )
        object.__setattr__(
            self,
            "coefficient_labels",
            tuple(str(value) for value in self.coefficient_labels),
        )
        object.__setattr__(
            self,
            "condition_columns",
            tuple(str(value) for value in self.condition_columns),
        )
        object.__setattr__(self, "covariates", covariates)
        object.__setattr__(
            self,
            "paired_design_policy",
            str(self.paired_design_policy),
        )
        object.__setattr__(
            self,
            "block_id_field_name",
            str(self.block_id_field_name),
        )
        object.__setattr__(self, "block_count", int(self.block_count))
        object.__setattr__(
            self,
            "block_levels",
            tuple(str(value) for value in self.block_levels),
        )
        object.__setattr__(
            self,
            "block_levels_included",
            tuple(str(value) for value in self.block_levels_included),
        )
        object.__setattr__(
            self,
            "block_reference_level",
            (
                None
                if self.block_reference_level is None
                else str(self.block_reference_level)
            ),
        )
        object.__setattr__(
            self,
            "block_columns",
            tuple((str(level), str(column)) for level, column in self.block_columns),
        )
        object.__setattr__(
            self,
            "block_column_names",
            tuple(str(value) for value in self.block_column_names),
        )
        object.__setattr__(
            self,
            "condition_coverage_rule",
            str(self.condition_coverage_rule),
        )
        object.__setattr__(
            self,
            "limitations",
            tuple(str(value) for value in self.limitations),
        )
        object.__setattr__(
            self,
            "rank_validation_status",
            str(self.rank_validation_status),
        )
        object.__setattr__(
            self,
            "estimability_validation_status",
            str(self.estimability_validation_status),
        )


@dataclass(frozen=True, slots=True)
class DifferentialContrastDefinition:
    """Structured differential contrast definition."""

    name: str
    numerator_condition: str
    denominator_condition: str
    coefficients: tuple[tuple[str, float], ...]
    description: str = ""

    def __post_init__(self) -> None:
        if not self.name:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts[].name must be non-empty"
            )
        if not self.numerator_condition:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts[].numerator_condition must "
                "be non-empty"
            )
        if not self.denominator_condition:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts[].denominator_condition "
                "must be non-empty"
            )
        if not self.coefficients:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts[].coefficients must be "
                "non-empty"
            )
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(
            self,
            "numerator_condition",
            str(self.numerator_condition),
        )
        object.__setattr__(
            self,
            "denominator_condition",
            str(self.denominator_condition),
        )
        object.__setattr__(
            self,
            "coefficients",
            tuple((str(name), float(value)) for name, value in self.coefficients),
        )
        object.__setattr__(self, "description", str(self.description))


@dataclass(frozen=True, slots=True)
class DifferentialTechnicalReplicateGroup:
    """Structured technical-replicate lineage for one resolved group."""

    condition: str
    biological_replicate_id: str
    output_sample_id: str
    input_sample_ids: tuple[str, ...]
    technical_replicate_ids: tuple[str, ...]
    n_technical_replicates: int

    def __post_init__(self) -> None:
        if not self.condition:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".condition must be non-empty"
            )
        if not self.biological_replicate_id:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".biological_replicate_id must be non-empty"
            )
        if not self.output_sample_id:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".output_sample_id must be non-empty"
            )
        if not self.input_sample_ids:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".input_sample_ids must be non-empty"
            )
        if self.n_technical_replicates < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".n_technical_replicates must be >= 1"
            )


@dataclass(frozen=True, slots=True)
class DifferentialReplicatePolicyProvenance:
    """Structured replicate/group requirements for differential analysis."""

    minimum_condition_replicates: int
    technical_replicate_policy: str
    condition_replicate_counts: tuple[tuple[str, int], ...]
    technical_replicate_groups: tuple[DifferentialTechnicalReplicateGroup, ...] = ()

    def __post_init__(self) -> None:
        if self.minimum_condition_replicates < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates."
                "minimum_condition_replicates must be >= 1"
            )
        if not self.technical_replicate_policy:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates."
                "technical_replicate_policy must be non-empty"
            )
        if not self.condition_replicate_counts:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.condition_replicate_counts "
                "must be non-empty"
            )


@dataclass(frozen=True, slots=True)
class DifferentialEmpiricalBayesProvenance:
    """Structured empirical-Bayes moderation settings."""

    method: str
    robust: bool
    trend: bool
    winsor_tail_p: tuple[float, float]

    def __post_init__(self) -> None:
        if not self.method:
            raise PhosPyInputError(
                "differential_policy_provenance.empirical_bayes.method must be "
                "non-empty"
            )


@dataclass(frozen=True, slots=True)
class DifferentialStatisticalTestingProvenance:
    """Structured p-value and multiple-testing adjustment settings."""

    test_statistic: str
    p_value_method: str
    adjusted_p_value_method: str

    def __post_init__(self) -> None:
        if not self.test_statistic:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing.test_statistic "
                "must be non-empty"
            )
        if not self.p_value_method:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing.p_value_method "
                "must be non-empty"
            )
        if not self.adjusted_p_value_method:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing."
                "adjusted_p_value_method must be non-empty"
            )


@dataclass(frozen=True, slots=True)
class DifferentialMissingValuePolicyProvenance:
    """Structured missing-value handling policy for differential execution."""

    policy: str
    stage: str
    imputed_value_policy: str = "reject"
    imputed_value_max_fraction: float = 0.0
    imputation_metadata_required: bool = False
    adjusted_p_value_scope: str = "all_tested_features"
    limitations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.policy:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values.policy must be non-empty"
            )
        if not self.stage:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values.stage must be non-empty"
            )
        if not self.imputed_value_policy:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values."
                "imputed_value_policy must be non-empty"
            )
        object.__setattr__(self, "policy", str(self.policy))
        object.__setattr__(self, "stage", str(self.stage))
        object.__setattr__(
            self,
            "imputed_value_policy",
            str(self.imputed_value_policy),
        )
        object.__setattr__(
            self,
            "imputed_value_max_fraction",
            float(self.imputed_value_max_fraction),
        )
        object.__setattr__(
            self,
            "imputation_metadata_required",
            bool(self.imputation_metadata_required),
        )
        object.__setattr__(
            self,
            "adjusted_p_value_scope",
            str(self.adjusted_p_value_scope),
        )
        object.__setattr__(
            self,
            "limitations",
            tuple(str(value) for value in self.limitations),
        )


@dataclass(frozen=True, slots=True)
class DifferentialUnsupportedDesignPolicyProvenance:
    """Structured record of unsupported differential-design features."""

    intentionally_rejected_features: tuple[str, ...]
    enforcement_stage: str
    policy: str = "reject_unsupported_design_features_before_execution"

    def __post_init__(self) -> None:
        if not self.intentionally_rejected_features:
            raise PhosPyInputError(
                "differential_policy_provenance.unsupported_design."
                "intentionally_rejected_features must be non-empty"
            )
        if not self.enforcement_stage:
            raise PhosPyInputError(
                "differential_policy_provenance.unsupported_design.enforcement_stage "
                "must be non-empty"
            )
        if not self.policy:
            raise PhosPyInputError(
                "differential_policy_provenance.unsupported_design.policy must be "
                "non-empty"
            )
        object.__setattr__(
            self,
            "intentionally_rejected_features",
            tuple(str(value) for value in self.intentionally_rejected_features),
        )
        object.__setattr__(self, "enforcement_stage", str(self.enforcement_stage))
        object.__setattr__(self, "policy", str(self.policy))


@dataclass(frozen=True, slots=True)
class DifferentialPolicyProvenance:
    """Structured differential-analysis statistical policy provenance."""

    design: DifferentialDesignMatrixSummary
    contrasts: tuple[DifferentialContrastDefinition, ...]
    replicates: DifferentialReplicatePolicyProvenance
    empirical_bayes: DifferentialEmpiricalBayesProvenance
    statistical_testing: DifferentialStatisticalTestingProvenance
    missing_values: DifferentialMissingValuePolicyProvenance
    unsupported_design: DifferentialUnsupportedDesignPolicyProvenance

    def __post_init__(self) -> None:
        if not self.contrasts:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts must be non-empty"
            )


@dataclass(frozen=True, slots=True, init=False)
class DifferentialComputationResult:
    """Internal stat-only output from differential model computation.

    Contrast tables are indexed by the input matrix feature index and contain
    moderated statistic columns. Workflow layers attach biological identity
    metadata before constructing the public ``DifferentialAnalysisResult``.
    """

    residual_variance: pd.Series
    posterior_residual_variance: pd.Series
    prior_residual_variance: pd.Series
    prior_degrees_of_freedom_series_value: pd.Series
    prior_variance: float
    prior_degrees_of_freedom: float
    residual_degrees_of_freedom: float
    empirical_bayes_method: str
    empirical_bayes_robust: bool
    empirical_bayes_trend: bool
    prior_diagnostics: EmpiricalBayesPriorDiagnostics
    mean_variance_trend_diagnostics: MeanVarianceTrendDiagnostics | None
    _contrast_tables: Mapping[str, pd.DataFrame]

    def __init__(
        self,
        *,
        residual_variance: pd.Series,
        posterior_residual_variance: pd.Series,
        prior_residual_variance: pd.Series,
        prior_degrees_of_freedom_series_value: pd.Series,
        prior_variance: float,
        prior_degrees_of_freedom: float,
        residual_degrees_of_freedom: float,
        empirical_bayes_method: str,
        empirical_bayes_robust: bool,
        empirical_bayes_trend: bool,
        prior_diagnostics: EmpiricalBayesPriorDiagnostics,
        mean_variance_trend_diagnostics: MeanVarianceTrendDiagnostics | None,
        contrast_tables: Mapping[str, pd.DataFrame],
        _assume_owned: bool = False,
    ) -> None:
        residual_variance = own_series(
            residual_variance,
            field_name="differential_computation_result.residual_variance",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        posterior_residual_variance = own_series(
            posterior_residual_variance,
            field_name="differential_computation_result.posterior_residual_variance",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        prior_residual_variance = own_series(
            prior_residual_variance,
            field_name="differential_computation_result.prior_residual_variance",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        prior_degrees_of_freedom_series_value = own_series(
            prior_degrees_of_freedom_series_value,
            field_name=(
                "differential_computation_result.prior_degrees_of_freedom_series"
            ),
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        if not residual_variance.index.equals(posterior_residual_variance.index):
            raise PhosPyInputError(
                "differential_computation_result.posterior_residual_variance index "
                "must match differential_computation_result.residual_variance index"
            )
        if not residual_variance.index.equals(prior_residual_variance.index):
            raise PhosPyInputError(
                "differential_computation_result.prior_residual_variance index must "
                "match differential_computation_result.residual_variance index"
            )
        if not residual_variance.index.equals(
            prior_degrees_of_freedom_series_value.index
        ):
            raise PhosPyInputError(
                "differential_computation_result.prior_degrees_of_freedom_series "
                "index must match differential_computation_result.residual_variance "
                "index"
            )
        if not prior_diagnostics.prior_variance.index.equals(residual_variance.index):
            raise PhosPyInputError(
                "differential_computation_result.prior_diagnostics.prior_variance "
                "index must match matrix feature index"
            )
        if not prior_diagnostics.prior_degrees_of_freedom.index.equals(
            residual_variance.index
        ):
            raise PhosPyInputError(
                "differential_computation_result.prior_diagnostics."
                "prior_degrees_of_freedom index must match matrix feature index"
            )
        if (
            mean_variance_trend_diagnostics is not None
            and not mean_variance_trend_diagnostics.mean_intensity.index.equals(
                residual_variance.index
            )
        ):
            raise PhosPyInputError(
                "differential_computation_result.mean_variance_trend_diagnostics "
                "index must match matrix feature index"
            )
        if not contrast_tables:
            raise PhosPyInputError(
                "differential_computation_result.contrast_tables must include at "
                "least one contrast"
            )

        owned_tables: dict[str, pd.DataFrame] = {}
        for contrast_name, table in contrast_tables.items():
            if not isinstance(cast(object, contrast_name), str) or not contrast_name:
                raise PhosPyInputError(
                    "differential_computation_result.contrast_tables keys must be "
                    "non-empty strings"
                )
            owned_table = own_dataframe(
                table,
                field_name=(
                    "differential_computation_result.contrast_tables"
                    f"[{contrast_name!r}]"
                ),
                error_type=PhosPyInputError,
                assume_owned=_assume_owned,
            )
            _validate_computation_result_table(
                owned_table,
                field_name=(
                    "differential_computation_result.contrast_tables"
                    f"[{contrast_name!r}]"
                ),
            )
            if not owned_table.index.equals(residual_variance.index):
                raise PhosPyInputError(
                    "differential computation result table index must match matrix "
                    "feature index"
                )
            owned_tables[contrast_name] = owned_table

        object.__setattr__(self, "residual_variance", residual_variance)
        object.__setattr__(
            self,
            "posterior_residual_variance",
            posterior_residual_variance,
        )
        object.__setattr__(
            self,
            "prior_residual_variance",
            prior_residual_variance,
        )
        object.__setattr__(
            self,
            "prior_degrees_of_freedom_series_value",
            prior_degrees_of_freedom_series_value,
        )
        object.__setattr__(self, "prior_variance", float(prior_variance))
        object.__setattr__(
            self,
            "prior_degrees_of_freedom",
            float(prior_degrees_of_freedom),
        )
        object.__setattr__(
            self,
            "residual_degrees_of_freedom",
            float(residual_degrees_of_freedom),
        )
        object.__setattr__(self, "empirical_bayes_method", str(empirical_bayes_method))
        object.__setattr__(self, "empirical_bayes_robust", bool(empirical_bayes_robust))
        object.__setattr__(self, "empirical_bayes_trend", bool(empirical_bayes_trend))
        object.__setattr__(self, "prior_diagnostics", prior_diagnostics)
        object.__setattr__(
            self,
            "mean_variance_trend_diagnostics",
            mean_variance_trend_diagnostics,
        )
        object.__setattr__(self, "_contrast_tables", owned_tables)

    @property
    def contrast_tables(self) -> dict[str, pd.DataFrame]:
        return {
            contrast_name: export_dataframe(table)
            for contrast_name, table in self._contrast_tables.items()
        }

    def table_for(self, contrast_name: str) -> pd.DataFrame:
        if contrast_name not in self._contrast_tables:
            available = ", ".join(sorted(self._contrast_tables))
            raise KeyError(
                f"unknown contrast {contrast_name!r}; available: {available}"
            )
        return export_dataframe(self._contrast_tables[contrast_name])

    def residual_variance_series(self) -> pd.Series:
        return export_series(self.residual_variance)

    def posterior_residual_variance_series(self) -> pd.Series:
        return export_series(self.posterior_residual_variance)

    def prior_residual_variance_series(self) -> pd.Series:
        return export_series(self.prior_residual_variance)

    def prior_degrees_of_freedom_series(self) -> pd.Series:
        return export_series(self.prior_degrees_of_freedom_series_value)

    @classmethod
    def _from_owned(
        cls,
        *,
        residual_variance: pd.Series,
        posterior_residual_variance: pd.Series,
        prior_residual_variance: pd.Series,
        prior_degrees_of_freedom_series_value: pd.Series,
        prior_variance: float,
        prior_degrees_of_freedom: float,
        residual_degrees_of_freedom: float,
        empirical_bayes_method: str,
        empirical_bayes_robust: bool,
        empirical_bayes_trend: bool,
        prior_diagnostics: EmpiricalBayesPriorDiagnostics,
        mean_variance_trend_diagnostics: MeanVarianceTrendDiagnostics | None,
        contrast_tables: Mapping[str, pd.DataFrame],
    ) -> DifferentialComputationResult:
        return cls(
            residual_variance=residual_variance,
            posterior_residual_variance=posterior_residual_variance,
            prior_residual_variance=prior_residual_variance,
            prior_degrees_of_freedom_series_value=prior_degrees_of_freedom_series_value,
            prior_variance=prior_variance,
            prior_degrees_of_freedom=prior_degrees_of_freedom,
            residual_degrees_of_freedom=residual_degrees_of_freedom,
            empirical_bayes_method=empirical_bayes_method,
            empirical_bayes_robust=empirical_bayes_robust,
            empirical_bayes_trend=empirical_bayes_trend,
            prior_diagnostics=prior_diagnostics,
            mean_variance_trend_diagnostics=mean_variance_trend_diagnostics,
            contrast_tables=contrast_tables,
            _assume_owned=True,
        )


@dataclass(frozen=True, slots=True, init=False)
class DifferentialAnalysisResult:
    """Differential-analysis output with per-contrast moderated tables.

    Public contrast tables are indexed by protein-scoped ``site_key`` values and
    must include ``site_key``, ``display_id``, ``gene_symbol``, and ``site``
    columns.
    """

    residual_variance: pd.Series
    posterior_residual_variance: pd.Series
    prior_residual_variance: pd.Series
    prior_degrees_of_freedom_series_value: pd.Series
    prior_variance: float
    prior_degrees_of_freedom: float
    residual_degrees_of_freedom: float
    empirical_bayes_method: str
    empirical_bayes_robust: bool
    empirical_bayes_trend: bool
    prior_diagnostics: EmpiricalBayesPriorDiagnostics
    mean_variance_trend_diagnostics: MeanVarianceTrendDiagnostics | None
    policy_provenance: DifferentialPolicyProvenance | None
    workflow_provenance: Mapping[str, object] | None
    input_dataset_preprocessing_report: DatasetPreprocessingReport | None
    _contrast_tables: Mapping[str, pd.DataFrame]

    def __init__(
        self,
        *,
        residual_variance: pd.Series,
        posterior_residual_variance: pd.Series,
        prior_residual_variance: pd.Series,
        prior_degrees_of_freedom_series_value: pd.Series,
        prior_variance: float,
        prior_degrees_of_freedom: float,
        residual_degrees_of_freedom: float,
        empirical_bayes_method: str,
        empirical_bayes_robust: bool,
        empirical_bayes_trend: bool,
        prior_diagnostics: EmpiricalBayesPriorDiagnostics,
        mean_variance_trend_diagnostics: MeanVarianceTrendDiagnostics | None,
        contrast_tables: Mapping[str, pd.DataFrame],
        policy_provenance: DifferentialPolicyProvenance | None = None,
        workflow_provenance: Mapping[str, object] | None = None,
        input_dataset_preprocessing_report: DatasetPreprocessingReport | None = None,
        _assume_owned: bool = False,
    ) -> None:
        residual_variance = own_series(
            residual_variance,
            field_name="differential_result.residual_variance",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        posterior_residual_variance = own_series(
            posterior_residual_variance,
            field_name="differential_result.posterior_residual_variance",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        prior_residual_variance = own_series(
            prior_residual_variance,
            field_name="differential_result.prior_residual_variance",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        prior_degrees_of_freedom_series_value = own_series(
            prior_degrees_of_freedom_series_value,
            field_name="differential_result.prior_degrees_of_freedom_series",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        if not residual_variance.index.equals(posterior_residual_variance.index):
            raise PhosPyInputError(
                "differential_result.posterior_residual_variance index must match "
                "differential_result.residual_variance index"
            )
        if not residual_variance.index.equals(prior_residual_variance.index):
            raise PhosPyInputError(
                "differential_result.prior_residual_variance index must match "
                "differential_result.residual_variance index"
            )
        if not residual_variance.index.equals(
            prior_degrees_of_freedom_series_value.index
        ):
            raise PhosPyInputError(
                "differential_result.prior_degrees_of_freedom_series index must match "
                "differential_result.residual_variance index"
            )
        if not prior_diagnostics.prior_variance.index.equals(residual_variance.index):
            raise PhosPyInputError(
                "differential_result.prior_diagnostics.prior_variance index must match "
                "matrix feature index"
            )
        if not prior_diagnostics.prior_degrees_of_freedom.index.equals(
            residual_variance.index
        ):
            raise PhosPyInputError(
                "differential_result.prior_diagnostics.prior_degrees_of_freedom index "
                "must match matrix feature index"
            )
        if (
            mean_variance_trend_diagnostics is not None
            and not mean_variance_trend_diagnostics.mean_intensity.index.equals(
                residual_variance.index
            )
        ):
            raise PhosPyInputError(
                "differential_result.mean_variance_trend_diagnostics index must match "
                "matrix feature index"
            )
        if policy_provenance is not None and not isinstance(
            cast(object, policy_provenance), DifferentialPolicyProvenance
        ):
            raise PhosPyInputError(
                "differential_result.policy_provenance must be "
                "DifferentialPolicyProvenance or None"
            )
        if not contrast_tables:
            raise PhosPyInputError(
                "differential_result.contrast_tables must include at least one contrast"
            )
        if workflow_provenance is not None and not isinstance(
            cast(object, workflow_provenance), Mapping
        ):
            raise PhosPyInputError(
                "differential_result.workflow_provenance must be a mapping or None"
            )
        if (
            input_dataset_preprocessing_report is not None
            and not _is_dataset_preprocessing_report(input_dataset_preprocessing_report)
        ):
            raise PhosPyInputError(
                "differential_result.input_dataset_preprocessing_report must be "
                "DatasetPreprocessingReport or None"
            )
        owned_tables: dict[str, pd.DataFrame] = {}
        for contrast_name, table in contrast_tables.items():
            if not isinstance(cast(object, contrast_name), str) or not contrast_name:
                raise PhosPyInputError(
                    "differential_result.contrast_tables keys must be non-empty strings"
                )
            owned_table = own_dataframe(
                table,
                field_name=f"differential_result.contrast_tables[{contrast_name!r}]",
                error_type=PhosPyInputError,
                assume_owned=_assume_owned,
            )
            _validate_result_table(
                owned_table,
                field_name=f"differential_result.contrast_tables[{contrast_name!r}]",
            )
            if not owned_table.index.equals(residual_variance.index):
                raise PhosPyInputError(
                    "differential result table index must match matrix feature index"
                )
            owned_tables[contrast_name] = owned_table
        object.__setattr__(self, "residual_variance", residual_variance)
        object.__setattr__(
            self,
            "posterior_residual_variance",
            posterior_residual_variance,
        )
        object.__setattr__(
            self,
            "prior_residual_variance",
            prior_residual_variance,
        )
        object.__setattr__(
            self,
            "prior_degrees_of_freedom_series_value",
            prior_degrees_of_freedom_series_value,
        )
        object.__setattr__(self, "prior_variance", float(prior_variance))
        object.__setattr__(
            self,
            "prior_degrees_of_freedom",
            float(prior_degrees_of_freedom),
        )
        object.__setattr__(
            self,
            "residual_degrees_of_freedom",
            float(residual_degrees_of_freedom),
        )
        object.__setattr__(self, "empirical_bayes_method", str(empirical_bayes_method))
        object.__setattr__(self, "empirical_bayes_robust", bool(empirical_bayes_robust))
        object.__setattr__(self, "empirical_bayes_trend", bool(empirical_bayes_trend))
        object.__setattr__(self, "prior_diagnostics", prior_diagnostics)
        object.__setattr__(
            self,
            "mean_variance_trend_diagnostics",
            mean_variance_trend_diagnostics,
        )
        object.__setattr__(self, "policy_provenance", policy_provenance)
        object.__setattr__(
            self,
            "workflow_provenance",
            (
                None
                if workflow_provenance is None
                else {str(key): value for key, value in workflow_provenance.items()}
            ),
        )
        object.__setattr__(
            self,
            "input_dataset_preprocessing_report",
            input_dataset_preprocessing_report,
        )
        object.__setattr__(self, "_contrast_tables", owned_tables)

    @property
    def contrast_tables(self) -> dict[str, pd.DataFrame]:
        return {
            contrast_name: export_dataframe(table)
            for contrast_name, table in self._contrast_tables.items()
        }

    def table_for(self, contrast_name: str) -> pd.DataFrame:
        if contrast_name not in self._contrast_tables:
            available = ", ".join(sorted(self._contrast_tables))
            raise KeyError(
                f"unknown contrast {contrast_name!r}; available: {available}"
            )
        return export_dataframe(self._contrast_tables[contrast_name])

    def residual_variance_series(self) -> pd.Series:
        return export_series(self.residual_variance)

    def posterior_residual_variance_series(self) -> pd.Series:
        return export_series(self.posterior_residual_variance)

    def prior_residual_variance_series(self) -> pd.Series:
        return export_series(self.prior_residual_variance)

    def prior_degrees_of_freedom_series(self) -> pd.Series:
        return export_series(self.prior_degrees_of_freedom_series_value)

    @classmethod
    def _from_owned(
        cls,
        *,
        residual_variance: pd.Series,
        posterior_residual_variance: pd.Series,
        prior_residual_variance: pd.Series,
        prior_degrees_of_freedom_series_value: pd.Series,
        prior_variance: float,
        prior_degrees_of_freedom: float,
        residual_degrees_of_freedom: float,
        empirical_bayes_method: str,
        empirical_bayes_robust: bool,
        empirical_bayes_trend: bool,
        prior_diagnostics: EmpiricalBayesPriorDiagnostics,
        mean_variance_trend_diagnostics: MeanVarianceTrendDiagnostics | None,
        contrast_tables: Mapping[str, pd.DataFrame],
        policy_provenance: DifferentialPolicyProvenance | None = None,
        workflow_provenance: Mapping[str, object] | None = None,
        input_dataset_preprocessing_report: DatasetPreprocessingReport | None = None,
    ) -> DifferentialAnalysisResult:
        return cls(
            residual_variance=residual_variance,
            posterior_residual_variance=posterior_residual_variance,
            prior_residual_variance=prior_residual_variance,
            prior_degrees_of_freedom_series_value=prior_degrees_of_freedom_series_value,
            prior_variance=prior_variance,
            prior_degrees_of_freedom=prior_degrees_of_freedom,
            residual_degrees_of_freedom=residual_degrees_of_freedom,
            empirical_bayes_method=empirical_bayes_method,
            empirical_bayes_robust=empirical_bayes_robust,
            empirical_bayes_trend=empirical_bayes_trend,
            prior_diagnostics=prior_diagnostics,
            mean_variance_trend_diagnostics=mean_variance_trend_diagnostics,
            policy_provenance=policy_provenance,
            contrast_tables=contrast_tables,
            workflow_provenance=workflow_provenance,
            input_dataset_preprocessing_report=input_dataset_preprocessing_report,
            _assume_owned=True,
        )


def _is_dataset_preprocessing_report(value: object) -> bool:
    """Runtime guard without module-load import cycle with dataset models."""

    try:
        from phospy.science.datasets.models import (
            DatasetPreprocessingReport as _DatasetPreprocessingReport,
        )
    except ImportError:
        return False

    return isinstance(value, _DatasetPreprocessingReport)


def _validate_numeric_matrix(frame: pd.DataFrame, *, field_name: str) -> None:
    require_dataframe(
        frame,
        field_name=field_name,
        allow_empty=False,
        error_type=PhosPyInputError,
    )
    require_non_empty_dataframe(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_unique_index(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_unique_columns(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_finite_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
        allow_missing=False,
    )


def _validate_result_table(
    table: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    _validate_result_table_statistics(
        table=table,
        field_name=field_name,
        allow_imputation_withheld_status=True,
    )
    enforce_phosphosite_identity_contract(
        site_metadata=table,
        field_name=field_name,
        contract=RESULT_TABLE_IDENTITY_CONTRACT,
        error_type=PhosPyInputError,
        compare_raw_site_key_column_before_decode=True,
    )
    enforce_required_identity_text_columns(
        table=table,
        field_name=field_name,
        columns=_PUBLIC_RESULT_IDENTITY_COLUMNS,
        error_type=PhosPyInputError,
    )
    enforce_result_identity_metadata_coherence(
        table=table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )


def _validate_computation_result_table(
    table: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    _validate_result_table_statistics(
        table=table,
        field_name=field_name,
        allow_imputation_withheld_status=False,
    )
    present_identity = [
        column for column in _PUBLIC_RESULT_IDENTITY_COLUMNS if column in table.columns
    ]
    if present_identity:
        joined = ", ".join(present_identity)
        raise PhosPyInputError(
            f"{field_name} must be stat-only and must not include identity columns: "
            f"{joined}"
        )


def _validate_result_table_statistics(
    *,
    table: pd.DataFrame,
    field_name: str,
    allow_imputation_withheld_status: bool,
) -> None:
    require_dataframe(
        table,
        field_name=field_name,
        allow_empty=False,
        error_type=PhosPyInputError,
    )
    require_non_empty_dataframe(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_string_index(
        table.index,
        field_name=f"{field_name}.index",
        error_type=PhosPyInputError,
    )
    require_unique_index(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_unique_columns(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    missing = [
        column for column in _RESULT_STATISTIC_COLUMNS if column not in table.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise PhosPyInputError(f"{field_name} is missing required columns: {joined}")
    stat_table = cast(pd.DataFrame, table[list(_RESULT_STATISTIC_COLUMNS)])
    require_numeric_dataframe(
        stat_table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    if (
        allow_imputation_withheld_status
        and DIFFERENTIAL_RESULT_STATUS_COLUMN in table.columns
    ):
        _validate_imputation_status_statistics(
            table=table,
            stat_table=stat_table,
            field_name=field_name,
        )
    else:
        require_finite_numeric_dataframe(
            stat_table,
            field_name=field_name,
            error_type=PhosPyInputError,
            allow_missing=False,
        )
    _validate_unit_interval_column(
        table=table,
        column_name="P.Value",
        field_name=field_name,
    )
    _validate_unit_interval_column(
        table=table,
        column_name="adj.P.Val",
        field_name=field_name,
    )


def _validate_imputation_status_statistics(
    *,
    table: pd.DataFrame,
    stat_table: pd.DataFrame,
    field_name: str,
) -> None:
    status_values = table[DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str)
    allowed_statuses = {
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        *DIFFERENTIAL_RESULT_WITHHELD_STATUSES,
    }
    unknown_statuses = sorted(set(status_values.tolist()) - allowed_statuses)
    if unknown_statuses:
        raise PhosPyInputError(
            f"{field_name}.{DIFFERENTIAL_RESULT_STATUS_COLUMN} contains unsupported "
            "values: " + ", ".join(repr(value) for value in unknown_statuses)
        )

    status_array = status_values.to_numpy(dtype=str)
    tested_mask: npt.NDArray[np.bool_] = np.asarray(
        status_array == DIFFERENTIAL_RESULT_STATUS_TESTED,
        dtype=bool,
    )
    withheld_mask: npt.NDArray[np.bool_] = np.isin(
        status_array,
        DIFFERENTIAL_RESULT_WITHHELD_STATUSES,
    )
    stat_values_float: npt.NDArray[np.float64] = np.asarray(
        stat_table.to_numpy(dtype=float),
        dtype=np.float64,
    )
    invalid_tested_mask = ~np.isfinite(stat_values_float[tested_mask, :])
    if bool(invalid_tested_mask.any()):
        raise PhosPyInputError(
            f"{field_name} rows with "
            f"{DIFFERENTIAL_RESULT_STATUS_COLUMN}="
            f"{DIFFERENTIAL_RESULT_STATUS_TESTED!r} must contain finite "
            "numeric logFC, t, P.Value, and adj.P.Val values"
        )

    withheld_row_positions = np.flatnonzero(withheld_mask)
    if not int(withheld_row_positions.size):
        return
    withheld_values: npt.NDArray[np.object_] = np.asarray(
        stat_table.to_numpy(dtype=object)[withheld_mask, :],
        dtype=object,
    )
    non_missing_mask: npt.NDArray[np.bool_] = np.asarray(
        ~pd.isna(withheld_values),
        dtype=bool,
    )
    if bool(non_missing_mask.any()):
        invalid_positions = np.argwhere(non_missing_mask)
        previews: list[str] = []
        for row_position, column_position in invalid_positions[:3]:
            source_row_position = int(withheld_row_positions[int(row_position)])
            previews.append(
                f"({stat_table.index[source_row_position]!r}, "
                f"{stat_table.columns[int(column_position)]!r})"
            )
        suffix = (
            ""
            if int(invalid_positions.shape[0]) <= 3
            else f", +{int(invalid_positions.shape[0] - 3)} more"
        )
        raise PhosPyInputError(
            f"{field_name} withheld imputation-policy rows must contain missing "
            "values for logFC, t, P.Value, and adj.P.Val; invalid values: "
            + ", ".join(previews)
            + suffix
        )


def _validate_unit_interval_column(
    *,
    table: pd.DataFrame,
    column_name: str,
    field_name: str,
) -> None:
    column = table[column_name]
    values = column.to_numpy(dtype=float)
    invalid_mask = (values < 0.0) | (values > 1.0)
    if not np.any(invalid_mask):
        return

    invalid_positions = np.flatnonzero(invalid_mask)
    preview: list[str] = []
    for position in invalid_positions[:3]:
        preview.append(f"({table.index[position]!r}, {values[position]:.6g})")
    suffix = (
        ""
        if invalid_positions.size <= 3
        else f", +{int(invalid_positions.size - 3)} more"
    )
    examples = ", ".join(preview)
    raise PhosPyInputError(
        f"{field_name}.{column_name} must be within [0, 1] for each feature; "
        f"invalid values: {examples}{suffix}; "
        f"invalid_entry_count={int(invalid_positions.size)}"
    )
