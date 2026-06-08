"""Public models for differential analysis."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import (
    export_dataframe,
    export_series,
    own_dataframe,
    own_series,
)
from phospy.science.sites.validation import require_site_key_index
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_string_index,
    require_unique_columns,
    require_unique_index,
)

if TYPE_CHECKING:
    from phospy.science.datasets.models import DatasetPreprocessingReport

EMPIRICAL_BAYES_METHOD_STANDARD = "standard"
EMPIRICAL_BAYES_METHOD_ROBUST = "robust"
SUPPORTED_EMPIRICAL_BAYES_METHODS: tuple[str, ...] = (
    EMPIRICAL_BAYES_METHOD_STANDARD,
    EMPIRICAL_BAYES_METHOD_ROBUST,
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
class DifferentialDesignMatrixSummary:
    """Structured summary of the resolved differential design matrix."""

    formula: str
    sample_labels: tuple[str, ...]
    coefficient_labels: tuple[str, ...]
    sample_count: int
    coefficient_count: int
    rank: int
    residual_degrees_of_freedom: float

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


@dataclass(frozen=True, slots=True)
class DifferentialContrastDefinition:
    """Structured differential contrast definition."""

    name: str
    numerator_condition: str
    denominator_condition: str
    coefficients: tuple[tuple[str, float], ...]

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

    def __post_init__(self) -> None:
        if not self.policy:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values.policy must be non-empty"
            )
        if not self.stage:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values.stage must be non-empty"
            )


@dataclass(frozen=True, slots=True)
class DifferentialUnsupportedDesignPolicyProvenance:
    """Structured record of unsupported differential-design features."""

    intentionally_rejected_features: tuple[str, ...]
    enforcement_stage: str

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
        _require_identity_columns: bool = True,
        _assume_owned: bool = False,
    ) -> None:
        if not _require_identity_columns and not _assume_owned:
            raise PhosPyInputError(
                "differential_result._require_identity_columns=False is reserved for "
                "internal stat-only compatibility construction"
            )
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
                require_identity_columns=_require_identity_columns,
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
            contrast_name: _export_public_contrast_table(table)
            for contrast_name, table in self._contrast_tables.items()
        }

    def table_for(self, contrast_name: str) -> pd.DataFrame:
        if contrast_name not in self._contrast_tables:
            available = ", ".join(sorted(self._contrast_tables))
            raise KeyError(
                f"unknown contrast {contrast_name!r}; available: {available}"
            )
        table = self._contrast_tables[contrast_name]
        return _export_public_contrast_table(table)

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
        require_identity_columns: bool = True,
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
            _require_identity_columns=require_identity_columns,
            _assume_owned=True,
        )

    @classmethod
    def _from_owned_stat_only_tables(
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
        """Internal compatibility path for computation-only result tables.

        This deliberately allows stat-only contrast tables without
        ``site_key``/``display_id``/``gene_symbol``/``site`` metadata. Workflow
        result assembly must use ``_from_owned`` or the public constructor so the
        public identity contract stays strict.
        """

        return cls._from_owned(
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
            require_identity_columns=False,
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
    require_identity_columns: bool,
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
    required_stat_columns = ("logFC", "t", "P.Value", "adj.P.Val")
    missing = [
        column for column in required_stat_columns if column not in table.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise PhosPyInputError(f"{field_name} is missing required columns: {joined}")
    stat_table = cast(pd.DataFrame, table[list(required_stat_columns)])
    require_numeric_dataframe(
        stat_table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_finite_numeric_dataframe(
        stat_table,
        field_name=field_name,
        error_type=PhosPyInputError,
        allow_missing=False,
    )
    enforce_identity_columns = require_identity_columns or (
        "site_key" in table.columns or "display_id" in table.columns
    )
    if enforce_identity_columns:
        identity_required = ("site_key", "display_id", "gene_symbol", "site")
        missing_identity = [
            column for column in identity_required if column not in table.columns
        ]
        if missing_identity:
            joined = ", ".join(missing_identity)
            raise PhosPyInputError(
                f"{field_name} is missing required columns: {joined}"
            )
        for column_name in identity_required:
            require_non_empty_string_column(
                table,
                field_name=field_name,
                column_name=column_name,
                error_type=PhosPyInputError,
            )
        _validate_site_key_column_matches_index(table=table, field_name=field_name)
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


def _validate_site_key_column_matches_index(
    *,
    table: pd.DataFrame,
    field_name: str,
) -> None:
    require_site_key_index(
        table.index,
        field_name=f"{field_name}.index",
        error_type=PhosPyInputError,
    )
    site_key_column = table["site_key"]
    site_key_values = [str(value) for value in site_key_column.tolist()]
    index_values = [str(value) for value in table.index.tolist()]
    mismatches = [
        idx
        for idx, site_key in zip(index_values, site_key_values, strict=True)
        if idx != site_key
    ]
    if not mismatches:
        return
    preview = ", ".join(repr(value) for value in mismatches[:5])
    suffix = "" if len(mismatches) <= 5 else " ..."
    raise PhosPyInputError(
        f"{field_name}.site_key must exactly match {field_name}.index; "
        f"mismatched_labels={preview}{suffix}"
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


def _export_public_contrast_table(table: pd.DataFrame) -> pd.DataFrame:
    exported = export_dataframe(table)
    if {"site_key", "display_id", "gene_symbol", "site"}.issubset(exported.columns):
        return exported
    if _index_uses_site_key_identity(exported.index):
        return exported
    return cast(pd.DataFrame, exported[["logFC", "t", "P.Value", "adj.P.Val"]])


def _index_uses_site_key_identity(index: pd.Index) -> bool:
    try:
        require_site_key_index(
            index,
            field_name="differential_result.contrast_table.index",
            error_type=PhosPyInputError,
        )
        return True
    except PhosPyInputError:
        return False
