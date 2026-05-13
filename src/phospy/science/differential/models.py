"""Public models for differential analysis."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

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
    require_unique_columns,
    require_unique_index,
)

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
        if not isinstance(self.trend, bool):
            raise PhosPyInputError("empirical_bayes.trend must be a bool")
        winsor_tail_p = self.winsor_tail_p
        if (
            not isinstance(winsor_tail_p, tuple)
            or len(winsor_tail_p) != 2
            or not all(isinstance(value, int | float) for value in winsor_tail_p)
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
        if not isinstance(design, DesignMatrix):
            raise PhosPyInputError(
                "differential.design must be a DesignMatrix or pandas DataFrame"
            )
        contrasts = self.contrasts
        if isinstance(contrasts, pd.DataFrame):
            contrasts = ContrastMatrix(contrasts)
        if not isinstance(contrasts, ContrastMatrix):
            raise PhosPyInputError(
                "differential.contrasts must be a ContrastMatrix or pandas DataFrame"
            )
        if not isinstance(self.empirical_bayes, EmpiricalBayesConfig):
            raise PhosPyInputError(
                "differential.empirical_bayes must be an EmpiricalBayesConfig"
            )
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "design", design)
        object.__setattr__(self, "contrasts", contrasts)


@dataclass(frozen=True, slots=True, init=False)
class DifferentialAnalysisResult:
    """Differential-analysis output with per-contrast moderated tables."""

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
    workflow_provenance: Mapping[str, object] | None
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
        workflow_provenance: Mapping[str, object] | None = None,
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
        if not contrast_tables:
            raise PhosPyInputError(
                "differential_result.contrast_tables must include at least one contrast"
            )
        if workflow_provenance is not None and not isinstance(
            workflow_provenance,
            Mapping,
        ):
            raise PhosPyInputError(
                "differential_result.workflow_provenance must be a mapping or None"
            )
        owned_tables: dict[str, pd.DataFrame] = {}
        for contrast_name, table in contrast_tables.items():
            if not isinstance(contrast_name, str) or not contrast_name:
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
        object.__setattr__(
            self,
            "workflow_provenance",
            (
                None
                if workflow_provenance is None
                else {str(key): value for key, value in workflow_provenance.items()}
            ),
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
        workflow_provenance: Mapping[str, object] | None = None,
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
            contrast_tables=contrast_tables,
            workflow_provenance=workflow_provenance,
            _assume_owned=True,
        )


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


def _validate_result_table(table: pd.DataFrame, *, field_name: str) -> None:
    _validate_numeric_matrix(table, field_name=field_name)
    required_columns = ("logFC", "t", "P.Value", "adj.P.Val")
    missing = [column for column in required_columns if column not in table.columns]
    if missing:
        joined = ", ".join(missing)
        raise PhosPyInputError(f"{field_name} is missing required columns: {joined}")
