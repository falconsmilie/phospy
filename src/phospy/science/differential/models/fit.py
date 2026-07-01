"""Fit artifact models for differential analysis."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import (
    export_dataframe,
    export_series,
    own_dataframe,
    own_series,
)
from phospy.science.differential.models.diagnostics import (
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.models.tables import (
    validate_computation_result_table_contract,
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
            validate_computation_result_table_contract(
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


__all__ = ["DifferentialComputationResult"]
