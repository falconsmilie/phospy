"""Public result DTOs for differential analysis."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING, cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.comparison import (
    dataframe_mapping_equals,
    optional_dataframe_equals,
    series_equals,
)
from phospy.frames.ownership import (
    export_dataframe,
    export_series,
    own_dataframe,
    own_series,
)
from phospy.provenance.immutability import (
    freeze_optional_json_mapping,
    thaw_json_mapping,
)
from phospy.science.differential.models.diagnostics import (
    DifferentialModelDiagnostics,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.models.provenance import DifferentialPolicyProvenance
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_WITHHELD_STATUSES,
    validate_result_table_contract,
)
from phospy.science.result_caveats import ResultCaveat, validate_result_caveats

if TYPE_CHECKING:
    from phospy.science.datasets.models import DatasetPreprocessingReport


@dataclass(frozen=True, slots=True, init=False, eq=False)
class DifferentialAnalysisResult:
    """Differential-analysis output with per-contrast moderated tables.

    Public contrast tables are indexed by protein-scoped ``site_key`` values and
    must include ``site_key``, ``display_id``, ``gene_symbol``, and ``site``
    columns.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit differential-result content
    comparison.
    """

    __hash__ = object.__hash__

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
    diagnostics: DifferentialModelDiagnostics
    policy_provenance: DifferentialPolicyProvenance | None
    workflow_provenance: Mapping[str, object] | None
    caveats: tuple[ResultCaveat, ...]
    input_dataset_preprocessing_report: DatasetPreprocessingReport | None
    _feature_eligibility: pd.DataFrame | None
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
        diagnostics: DifferentialModelDiagnostics | None = None,
        policy_provenance: DifferentialPolicyProvenance | None = None,
        workflow_provenance: Mapping[str, object] | None = None,
        caveats: tuple[ResultCaveat, ...] = (),
        input_dataset_preprocessing_report: DatasetPreprocessingReport | None = None,
        feature_eligibility: pd.DataFrame | None = None,
    ) -> None:
        self._init_differential_result(
            residual_variance=residual_variance,
            posterior_residual_variance=posterior_residual_variance,
            prior_residual_variance=prior_residual_variance,
            prior_degrees_of_freedom_series_value=(
                prior_degrees_of_freedom_series_value
            ),
            prior_variance=prior_variance,
            prior_degrees_of_freedom=prior_degrees_of_freedom,
            residual_degrees_of_freedom=residual_degrees_of_freedom,
            empirical_bayes_method=empirical_bayes_method,
            empirical_bayes_robust=empirical_bayes_robust,
            empirical_bayes_trend=empirical_bayes_trend,
            prior_diagnostics=prior_diagnostics,
            mean_variance_trend_diagnostics=mean_variance_trend_diagnostics,
            contrast_tables=contrast_tables,
            diagnostics=diagnostics,
            policy_provenance=policy_provenance,
            workflow_provenance=workflow_provenance,
            caveats=caveats,
            input_dataset_preprocessing_report=input_dataset_preprocessing_report,
            feature_eligibility=feature_eligibility,
            assume_owned=False,
        )

    def _init_differential_result(
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
        diagnostics: DifferentialModelDiagnostics | None = None,
        policy_provenance: DifferentialPolicyProvenance | None = None,
        workflow_provenance: Mapping[str, object] | None = None,
        caveats: tuple[ResultCaveat, ...] = (),
        input_dataset_preprocessing_report: DatasetPreprocessingReport | None = None,
        feature_eligibility: pd.DataFrame | None = None,
        assume_owned: bool,
    ) -> None:
        residual_variance = own_series(
            residual_variance,
            field_name="differential_result.residual_variance",
            error_type=PhosPyInputError,
            assume_owned=assume_owned,
        )
        posterior_residual_variance = own_series(
            posterior_residual_variance,
            field_name="differential_result.posterior_residual_variance",
            error_type=PhosPyInputError,
            assume_owned=assume_owned,
        )
        prior_residual_variance = own_series(
            prior_residual_variance,
            field_name="differential_result.prior_residual_variance",
            error_type=PhosPyInputError,
            assume_owned=assume_owned,
        )
        prior_degrees_of_freedom_series_value = own_series(
            prior_degrees_of_freedom_series_value,
            field_name="differential_result.prior_degrees_of_freedom_series",
            error_type=PhosPyInputError,
            assume_owned=assume_owned,
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
        if diagnostics is None:
            diagnostics = _build_default_model_diagnostics(
                residual_variance=residual_variance,
                residual_degrees_of_freedom=residual_degrees_of_freedom,
                empirical_bayes_method=empirical_bayes_method,
                empirical_bayes_robust=empirical_bayes_robust,
                empirical_bayes_trend=empirical_bayes_trend,
                policy_provenance=policy_provenance,
            )
        if not isinstance(cast(object, diagnostics), DifferentialModelDiagnostics):
            raise PhosPyInputError(
                "differential_result.diagnostics must be DifferentialModelDiagnostics"
            )
        frozen_workflow_provenance = freeze_optional_json_mapping(
            workflow_provenance,
            field_name="differential_result.workflow_provenance",
        )
        caveats = validate_result_caveats(
            caveats,
            field_name="differential_result.caveats",
        )
        if (
            input_dataset_preprocessing_report is not None
            and not _is_dataset_preprocessing_report(input_dataset_preprocessing_report)
        ):
            raise PhosPyInputError(
                "differential_result.input_dataset_preprocessing_report must be "
                "DatasetPreprocessingReport or None"
            )
        owned_feature_eligibility: pd.DataFrame | None = None
        if feature_eligibility is not None:
            owned_feature_eligibility = own_dataframe(
                feature_eligibility,
                field_name="differential_result.feature_eligibility",
                error_type=PhosPyInputError,
                assume_owned=assume_owned,
            )
            _validate_feature_eligibility_table(
                owned_feature_eligibility,
                expected_index=residual_variance.index,
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
                assume_owned=assume_owned,
            )
            validate_result_table_contract(
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
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "policy_provenance", policy_provenance)
        object.__setattr__(
            self,
            "workflow_provenance",
            frozen_workflow_provenance,
        )
        object.__setattr__(self, "caveats", caveats)
        object.__setattr__(
            self,
            "input_dataset_preprocessing_report",
            input_dataset_preprocessing_report,
        )
        object.__setattr__(
            self,
            "_feature_eligibility",
            owned_feature_eligibility,
        )
        object.__setattr__(self, "_contrast_tables", owned_tables)

    @property
    def contrast_tables(self) -> dict[str, pd.DataFrame]:
        return {
            contrast_name: export_dataframe(table)
            for contrast_name, table in self._contrast_tables.items()
        }

    @property
    def feature_eligibility(self) -> pd.DataFrame | None:
        if self._feature_eligibility is None:
            return None
        return export_dataframe(self._feature_eligibility)

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

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible differential result payload."""

        return {
            "caveats": [caveat.to_payload() for caveat in self.caveats],
            "diagnostics": self.diagnostics.to_payload(),
            "workflow_provenance": (
                None
                if self.workflow_provenance is None
                else thaw_json_mapping(
                    self.workflow_provenance,
                    field_name="differential_result.workflow_provenance",
                )
            ),
            "policy_provenance": (
                None
                if self.policy_provenance is None
                else _json_payload(asdict(self.policy_provenance))
            ),
            "empirical_bayes": {
                "method": self.empirical_bayes_method,
                "robust": self.empirical_bayes_robust,
                "trend": self.empirical_bayes_trend,
                "prior_variance": _json_scalar(self.prior_variance),
                "prior_degrees_of_freedom": _json_scalar(self.prior_degrees_of_freedom),
                "residual_degrees_of_freedom": _json_scalar(
                    self.residual_degrees_of_freedom
                ),
            },
            "contrast_tables": {
                contrast_name: _dataframe_records_payload(table)
                for contrast_name, table in self._contrast_tables.items()
            },
            **(
                {}
                if self._feature_eligibility is None
                else {
                    "feature_eligibility": _dataframe_records_payload(
                        self._feature_eligibility
                    )
                }
            ),
        }

    def scientifically_equals(
        self,
        other: object,
        *,
        include_provenance: bool = True,
    ) -> bool:
        """Return ``True`` when another differential result has the same content."""

        if not isinstance(other, DifferentialAnalysisResult):
            return False
        same_content = (
            series_equals(self.residual_variance, other.residual_variance)
            and series_equals(
                self.posterior_residual_variance,
                other.posterior_residual_variance,
            )
            and series_equals(
                self.prior_residual_variance,
                other.prior_residual_variance,
            )
            and series_equals(
                self.prior_degrees_of_freedom_series_value,
                other.prior_degrees_of_freedom_series_value,
            )
            and self.prior_variance == other.prior_variance
            and self.prior_degrees_of_freedom == other.prior_degrees_of_freedom
            and self.residual_degrees_of_freedom == other.residual_degrees_of_freedom
            and self.empirical_bayes_method == other.empirical_bayes_method
            and self.empirical_bayes_robust == other.empirical_bayes_robust
            and self.empirical_bayes_trend == other.empirical_bayes_trend
            and self.prior_diagnostics.scientifically_equals(other.prior_diagnostics)
            and _optional_trend_diagnostics_equals(
                self.mean_variance_trend_diagnostics,
                other.mean_variance_trend_diagnostics,
            )
            and self.diagnostics == other.diagnostics
            and self.policy_provenance == other.policy_provenance
            and self.caveats == other.caveats
            and _optional_preprocessing_report_equals(
                self.input_dataset_preprocessing_report,
                other.input_dataset_preprocessing_report,
            )
            and optional_dataframe_equals(
                self._feature_eligibility,
                other._feature_eligibility,
            )
            and dataframe_mapping_equals(
                self._contrast_tables,
                other._contrast_tables,
            )
        )
        if not same_content:
            return False
        if include_provenance and self.workflow_provenance != other.workflow_provenance:
            return False
        return True

    @classmethod
    def from_trusted_owned(
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
        diagnostics: DifferentialModelDiagnostics | None = None,
        policy_provenance: DifferentialPolicyProvenance | None = None,
        workflow_provenance: Mapping[str, object] | None = None,
        caveats: tuple[ResultCaveat, ...] = (),
        input_dataset_preprocessing_report: DatasetPreprocessingReport | None = None,
        feature_eligibility: pd.DataFrame | None = None,
    ) -> DifferentialAnalysisResult:
        """Construct from already-owned tables at trusted workflow boundaries."""

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
            contrast_tables=contrast_tables,
            diagnostics=diagnostics,
            policy_provenance=policy_provenance,
            workflow_provenance=workflow_provenance,
            caveats=caveats,
            input_dataset_preprocessing_report=input_dataset_preprocessing_report,
            feature_eligibility=feature_eligibility,
        )

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
        diagnostics: DifferentialModelDiagnostics | None = None,
        policy_provenance: DifferentialPolicyProvenance | None = None,
        workflow_provenance: Mapping[str, object] | None = None,
        caveats: tuple[ResultCaveat, ...] = (),
        input_dataset_preprocessing_report: DatasetPreprocessingReport | None = None,
        feature_eligibility: pd.DataFrame | None = None,
    ) -> DifferentialAnalysisResult:
        result = object.__new__(cls)
        DifferentialAnalysisResult._init_differential_result(
            result,
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
            diagnostics=diagnostics,
            policy_provenance=policy_provenance,
            contrast_tables=contrast_tables,
            workflow_provenance=workflow_provenance,
            caveats=caveats,
            input_dataset_preprocessing_report=input_dataset_preprocessing_report,
            feature_eligibility=feature_eligibility,
            assume_owned=True,
        )
        return result


def _is_dataset_preprocessing_report(value: object) -> bool:
    """Runtime guard without module-load import cycle with dataset models."""

    try:
        from phospy.science.datasets.models import (
            DatasetPreprocessingReport as _DatasetPreprocessingReport,
        )
    except ImportError:
        return False

    return isinstance(value, _DatasetPreprocessingReport)


def _optional_trend_diagnostics_equals(
    left: MeanVarianceTrendDiagnostics | None,
    right: MeanVarianceTrendDiagnostics | None,
) -> bool:
    if left is None or right is None:
        return left is right
    return left.scientifically_equals(right)


def _optional_preprocessing_report_equals(
    left: object | None,
    right: object | None,
) -> bool:
    if left is None or right is None:
        return left is right
    method = getattr(left, "scientifically_equals", None)
    if callable(method):
        return bool(method(right))
    return left == right


def _validate_feature_eligibility_table(
    table: pd.DataFrame,
    *,
    expected_index: pd.Index,
) -> None:
    if not table.index.equals(expected_index):
        raise PhosPyInputError(
            "differential_result.feature_eligibility index must match matrix feature "
            "index"
        )
    required_columns = (
        "site_key",
        DIFFERENTIAL_RESULT_STATUS_COLUMN,
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    )
    missing = [column for column in required_columns if column not in table.columns]
    if missing:
        raise PhosPyInputError(
            "differential_result.feature_eligibility is missing required columns: "
            + ", ".join(missing)
        )
    site_key_column = table["site_key"]
    site_key_values = [str(value) for value in site_key_column.tolist()]
    expected_site_key_values = [str(value) for value in expected_index.tolist()]
    if site_key_values != expected_site_key_values:
        raise PhosPyInputError(
            "differential_result.feature_eligibility.site_key must exactly match "
            "the feature index"
        )
    status_column = table[DIFFERENTIAL_RESULT_STATUS_COLUMN]
    status_values = np.asarray(
        [str(value) for value in status_column.tolist()],
        dtype=str,
    )
    allowed_statuses = {
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        *DIFFERENTIAL_RESULT_WITHHELD_STATUSES,
    }
    unknown_statuses = sorted(set(status_values.tolist()) - allowed_statuses)
    if unknown_statuses:
        raise PhosPyInputError(
            "differential_result.feature_eligibility.result_status contains "
            "unsupported values: "
            + ", ".join(repr(value) for value in unknown_statuses)
        )
    withheld_mask = np.isin(status_values, DIFFERENTIAL_RESULT_WITHHELD_STATUSES)
    reason_column = table[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
    empty_reason_mask = np.asarray(
        [str(value).strip() == "" for value in reason_column.tolist()],
        dtype=bool,
    )
    invalid_positions = np.flatnonzero(withheld_mask & empty_reason_mask)
    if int(invalid_positions.size):
        invalid_labels = [
            str(table.index[int(position)]) for position in invalid_positions[:3]
        ]
        raise PhosPyInputError(
            "differential_result.feature_eligibility withheld rows must include "
            "non-empty result_status_reason values: " + ", ".join(invalid_labels[:3])
        )


def _build_default_model_diagnostics(
    *,
    residual_variance: pd.Series,
    residual_degrees_of_freedom: float,
    empirical_bayes_method: str,
    empirical_bayes_robust: bool,
    empirical_bayes_trend: bool,
    policy_provenance: DifferentialPolicyProvenance | None,
) -> DifferentialModelDiagnostics:
    if policy_provenance is not None:
        design = policy_provenance.design
        covariate_terms = tuple(
            column for covariate in design.covariates for column in covariate.columns
        )
        batch_or_covariate_terms = (*covariate_terms, *design.block_column_names)
        unsupported_assumptions = (
            *policy_provenance.unsupported_design.intentionally_rejected_features,
            *design.limitations,
        )
        return DifferentialModelDiagnostics(
            model_type="moderated_ols_fixed_effect",
            design_columns=design.coefficient_labels,
            contrast_definitions=policy_provenance.contrasts,
            rank=design.rank,
            n_samples=design.sample_count,
            n_sites=int(residual_variance.size),
            residual_degrees_of_freedom=float(residual_degrees_of_freedom),
            decomposition_method=design.decomposition_method,
            solver=design.solver,
            column_scale_method=design.column_scale_method,
            rank_tolerance_policy=design.rank_tolerance_policy,
            rank_tolerance=design.rank_tolerance,
            condition_number=design.condition_number,
            max_condition_number=design.max_condition_number,
            singular_values=design.singular_values,
            variance_method="ordinary_least_squares_residual_variance",
            moderation_method=_moderation_method(
                empirical_bayes_method,
                robust=empirical_bayes_robust,
                trend=empirical_bayes_trend,
            ),
            multiple_testing_method=(
                policy_provenance.statistical_testing.adjusted_p_value_method
            ),
            imputation_policy=policy_provenance.missing_values.imputed_value_policy,
            missing_value_policy=policy_provenance.missing_values.policy,
            intensity_scale=(
                policy_provenance.statistical_testing.input_intensity_scale
            ),
            normalisation_state="not_recorded",
            batch_or_covariate_terms=batch_or_covariate_terms,
            unsupported_assumptions=_unique_text(unsupported_assumptions),
            warnings=(),
        )

    warning = (
        "DifferentialAnalysisResult was constructed without workflow policy "
        "provenance; design scope diagnostics are not recorded."
    )
    return DifferentialModelDiagnostics(
        model_type="not_recorded",
        design_columns=(),
        contrast_definitions=(),
        rank=0,
        n_samples=0,
        n_sites=int(residual_variance.size),
        residual_degrees_of_freedom=float(residual_degrees_of_freedom),
        variance_method="not_recorded",
        moderation_method=_moderation_method(
            empirical_bayes_method,
            robust=empirical_bayes_robust,
            trend=empirical_bayes_trend,
        ),
        multiple_testing_method="not_recorded",
        imputation_policy="not_recorded",
        missing_value_policy="not_recorded",
        intensity_scale="not_recorded",
        normalisation_state="not_recorded",
        batch_or_covariate_terms=(),
        unsupported_assumptions=(warning,),
        warnings=(warning,),
    )


def _moderation_method(method: str, *, robust: bool, trend: bool) -> str:
    parts = ["empirical_bayes", str(method)]
    if robust and str(method) != "robust":
        parts.append("robust")
    if trend:
        parts.append("trend")
    return "_".join(parts)


def _unique_text(values: tuple[str, ...]) -> tuple[str, ...]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return tuple(ordered)


def _dataframe_records_payload(frame: pd.DataFrame) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    columns = [str(column) for column in frame.columns]
    values = cast(npt.NDArray[np.object_], frame.to_numpy(dtype=object))
    for row_position in range(int(values.shape[0])):
        record: dict[str, object] = {}
        for column_position, column in enumerate(columns):
            record[column] = _json_scalar(
                cast(object, values[row_position, column_position])
            )
        records.append(record)
    return records


def _json_scalar(value: object) -> object:
    if value is None:
        return None
    item: object = (
        cast(object, value.item()) if isinstance(value, np.generic) else value
    )
    if item is pd.NA or item is pd.NaT:
        return None
    if isinstance(item, float) and not math.isfinite(item):
        return None
    return item


def _json_payload(value: object) -> object:
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, object], value)
        return {str(key): _json_payload(item) for key, item in mapping.items()}
    if isinstance(value, tuple):
        values = cast(tuple[object, ...], value)
        return [_json_payload(item) for item in values]
    if isinstance(value, list):
        values = cast(list[object], value)
        return [_json_payload(item) for item in values]
    return _json_scalar(value)


__all__ = ["DifferentialAnalysisResult"]
