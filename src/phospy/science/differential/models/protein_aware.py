"""Private protein-aware differential computation contracts."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import (
    export_dataframe,
    export_series,
    own_dataframe,
    own_series,
)
from phospy.frames.validation import (
    require_columns,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_string_index,
    require_unique_columns,
    require_unique_index,
)
from phospy.science.configs.differential import (
    DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1,
    SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS,
    DifferentialProteinAwareModelMethod,
)
from phospy.science.differential.models.design import ContrastMatrix, DesignMatrix
from phospy.science.differential.models.diagnostics import (
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.models.empirical_bayes_config import (
    EmpiricalBayesConfig,
)
from phospy.science.differential.models.tables import (
    validate_computation_result_table_contract,
)
from phospy.science.statistics.multiple_testing import (
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    SUPPORTED_MULTIPLE_TESTING_CORRECTIONS,
    MultipleTestingCorrection,
)

PROTEIN_AWARE_DIFFERENTIAL_SITE_DIAGNOSTIC_COLUMNS = (
    "site_key",
    "total_protein_row_key",
    "protein_coefficient",
)
PROTEIN_AWARE_DIFFERENTIAL_AUGMENTED_DESIGN_DIAGNOSTIC_COLUMNS = (
    "total_protein_row_key",
    "sample_count",
    "coefficient_count",
    "rank",
    "residual_degrees_of_freedom",
    "condition_number",
    "max_condition_number",
)
PROTEIN_AWARE_DIFFERENTIAL_MATCHED_PAIR_COLUMNS = (
    "site_key",
    "protein_identifier",
    "total_protein_row_key",
)


@dataclass(frozen=True, slots=True, eq=False)
class ProteinAwareDifferentialComputationRequest:
    """Execution input for the protein-covariate-adjusted kernel."""

    __hash__ = object.__hash__

    phosphosite_matrix: pd.DataFrame
    base_design: DesignMatrix | pd.DataFrame
    base_contrasts: ContrastMatrix | pd.DataFrame
    sample_order: Iterable[object]
    matched_pairs: pd.DataFrame
    resolved_protein_covariates: pd.DataFrame
    empirical_bayes: EmpiricalBayesConfig = field(default_factory=EmpiricalBayesConfig)
    multiple_testing_method: MultipleTestingCorrection = (
        MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
    )
    method_id: DifferentialProteinAwareModelMethod = DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1

    def __post_init__(self) -> None:
        self._validate_and_store(assume_owned=False)

    @classmethod
    def _from_owned(
        cls,
        *,
        phosphosite_matrix: pd.DataFrame,
        base_design: DesignMatrix | pd.DataFrame,
        base_contrasts: ContrastMatrix | pd.DataFrame,
        sample_order: Iterable[object],
        matched_pairs: pd.DataFrame,
        resolved_protein_covariates: pd.DataFrame,
        empirical_bayes: EmpiricalBayesConfig | None = None,
        multiple_testing_method: MultipleTestingCorrection = (
            MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
        ),
        method_id: DifferentialProteinAwareModelMethod = (
            DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        ),
    ) -> ProteinAwareDifferentialComputationRequest:
        request = object.__new__(cls)
        object.__setattr__(request, "phosphosite_matrix", phosphosite_matrix)
        object.__setattr__(request, "base_design", base_design)
        object.__setattr__(request, "base_contrasts", base_contrasts)
        object.__setattr__(request, "sample_order", sample_order)
        object.__setattr__(request, "matched_pairs", matched_pairs)
        object.__setattr__(
            request,
            "resolved_protein_covariates",
            resolved_protein_covariates,
        )
        object.__setattr__(
            request,
            "empirical_bayes",
            empirical_bayes if empirical_bayes is not None else EmpiricalBayesConfig(),
        )
        object.__setattr__(
            request,
            "multiple_testing_method",
            multiple_testing_method,
        )
        object.__setattr__(request, "method_id", method_id)
        request._validate_and_store(assume_owned=True)
        return request

    def _validate_and_store(self, *, assume_owned: bool) -> None:
        sample_order = _normalize_label_sequence(
            self.sample_order,
            field_name="protein_aware_differential_request.sample_order",
        )
        sample_index = pd.Index(sample_order)
        phosphosite_matrix = own_dataframe(
            self.phosphosite_matrix,
            field_name="protein_aware_differential_request.phosphosite_matrix",
            error_type=PhosPyInputError,
            assume_owned=assume_owned,
        )
        _validate_finite_numeric_matrix(
            phosphosite_matrix,
            field_name="protein_aware_differential_request.phosphosite_matrix",
        )
        _require_exact_index(
            left=phosphosite_matrix.columns,
            right=sample_index,
            left_name="protein_aware_differential_request.phosphosite_matrix.columns",
            right_name="protein_aware_differential_request.sample_order",
        )

        base_design = _coerce_design_matrix(self.base_design)
        _validate_matrix_labels(
            base_design.frame,
            field_name="protein_aware_differential_request.base_design",
        )
        _require_exact_index(
            left=base_design.frame.index,
            right=sample_index,
            left_name="protein_aware_differential_request.base_design.index",
            right_name="protein_aware_differential_request.sample_order",
        )

        base_contrasts = _coerce_contrast_matrix(self.base_contrasts)
        _validate_matrix_labels(
            base_contrasts.frame,
            field_name="protein_aware_differential_request.base_contrasts",
        )
        _require_exact_index(
            left=base_contrasts.frame.index,
            right=base_design.frame.columns,
            left_name="protein_aware_differential_request.base_contrasts.index",
            right_name="protein_aware_differential_request.base_design.columns",
        )

        matched_pairs = own_dataframe(
            self.matched_pairs,
            field_name="protein_aware_differential_request.matched_pairs",
            error_type=PhosPyInputError,
            assume_owned=assume_owned,
        )
        _validate_matched_pairs(
            matched_pairs,
            expected_site_index=phosphosite_matrix.index,
        )

        resolved_protein_covariates = own_dataframe(
            self.resolved_protein_covariates,
            field_name="protein_aware_differential_request.resolved_protein_covariates",
            error_type=PhosPyInputError,
            assume_owned=assume_owned,
        )
        _validate_finite_numeric_matrix(
            resolved_protein_covariates,
            field_name="protein_aware_differential_request.resolved_protein_covariates",
        )
        _require_exact_index(
            left=resolved_protein_covariates.columns,
            right=sample_index,
            left_name=(
                "protein_aware_differential_request.resolved_protein_covariates.columns"
            ),
            right_name="protein_aware_differential_request.sample_order",
        )
        _require_matched_protein_rows(
            matched_pairs=matched_pairs,
            resolved_protein_covariates=resolved_protein_covariates,
        )
        if not isinstance(cast(object, self.empirical_bayes), EmpiricalBayesConfig):
            raise PhosPyInputError(
                "protein_aware_differential_request.empirical_bayes must be "
                "EmpiricalBayesConfig"
            )

        object.__setattr__(self, "phosphosite_matrix", phosphosite_matrix)
        object.__setattr__(self, "base_design", base_design)
        object.__setattr__(self, "base_contrasts", base_contrasts)
        object.__setattr__(self, "sample_order", sample_order)
        object.__setattr__(self, "matched_pairs", matched_pairs)
        object.__setattr__(
            self,
            "resolved_protein_covariates",
            resolved_protein_covariates,
        )
        object.__setattr__(
            self,
            "multiple_testing_method",
            _normalize_multiple_testing_method(self.multiple_testing_method),
        )
        object.__setattr__(
            self,
            "method_id",
            _normalize_method_id(
                self.method_id,
                field_name="protein_aware_differential_request.method_id",
            ),
        )


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ProteinAwareDifferentialComputationResult:
    """Internal output from protein-aware differential computation."""

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
    protein_coefficient: pd.Series
    site_diagnostics: pd.DataFrame
    augmented_design_diagnostics: pd.DataFrame
    tested_site_ids: tuple[str, ...]
    method_id: DifferentialProteinAwareModelMethod
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
        protein_coefficient: pd.Series,
        site_diagnostics: pd.DataFrame,
        augmented_design_diagnostics: pd.DataFrame,
        tested_site_ids: Iterable[object],
        method_id: DifferentialProteinAwareModelMethod = (
            DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        ),
        _assume_owned: bool = False,
    ) -> None:
        tested_site_ids = _normalize_label_sequence(
            tested_site_ids,
            field_name="protein_aware_differential_result.tested_site_ids",
        )
        tested_index = pd.Index(tested_site_ids, name="site_key")
        residual_variance = _owned_numeric_series(
            residual_variance,
            expected_index=tested_index,
            field_name="protein_aware_differential_result.residual_variance",
            assume_owned=_assume_owned,
        )
        posterior_residual_variance = _owned_numeric_series(
            posterior_residual_variance,
            expected_index=tested_index,
            field_name=(
                "protein_aware_differential_result.posterior_residual_variance"
            ),
            assume_owned=_assume_owned,
        )
        prior_residual_variance = _owned_numeric_series(
            prior_residual_variance,
            expected_index=tested_index,
            field_name="protein_aware_differential_result.prior_residual_variance",
            assume_owned=_assume_owned,
        )
        prior_degrees_of_freedom_series_value = _owned_numeric_series(
            prior_degrees_of_freedom_series_value,
            expected_index=tested_index,
            field_name=(
                "protein_aware_differential_result.prior_degrees_of_freedom_series"
            ),
            assume_owned=_assume_owned,
        )
        protein_coefficient = _owned_numeric_series(
            protein_coefficient,
            expected_index=tested_index,
            field_name="protein_aware_differential_result.protein_coefficient",
            assume_owned=_assume_owned,
        )
        if not isinstance(
            cast(object, prior_diagnostics),
            EmpiricalBayesPriorDiagnostics,
        ):
            raise PhosPyInputError(
                "protein_aware_differential_result.prior_diagnostics must be "
                "EmpiricalBayesPriorDiagnostics"
            )
        _require_exact_index(
            left=prior_diagnostics.prior_variance.index,
            right=tested_index,
            left_name=(
                "protein_aware_differential_result.prior_diagnostics.prior_variance.index"
            ),
            right_name="protein_aware_differential_result.tested_site_ids",
        )
        _require_exact_index(
            left=prior_diagnostics.prior_degrees_of_freedom.index,
            right=tested_index,
            left_name=(
                "protein_aware_differential_result.prior_diagnostics."
                "prior_degrees_of_freedom.index"
            ),
            right_name="protein_aware_differential_result.tested_site_ids",
        )
        if mean_variance_trend_diagnostics is not None:
            if not isinstance(
                cast(object, mean_variance_trend_diagnostics),
                MeanVarianceTrendDiagnostics,
            ):
                raise PhosPyInputError(
                    "protein_aware_differential_result."
                    "mean_variance_trend_diagnostics must be "
                    "MeanVarianceTrendDiagnostics or None"
                )
            _require_exact_index(
                left=mean_variance_trend_diagnostics.mean_intensity.index,
                right=tested_index,
                left_name=(
                    "protein_aware_differential_result."
                    "mean_variance_trend_diagnostics.mean_intensity.index"
                ),
                right_name="protein_aware_differential_result.tested_site_ids",
            )

        residual_dof = _require_positive_finite_float(
            residual_degrees_of_freedom,
            field_name="protein_aware_differential_result.residual_degrees_of_freedom",
        )

        site_diagnostics = own_dataframe(
            site_diagnostics,
            field_name="protein_aware_differential_result.site_diagnostics",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        _validate_site_diagnostics(
            site_diagnostics,
            expected_index=tested_index,
            protein_coefficient=protein_coefficient,
        )
        augmented_design_diagnostics = own_dataframe(
            augmented_design_diagnostics,
            field_name=(
                "protein_aware_differential_result.augmented_design_diagnostics"
            ),
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        _validate_augmented_design_diagnostics(
            augmented_design_diagnostics,
            residual_degrees_of_freedom=residual_dof,
        )
        _require_site_groups_in_augmented_design_diagnostics(
            site_diagnostics=site_diagnostics,
            augmented_design_diagnostics=augmented_design_diagnostics,
        )

        owned_tables: dict[str, pd.DataFrame] = {}
        if not contrast_tables:
            raise PhosPyInputError(
                "protein_aware_differential_result.contrast_tables must include "
                "at least one contrast"
            )
        for contrast_name, table in contrast_tables.items():
            if not isinstance(cast(object, contrast_name), str) or not contrast_name:
                raise PhosPyInputError(
                    "protein_aware_differential_result.contrast_tables keys must "
                    "be non-empty strings"
                )
            owned_table = own_dataframe(
                table,
                field_name=(
                    "protein_aware_differential_result.contrast_tables"
                    f"[{contrast_name!r}]"
                ),
                error_type=PhosPyInputError,
                assume_owned=_assume_owned,
            )
            validate_computation_result_table_contract(
                owned_table,
                field_name=(
                    "protein_aware_differential_result.contrast_tables"
                    f"[{contrast_name!r}]"
                ),
            )
            _require_exact_index(
                left=owned_table.index,
                right=tested_index,
                left_name=(
                    "protein_aware_differential_result.contrast_tables"
                    f"[{contrast_name!r}].index"
                ),
                right_name="protein_aware_differential_result.tested_site_ids",
            )
            owned_tables[contrast_name] = owned_table

        object.__setattr__(self, "residual_variance", residual_variance)
        object.__setattr__(
            self,
            "posterior_residual_variance",
            posterior_residual_variance,
        )
        object.__setattr__(self, "prior_residual_variance", prior_residual_variance)
        object.__setattr__(
            self,
            "prior_degrees_of_freedom_series_value",
            prior_degrees_of_freedom_series_value,
        )
        object.__setattr__(
            self,
            "prior_variance",
            _require_finite_float(
                prior_variance,
                field_name="protein_aware_differential_result.prior_variance",
            ),
        )
        object.__setattr__(
            self,
            "prior_degrees_of_freedom",
            _require_positive_finite_float(
                prior_degrees_of_freedom,
                field_name=(
                    "protein_aware_differential_result.prior_degrees_of_freedom"
                ),
            ),
        )
        object.__setattr__(self, "residual_degrees_of_freedom", residual_dof)
        object.__setattr__(self, "empirical_bayes_method", str(empirical_bayes_method))
        object.__setattr__(self, "empirical_bayes_robust", bool(empirical_bayes_robust))
        object.__setattr__(self, "empirical_bayes_trend", bool(empirical_bayes_trend))
        object.__setattr__(self, "prior_diagnostics", prior_diagnostics)
        object.__setattr__(
            self,
            "mean_variance_trend_diagnostics",
            mean_variance_trend_diagnostics,
        )
        object.__setattr__(self, "protein_coefficient", protein_coefficient)
        object.__setattr__(self, "site_diagnostics", site_diagnostics)
        object.__setattr__(
            self,
            "augmented_design_diagnostics",
            augmented_design_diagnostics,
        )
        object.__setattr__(self, "tested_site_ids", tested_site_ids)
        object.__setattr__(
            self,
            "method_id",
            _normalize_method_id(
                method_id,
                field_name="protein_aware_differential_result.method_id",
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

    def protein_coefficient_series(self) -> pd.Series:
        return export_series(self.protein_coefficient)

    def site_diagnostics_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self.site_diagnostics)

    def augmented_design_diagnostics_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self.augmented_design_diagnostics)


def _coerce_design_matrix(value: object) -> DesignMatrix:
    if isinstance(value, DesignMatrix):
        return DesignMatrix(value.frame)
    if isinstance(value, pd.DataFrame):
        return DesignMatrix(value)
    raise PhosPyInputError(
        "protein_aware_differential_request.base_design must be a DesignMatrix "
        "or pandas DataFrame"
    )


def _coerce_contrast_matrix(value: object) -> ContrastMatrix:
    if isinstance(value, ContrastMatrix):
        return ContrastMatrix(value.frame)
    if isinstance(value, pd.DataFrame):
        return ContrastMatrix(value)
    raise PhosPyInputError(
        "protein_aware_differential_request.base_contrasts must be a ContrastMatrix "
        "or pandas DataFrame"
    )


def _normalize_label_sequence(value: object, *, field_name: str) -> tuple[str, ...]:
    if isinstance(value, str) or not isinstance(value, Iterable):
        raise PhosPyInputError(f"{field_name} must be a sequence of labels")
    raw_labels = tuple(cast(Iterable[object], value))
    if not raw_labels:
        raise PhosPyInputError(f"{field_name} must not be empty")
    invalid_count = sum(
        1
        for label in raw_labels
        if not isinstance(label, str) or not label.strip() or label != label.strip()
    )
    if invalid_count:
        raise PhosPyInputError(
            f"{field_name} must contain stripped non-empty string labels; "
            f"invalid_count={invalid_count}"
        )
    labels = cast(tuple[str, ...], raw_labels)
    duplicate_labels = [
        label for label in dict.fromkeys(labels) if labels.count(label) > 1
    ]
    if duplicate_labels:
        preview = ", ".join(repr(label) for label in duplicate_labels[:5])
        suffix = "" if len(duplicate_labels) <= 5 else " ..."
        raise PhosPyInputError(
            f"{field_name} must contain unique labels; duplicate_labels={preview}{suffix}"
        )
    return labels


def _validate_finite_numeric_matrix(frame: pd.DataFrame, *, field_name: str) -> None:
    require_non_empty_dataframe(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    _validate_matrix_labels(frame, field_name=field_name)
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


def _validate_matrix_labels(frame: pd.DataFrame, *, field_name: str) -> None:
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
    require_string_index(
        frame.index,
        field_name=f"{field_name}.index",
        error_type=PhosPyInputError,
    )
    require_string_index(
        frame.columns,
        field_name=f"{field_name}.columns",
        error_type=PhosPyInputError,
    )


def _validate_matched_pairs(
    matched_pairs: pd.DataFrame,
    *,
    expected_site_index: pd.Index,
) -> None:
    require_unique_columns(
        matched_pairs,
        field_name="protein_aware_differential_request.matched_pairs",
        error_type=PhosPyInputError,
    )
    require_columns(
        matched_pairs,
        field_name="protein_aware_differential_request.matched_pairs",
        required_columns=PROTEIN_AWARE_DIFFERENTIAL_MATCHED_PAIR_COLUMNS,
        error_type=PhosPyInputError,
    )
    site_keys = _require_non_empty_string_column_values(
        matched_pairs,
        column_name="site_key",
        field_name="protein_aware_differential_request.matched_pairs",
    )
    _require_non_empty_string_column_values(
        matched_pairs,
        column_name="protein_identifier",
        field_name="protein_aware_differential_request.matched_pairs",
    )
    _require_non_empty_string_column_values(
        matched_pairs,
        column_name="total_protein_row_key",
        field_name="protein_aware_differential_request.matched_pairs",
    )
    _require_exact_index(
        left=pd.Index(site_keys),
        right=expected_site_index,
        left_name="protein_aware_differential_request.matched_pairs.site_key",
        right_name="protein_aware_differential_request.phosphosite_matrix.index",
    )


def _require_non_empty_string_column_values(
    frame: pd.DataFrame,
    *,
    column_name: str,
    field_name: str,
) -> tuple[str, ...]:
    raw_values = tuple(frame.loc[:, column_name].tolist())
    invalid_count = sum(
        1
        for value in raw_values
        if not isinstance(value, str) or not value.strip() or value != value.strip()
    )
    if invalid_count:
        raise PhosPyInputError(
            f"{field_name}.{column_name} must contain stripped non-empty strings; "
            f"invalid_count={invalid_count}"
        )
    values = cast(tuple[str, ...], raw_values)
    duplicate_values = [
        value for value in dict.fromkeys(values) if values.count(value) > 1
    ]
    if column_name == "site_key" and duplicate_values:
        preview = ", ".join(repr(value) for value in duplicate_values[:5])
        suffix = "" if len(duplicate_values) <= 5 else " ..."
        raise PhosPyInputError(
            f"{field_name}.{column_name} must be unique; duplicate_values={preview}{suffix}"
        )
    return values


def _require_matched_protein_rows(
    *,
    matched_pairs: pd.DataFrame,
    resolved_protein_covariates: pd.DataFrame,
) -> None:
    row_keys = tuple(
        str(value) for value in matched_pairs.loc[:, "total_protein_row_key"].tolist()
    )
    missing = [
        row_key
        for row_key in dict.fromkeys(row_keys)
        if row_key not in resolved_protein_covariates.index
    ]
    if not missing:
        return
    preview = ", ".join(repr(value) for value in missing[:5])
    suffix = "" if len(missing) <= 5 else " ..."
    raise PhosPyInputError(
        "protein_aware_differential_request.resolved_protein_covariates.index must "
        "contain each matched total_protein_row_key; missing="
        f"{preview}{suffix}"
    )


def _normalize_method_id(
    value: object,
    *,
    field_name: str,
) -> DifferentialProteinAwareModelMethod:
    if value not in SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS:
        supported = ", ".join(
            repr(item) for item in SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS
        )
        raise PhosPyInputError(f"{field_name} must be one of: {supported}")
    return value


def _normalize_multiple_testing_method(value: object) -> MultipleTestingCorrection:
    if value not in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS:
        supported = ", ".join(
            repr(item) for item in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
        )
        raise PhosPyInputError(
            "protein_aware_differential_request.multiple_testing_method must be "
            f"one of: {supported}"
        )
    return value


def _owned_numeric_series(
    value: pd.Series,
    *,
    expected_index: pd.Index,
    field_name: str,
    assume_owned: bool,
) -> pd.Series:
    series = own_series(
        value,
        field_name=field_name,
        error_type=PhosPyInputError,
        assume_owned=assume_owned,
    )
    _require_exact_index(
        left=series.index,
        right=expected_index,
        left_name=f"{field_name}.index",
        right_name="protein_aware_differential_result.tested_site_ids",
    )
    if pd.api.types.is_bool_dtype(series) or not pd.api.types.is_numeric_dtype(series):
        raise PhosPyInputError(f"{field_name} must contain numeric values")
    values = np.asarray(series.to_numpy(dtype=float), dtype=float)
    if not np.isfinite(values).all():
        raise PhosPyInputError(f"{field_name} must contain finite numeric values")
    return series


def _validate_site_diagnostics(
    frame: pd.DataFrame,
    *,
    expected_index: pd.Index,
    protein_coefficient: pd.Series,
) -> None:
    require_columns(
        frame,
        field_name="protein_aware_differential_result.site_diagnostics",
        required_columns=PROTEIN_AWARE_DIFFERENTIAL_SITE_DIAGNOSTIC_COLUMNS,
        error_type=PhosPyInputError,
    )
    require_unique_columns(
        frame,
        field_name="protein_aware_differential_result.site_diagnostics",
        error_type=PhosPyInputError,
    )
    _require_exact_index(
        left=frame.index,
        right=expected_index,
        left_name="protein_aware_differential_result.site_diagnostics.index",
        right_name="protein_aware_differential_result.tested_site_ids",
    )
    site_key_values = tuple(str(value) for value in frame.loc[:, "site_key"].tolist())
    _require_exact_index(
        left=pd.Index(site_key_values),
        right=expected_index,
        left_name="protein_aware_differential_result.site_diagnostics.site_key",
        right_name="protein_aware_differential_result.tested_site_ids",
    )
    _require_non_empty_string_column_values(
        frame,
        column_name="total_protein_row_key",
        field_name="protein_aware_differential_result.site_diagnostics",
    )
    coefficients = pd.to_numeric(frame.loc[:, "protein_coefficient"], errors="coerce")
    values = np.asarray(coefficients.to_numpy(dtype=float), dtype=float)
    if not np.isfinite(values).all():
        raise PhosPyInputError(
            "protein_aware_differential_result.site_diagnostics.protein_coefficient "
            "must contain finite numeric values"
        )
    expected = protein_coefficient.to_numpy(dtype=float)
    if not np.array_equal(values, expected):
        raise PhosPyInputError(
            "protein_aware_differential_result.site_diagnostics.protein_coefficient "
            "must match protein_aware_differential_result.protein_coefficient"
        )


def _validate_augmented_design_diagnostics(
    frame: pd.DataFrame,
    *,
    residual_degrees_of_freedom: float,
) -> None:
    require_columns(
        frame,
        field_name="protein_aware_differential_result.augmented_design_diagnostics",
        required_columns=PROTEIN_AWARE_DIFFERENTIAL_AUGMENTED_DESIGN_DIAGNOSTIC_COLUMNS,
        error_type=PhosPyInputError,
    )
    require_non_empty_dataframe(
        frame,
        field_name="protein_aware_differential_result.augmented_design_diagnostics",
        error_type=PhosPyInputError,
    )
    require_unique_index(
        frame,
        field_name="protein_aware_differential_result.augmented_design_diagnostics",
        error_type=PhosPyInputError,
    )
    require_unique_columns(
        frame,
        field_name="protein_aware_differential_result.augmented_design_diagnostics",
        error_type=PhosPyInputError,
    )
    require_string_index(
        frame.index,
        field_name=(
            "protein_aware_differential_result.augmented_design_diagnostics.index"
        ),
        error_type=PhosPyInputError,
    )
    row_keys = tuple(str(value) for value in frame["total_protein_row_key"].tolist())
    _require_exact_index(
        left=pd.Index(row_keys),
        right=frame.index,
        left_name=(
            "protein_aware_differential_result.augmented_design_diagnostics."
            "total_protein_row_key"
        ),
        right_name=(
            "protein_aware_differential_result.augmented_design_diagnostics.index"
        ),
    )
    numeric_columns = [
        column
        for column in PROTEIN_AWARE_DIFFERENTIAL_AUGMENTED_DESIGN_DIAGNOSTIC_COLUMNS
        if column != "total_protein_row_key"
    ]
    numeric = frame.loc[:, numeric_columns]
    require_numeric_dataframe(
        numeric,
        field_name=("protein_aware_differential_result.augmented_design_diagnostics"),
        error_type=PhosPyInputError,
    )
    require_finite_numeric_dataframe(
        numeric,
        field_name=("protein_aware_differential_result.augmented_design_diagnostics"),
        error_type=PhosPyInputError,
        allow_missing=False,
    )
    group_dof = numeric.loc[:, "residual_degrees_of_freedom"].to_numpy(dtype=float)
    if not np.allclose(
        group_dof,
        float(residual_degrees_of_freedom),
        rtol=0.0,
        atol=1.0e-12,
    ):
        raise PhosPyInputError(
            "protein_aware_differential_result.augmented_design_diagnostics."
            "residual_degrees_of_freedom must match the common tested residual "
            "degrees of freedom"
        )


def _require_site_groups_in_augmented_design_diagnostics(
    *,
    site_diagnostics: pd.DataFrame,
    augmented_design_diagnostics: pd.DataFrame,
) -> None:
    group_index = augmented_design_diagnostics.index
    missing = [
        row_key
        for row_key in dict.fromkeys(
            str(value)
            for value in site_diagnostics.loc[:, "total_protein_row_key"].tolist()
        )
        if row_key not in group_index
    ]
    if not missing:
        return
    preview = ", ".join(repr(value) for value in missing[:5])
    suffix = "" if len(missing) <= 5 else " ..."
    raise PhosPyInputError(
        "protein_aware_differential_result.augmented_design_diagnostics.index must "
        "contain each per-site total_protein_row_key; missing="
        f"{preview}{suffix}"
    )


def _require_exact_index(
    *,
    left: pd.Index,
    right: pd.Index,
    left_name: str,
    right_name: str,
) -> None:
    require_exact_index_match(
        left=left,
        right=right,
        left_name=left_name,
        right_name=right_name,
        error_type=PhosPyInputError,
    )


def _require_finite_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric")
    number = float(value)
    if not math.isfinite(number):
        raise PhosPyInputError(f"{field_name} must be finite")
    return number


def _require_positive_finite_float(value: object, *, field_name: str) -> float:
    number = _require_finite_float(value, field_name=field_name)
    if number <= 0.0:
        raise PhosPyInputError(f"{field_name} must be > 0.0")
    return number


__all__ = [
    "PROTEIN_AWARE_DIFFERENTIAL_AUGMENTED_DESIGN_DIAGNOSTIC_COLUMNS",
    "PROTEIN_AWARE_DIFFERENTIAL_MATCHED_PAIR_COLUMNS",
    "PROTEIN_AWARE_DIFFERENTIAL_SITE_DIAGNOSTIC_COLUMNS",
    "ProteinAwareDifferentialComputationRequest",
    "ProteinAwareDifferentialComputationResult",
]
