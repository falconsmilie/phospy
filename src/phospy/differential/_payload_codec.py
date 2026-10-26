"""JSON-compatible payload encoding for differential-analysis results."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import asdict
from typing import Protocol, cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.differential.diagnostics import (
    DifferentialModelDiagnostics,
    MeanVarianceTrendDiagnostics,
    ProteinAwareDifferentialDiagnostics,
    QuantificationDepthTrendDiagnostics,
)
from phospy.differential.provenance.models import DifferentialPolicyProvenance
from phospy.provenance.immutability import thaw_json_mapping
from phospy.provenance.result_caveats import ResultCaveat


class _DifferentialResultPayloadSource(Protocol):
    @property
    def caveats(self) -> tuple[ResultCaveat, ...]: ...

    @property
    def diagnostics(self) -> DifferentialModelDiagnostics: ...

    @property
    def protein_aware_diagnostics(
        self,
    ) -> ProteinAwareDifferentialDiagnostics | None: ...

    @property
    def workflow_provenance(self) -> Mapping[str, object] | None: ...

    @property
    def policy_provenance(self) -> DifferentialPolicyProvenance | None: ...

    @property
    def empirical_bayes_method(self) -> str: ...

    @property
    def empirical_bayes_robust(self) -> bool: ...

    @property
    def empirical_bayes_trend(self) -> bool: ...

    @property
    def prior_variance(self) -> float: ...

    @property
    def prior_degrees_of_freedom(self) -> float: ...

    @property
    def residual_degrees_of_freedom(self) -> float: ...

    @property
    def prior_residual_variance(self) -> pd.Series: ...

    @property
    def prior_degrees_of_freedom_series_value(self) -> pd.Series: ...

    @property
    def mean_variance_trend_diagnostics(
        self,
    ) -> MeanVarianceTrendDiagnostics | None: ...

    @property
    def quantification_depth_trend_diagnostics(
        self,
    ) -> QuantificationDepthTrendDiagnostics | None: ...


def encode_differential_result_payload(
    result: _DifferentialResultPayloadSource,
    *,
    contrast_tables: Mapping[str, pd.DataFrame],
    feature_eligibility: pd.DataFrame | None,
) -> dict[str, object]:
    """Encode owned Differential result state without exporting frame copies."""

    return {
        "caveats": [caveat.to_payload() for caveat in result.caveats],
        "diagnostics": result.diagnostics.to_payload(),
        **(
            {}
            if result.protein_aware_diagnostics is None
            else {
                "protein_aware_diagnostics": (
                    result.protein_aware_diagnostics.to_payload()
                )
            }
        ),
        "workflow_provenance": (
            None
            if result.workflow_provenance is None
            else thaw_json_mapping(
                result.workflow_provenance,
                field_name="differential_result.workflow_provenance",
            )
        ),
        "policy_provenance": (
            None
            if result.policy_provenance is None
            else _policy_provenance_payload(result.policy_provenance)
        ),
        "empirical_bayes": {
            "method": result.empirical_bayes_method,
            "robust": result.empirical_bayes_robust,
            "trend": result.empirical_bayes_trend,
            "prior_variance": _json_scalar(result.prior_variance),
            "prior_degrees_of_freedom": _json_scalar(result.prior_degrees_of_freedom),
            "residual_degrees_of_freedom": _json_scalar(
                result.residual_degrees_of_freedom
            ),
            "prior_residual_variance_by_feature": _series_records_payload(
                result.prior_residual_variance
            ),
            "prior_degrees_of_freedom_by_feature": _series_records_payload(
                result.prior_degrees_of_freedom_series_value
            ),
        },
        **(
            {}
            if result.mean_variance_trend_diagnostics is None
            else {
                "mean_variance_trend_diagnostics": _mean_variance_trend_payload(
                    result.mean_variance_trend_diagnostics
                )
            }
        ),
        **(
            {}
            if result.quantification_depth_trend_diagnostics is None
            else {
                "quantification_depth_trend_diagnostics": (
                    _quantification_depth_trend_payload(
                        result.quantification_depth_trend_diagnostics
                    )
                )
            }
        ),
        "contrast_tables": {
            contrast_name: _dataframe_records_payload(table)
            for contrast_name, table in contrast_tables.items()
        },
        **(
            {}
            if feature_eligibility is None
            else {
                "feature_eligibility": _dataframe_records_payload(feature_eligibility)
            }
        ),
    }


def _mean_variance_trend_payload(
    diagnostics: MeanVarianceTrendDiagnostics,
) -> dict[str, object]:
    return {
        "mean_intensity": _series_records_payload(diagnostics.mean_intensity),
        "trend_covariate": _series_records_payload(diagnostics.trend_covariate),
        "trend_covariate_name": diagnostics.trend_covariate_name,
        "trend_covariate_transformation": diagnostics.trend_covariate_transformation,
        "log_residual_variance": _series_records_payload(
            diagnostics.log_residual_variance
        ),
        "fitted_log_prior_variance": _series_records_payload(
            diagnostics.fitted_log_prior_variance
        ),
    }


def _quantification_depth_trend_payload(
    diagnostics: QuantificationDepthTrendDiagnostics,
) -> dict[str, object]:
    return {
        "quantification_depth": _series_records_payload(
            diagnostics.quantification_depth
        ),
        "trend_covariate": _series_records_payload(diagnostics.trend_covariate),
        "trend_covariate_name": diagnostics.trend_covariate_name,
        "trend_covariate_transformation": diagnostics.trend_covariate_transformation,
        "quantification_depth_kind": diagnostics.quantification_depth_kind,
        "log_residual_variance": _series_records_payload(
            diagnostics.log_residual_variance
        ),
        "fitted_log_prior_variance": _series_records_payload(
            diagnostics.fitted_log_prior_variance
        ),
    }


def _series_records_payload(series: pd.Series) -> list[dict[str, object]]:
    return [
        {
            "feature_id": _json_scalar(index_value),
            "value": _json_scalar(value),
        }
        for index_value, value in series.items()
    ]


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


def _policy_provenance_payload(
    policy_provenance: DifferentialPolicyProvenance,
) -> object:
    payload = asdict(policy_provenance)
    if policy_provenance.protein_aware is None:
        payload.pop("protein_aware", None)
    return _json_payload(payload)
