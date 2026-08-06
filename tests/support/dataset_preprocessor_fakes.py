from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.science.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)
from phospy.science.transformations.quantitative_contracts import (
    QuantitativeContractState,
    initial_quantitative_contract_state,
)


@dataclass(frozen=True, slots=True)
class DatasetPreprocessorPreflightCall:
    plan: PreprocessingPlan
    initial_quantitative_scale_kind: IntensityScaleKind | None
    initial_quantitative_meaning: QuantitativeMeaning | None


@dataclass(frozen=True, slots=True)
class DatasetPreprocessorRunCall:
    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    plan: PreprocessingPlan
    corrected_preprocessing_output: CorrectedPreprocessingOutput | None
    initial_quantitative_scale_kind: IntensityScaleKind | None
    initial_quantitative_meaning: QuantitativeMeaning | None


class ConformingDatasetPreprocessorFake:
    """Shared test fake for the internal dataset preprocessor protocol."""

    def __init__(
        self,
        *,
        result: PreprocessedDatasetBuildTables | None = None,
        preflight_error: Exception | None = None,
    ) -> None:
        self.result = result
        self.preflight_error = preflight_error
        self.call_order: list[str] = []
        self.preflight_calls: list[DatasetPreprocessorPreflightCall] = []
        self.run_calls: list[DatasetPreprocessorRunCall] = []

    def validate_quantitative_contracts(
        self,
        *,
        plan: PreprocessingPlan,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> QuantitativeContractState:
        self.call_order.append("validate_quantitative_contracts")
        self.preflight_calls.append(
            DatasetPreprocessorPreflightCall(
                plan=plan,
                initial_quantitative_scale_kind=initial_quantitative_scale_kind,
                initial_quantitative_meaning=initial_quantitative_meaning,
            )
        )
        if self.preflight_error is not None:
            raise self.preflight_error
        return initial_quantitative_contract_state(
            declared_input_scale_kind=initial_quantitative_scale_kind,
            explicit_quantitative_meaning=initial_quantitative_meaning,
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        plan: PreprocessingPlan,
        corrected_preprocessing_output: CorrectedPreprocessingOutput | None = None,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> PreprocessedDatasetBuildTables:
        self.call_order.append("run")
        self.run_calls.append(
            DatasetPreprocessorRunCall(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                plan=plan,
                corrected_preprocessing_output=corrected_preprocessing_output,
                initial_quantitative_scale_kind=initial_quantitative_scale_kind,
                initial_quantitative_meaning=initial_quantitative_meaning,
            )
        )
        if self.result is not None:
            return self.result
        return PreprocessedDatasetBuildTables(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
        )
