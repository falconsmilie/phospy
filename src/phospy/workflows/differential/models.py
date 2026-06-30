"""Differential workflow stage-boundary models and contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.contracts.configs import DifferentialAnalysisConfig, PairedDesignPolicy
from phospy.science.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
from phospy.science.design.matrix_builder import DesignMatrixBuildResult
from phospy.science.design.models import Contrast, ExperimentalDesign
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisResult,
    DifferentialPolicyProvenance,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregationPlan,
)


@dataclass(frozen=True, slots=True)
class ValidatedDifferentialAnalysisRequest:
    """Validated differential request passed to interpretation."""

    dataset: AnalysisReadyPhosphoDataset
    design: ExperimentalDesign
    contrasts: tuple[Contrast, ...]
    analysis_sample_ids: tuple[str, ...]
    design_matrix: DesignMatrix
    contrast_matrix: ContrastMatrix
    config: DifferentialAnalysisConfig
    policy_provenance: DifferentialPolicyProvenance | None = None
    technical_replicate_aggregation_plan: TechnicalReplicateAggregationPlan | None = (
        None
    )
    workflow_provenance: Mapping[str, object] | None = None
    dataset_preprocessing_report: DatasetPreprocessingReport | None = None
    design_build_result: DesignMatrixBuildResult | None = None


@dataclass(frozen=True, slots=True)
class DifferentialImputationPolicyInputs:
    """Aligned imputation-policy inputs resolved before execution."""

    feature_metadata: pd.DataFrame
    result_status: pd.Series
    testable_feature_ids: tuple[str, ...]
    policy: str
    max_fraction: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "policy", str(self.policy))
        object.__setattr__(self, "max_fraction", float(self.max_fraction))
        object.__setattr__(
            self,
            "testable_feature_ids",
            tuple(str(value) for value in self.testable_feature_ids),
        )


@dataclass(frozen=True, slots=True)
class DifferentialCovariateColumnMetadata:
    """Resolved fixed-effect covariate encoding metadata for execution."""

    name: str
    kind: str
    columns: tuple[str, ...]
    levels: tuple[str, ...] = ()
    reference_level: str | None = None
    unused_levels: tuple[str, ...] = ()

    def __post_init__(self) -> None:
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
class DifferentialBlockColumnMetadata:
    """Resolved fixed-block encoding metadata for execution."""

    levels: tuple[str, ...]
    reference_level: str | None
    columns: tuple[tuple[str, str], ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "levels",
            tuple(str(value) for value in self.levels),
        )
        object.__setattr__(
            self,
            "reference_level",
            None if self.reference_level is None else str(self.reference_level),
        )
        object.__setattr__(
            self,
            "columns",
            tuple((str(level), str(column)) for level, column in self.columns),
        )


@dataclass(frozen=True, slots=True)
class DifferentialConditionContrastVector:
    """Resolved condition contrast vector aligned to design coefficients."""

    name: str
    numerator_condition: str
    denominator_condition: str
    coefficients: tuple[tuple[str, float], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(self, "numerator_condition", str(self.numerator_condition))
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


@dataclass(frozen=True, slots=True)
class DifferentialExecutionDesignInputs:
    """Execution-ready fixed-effect design inputs resolved by interpretation."""

    design_matrix: DesignMatrix
    contrast_matrix: ContrastMatrix
    condition_contrast_vectors: tuple[DifferentialConditionContrastVector, ...]
    covariate_columns: tuple[DifferentialCovariateColumnMetadata, ...]
    formula: str
    description: str
    sample_order: tuple[str, ...]
    paired_design_policy: PairedDesignPolicy
    block_column_metadata: DifferentialBlockColumnMetadata | None
    condition_labels: tuple[str, ...]
    coefficient_labels: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "condition_contrast_vectors",
            tuple(self.condition_contrast_vectors),
        )
        object.__setattr__(
            self,
            "covariate_columns",
            tuple(self.covariate_columns),
        )
        object.__setattr__(self, "formula", str(self.formula))
        object.__setattr__(self, "description", str(self.description))
        object.__setattr__(
            self,
            "paired_design_policy",
            str(self.paired_design_policy),
        )
        object.__setattr__(
            self,
            "sample_order",
            tuple(str(value) for value in self.sample_order),
        )
        object.__setattr__(
            self,
            "condition_labels",
            tuple(str(value) for value in self.condition_labels),
        )
        object.__setattr__(
            self,
            "coefficient_labels",
            tuple(str(value) for value in self.coefficient_labels),
        )


@dataclass(frozen=True, slots=True)
class InterpretedDifferentialAnalysisRequest:
    """Execution-ready differential request produced by interpretation."""

    computation_request: DifferentialComputationRequest
    result_identity_metadata: pd.DataFrame
    config: DifferentialAnalysisConfig
    design_rank: int
    residual_degrees_of_freedom: float
    policy_provenance: DifferentialPolicyProvenance | None = None
    workflow_provenance: Mapping[str, object] | None = None
    dataset_preprocessing_report: DatasetPreprocessingReport | None = None
    execution_design: DifferentialExecutionDesignInputs | None = None
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None = None
    normalisation_state: str = "not_recorded"
    ruv_readiness_enabled: bool = False
    ruv_readiness_ready: bool = False


class DifferentialAnalysisValidatorContract(Protocol):
    """Internal contract for differential workflow validation."""

    def run(self, request: object) -> ValidatedDifferentialAnalysisRequest: ...


class DifferentialAnalysisInterpreterContract(Protocol):
    """Internal contract for differential workflow interpretation."""

    def run(
        self, request: ValidatedDifferentialAnalysisRequest
    ) -> InterpretedDifferentialAnalysisRequest: ...


class DifferentialAnalysisExecutorContract(Protocol):
    """Internal contract for differential workflow execution."""

    def run(
        self, request: InterpretedDifferentialAnalysisRequest
    ) -> DifferentialAnalysisResult: ...


__all__ = [
    "DifferentialBlockColumnMetadata",
    "DifferentialConditionContrastVector",
    "DifferentialCovariateColumnMetadata",
    "DifferentialExecutionDesignInputs",
    "DifferentialImputationPolicyInputs",
    "DifferentialAnalysisExecutorContract",
    "DifferentialAnalysisInterpreterContract",
    "DifferentialAnalysisValidatorContract",
    "InterpretedDifferentialAnalysisRequest",
    "ValidatedDifferentialAnalysisRequest",
]
