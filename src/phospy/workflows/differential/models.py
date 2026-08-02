"""Differential workflow stage-boundary models and contracts."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.contracts.configs import DifferentialAnalysisConfig
from phospy.contracts.result_caveats import ResultCaveat
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.configs.differential import (
    DifferentialImputedValuePolicy,
    MultipleTestingMethod,
    PairedDesignPolicy,
)
from phospy.science.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
from phospy.science.design.matrix_builder import DesignMatrixBuildResult
from phospy.science.design.models import Contrast, ExperimentalDesign
from phospy.science.differential.linear_model import (
    DifferentialDesignDecomposition,
    DifferentialDesignDecompositionError,
)
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisResult,
    DifferentialPolicyProvenance,
    EmpiricalBayesConfig,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.science.differential.policy_models import TechnicalReplicatePolicy
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregationPlan,
)

if TYPE_CHECKING:
    from phospy.science.datasets.internal_view import DatasetInternalView


@dataclass(frozen=True, slots=True)
class ValidatedDifferentialAnalysisRequest:
    """Validated differential request passed to interpretation."""

    dataset: AnalysisReadyPhosphoDataset
    design: ExperimentalDesign
    contrasts: tuple[Contrast, ...]
    analysis_sample_ids: tuple[str, ...]
    design_matrix: DesignMatrix
    contrast_matrix: ContrastMatrix
    design_decomposition: DifferentialDesignDecomposition
    config: DifferentialAnalysisConfig
    policy_provenance: DifferentialPolicyProvenance | None = None
    technical_replicate_aggregation_plan: TechnicalReplicateAggregationPlan | None = (
        None
    )
    workflow_provenance: Mapping[str, object] | None = None
    dataset_preprocessing_report: DatasetPreprocessingReport | None = None
    design_build_result: DesignMatrixBuildResult | None = None
    dataset_view: DatasetInternalView | None = field(
        default=None,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        _require_decomposition_matches_design_matrix(
            design_decomposition=self.design_decomposition,
            design_matrix=self.design_matrix,
            seam="differential.validator.design_decomposition_identity",
        )


@dataclass(frozen=True, slots=True)
class ResolvedDifferentialExecutionConfig:
    """Execution-ready differential policy resolved from public configuration."""

    technical_replicate_policy: TechnicalReplicatePolicy
    paired_design_policy: PairedDesignPolicy
    imputed_value_policy: DifferentialImputedValuePolicy
    imputed_value_max_fraction: float
    allow_design_subset: bool
    allow_suspicious_declared_input_scale: bool
    minimum_condition_replicates: int
    empirical_bayes: EmpiricalBayesConfig
    multiple_testing_method: MultipleTestingMethod

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "paired_design_policy",
            str(self.paired_design_policy),
        )
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
            "allow_design_subset",
            bool(self.allow_design_subset),
        )
        object.__setattr__(
            self,
            "allow_suspicious_declared_input_scale",
            bool(self.allow_suspicious_declared_input_scale),
        )
        object.__setattr__(
            self,
            "minimum_condition_replicates",
            int(self.minimum_condition_replicates),
        )


@dataclass(frozen=True, slots=True)
class DifferentialImputationPolicyInputs:
    """Aligned imputation-policy inputs resolved before execution."""

    feature_metadata: pd.DataFrame
    result_status: pd.Series
    result_status_reason: pd.Series
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
class DifferentialFeatureEligibilityInputs:
    """Aligned feature-level eligibility inputs resolved before execution."""

    feature_metadata: pd.DataFrame
    result_status: pd.Series
    testable_feature_ids: tuple[str, ...]
    attach_to_result_tables: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "testable_feature_ids",
            tuple(str(value) for value in self.testable_feature_ids),
        )
        object.__setattr__(
            self,
            "attach_to_result_tables",
            bool(self.attach_to_result_tables),
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
    design_decomposition: DifferentialDesignDecomposition

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
        _require_decomposition_matches_design_matrix(
            design_decomposition=self.design_decomposition,
            design_matrix=self.design_matrix,
            seam="differential.interpreter.execution_design_decomposition",
        )


@dataclass(frozen=True, slots=True)
class InterpretedDifferentialAnalysisRequest:
    """Execution-ready differential request produced by interpretation."""

    computation_request: DifferentialComputationRequest
    result_identity_metadata: pd.DataFrame
    config: DifferentialAnalysisConfig
    execution_config: ResolvedDifferentialExecutionConfig
    design_rank: int
    residual_degrees_of_freedom: float
    design_decomposition: DifferentialDesignDecomposition
    policy_provenance: DifferentialPolicyProvenance | None = None
    workflow_provenance: Mapping[str, object] | None = None
    caveats: tuple[ResultCaveat, ...] = ()
    dataset_preprocessing_report: DatasetPreprocessingReport | None = None
    execution_design: DifferentialExecutionDesignInputs | None = None
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None = None
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None = None
    normalisation_state: str = "not_recorded"
    ruv_readiness_enabled: bool = False
    ruv_readiness_ready: bool = False

    def __post_init__(self) -> None:
        if (
            self.computation_request.design_decomposition
            is not self.design_decomposition
        ):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.computation_decomposition_identity",
                next_action=(
                    "pass the interpreted design decomposition into the "
                    "differential computation request without rebuilding it"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if (
            self.execution_design is not None
            and self.execution_design.design_decomposition
            is not self.design_decomposition
        ):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.execution_design_decomposition_identity",
                next_action=(
                    "assemble execution design metadata from the same interpreted "
                    "design decomposition object"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if int(self.design_rank) != int(self.design_decomposition.rank):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.design_rank_consistency",
                next_action=(
                    "derive interpreted design rank from the shared design "
                    "decomposition"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if not math.isclose(
            float(self.residual_degrees_of_freedom),
            float(self.design_decomposition.residual_degrees_of_freedom),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.residual_dof_consistency",
                next_action=(
                    "derive interpreted residual degrees of freedom from the "
                    "shared design decomposition"
                ),
                message_prefix="differential workflow boundary validation failed",
            )


def _require_decomposition_matches_design_matrix(
    *,
    design_decomposition: DifferentialDesignDecomposition,
    design_matrix: DesignMatrix,
    seam: str,
) -> None:
    try:
        design_decomposition.assert_matches_design(
            design_matrix.frame.to_numpy(dtype=float),
            field_name="differential.design",
        )
    except DifferentialDesignDecompositionError as error:
        raise WorkflowBoundaryError(
            seam=seam,
            next_action=(
                "reuse the design decomposition built for the validated "
                "differential design matrix"
            ),
            details={"error": str(error)},
            message_prefix="differential workflow boundary validation failed",
        ) from error


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
    "DifferentialFeatureEligibilityInputs",
    "DifferentialImputationPolicyInputs",
    "ResolvedDifferentialExecutionConfig",
    "DifferentialAnalysisExecutorContract",
    "DifferentialAnalysisInterpreterContract",
    "DifferentialAnalysisValidatorContract",
    "InterpretedDifferentialAnalysisRequest",
    "ValidatedDifferentialAnalysisRequest",
]
