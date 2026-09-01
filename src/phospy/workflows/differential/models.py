"""Differential workflow stage-boundary models and contracts."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, NoReturn, Protocol, cast

import pandas as pd

from phospy.contracts.configs import DifferentialAnalysisConfig
from phospy.contracts.result_caveats import ResultCaveat
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.configs.differential import (
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS,
    DifferentialImputedValuePolicy,
    DifferentialProteinAwareModelMethod,
    DifferentialReliabilityProfile,
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
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationRequest,
)
from phospy.science.differential.policy_models import TechnicalReplicatePolicy
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregationPlan,
)

if TYPE_CHECKING:
    from phospy.science.datasets.internal_view import DatasetInternalView
    from phospy.validation.workflows.differential import (
        ValidatedExperimentalDesignContract,
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
    reliability_profile: DifferentialReliabilityProfile
    minimum_condition_replicates: int
    empirical_bayes: EmpiricalBayesConfig
    multiple_testing_method: MultipleTestingMethod
    protein_aware_method: DifferentialProteinAwareModelMethod | None = None

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
            "reliability_profile",
            str(self.reliability_profile),
        )
        object.__setattr__(
            self,
            "minimum_condition_replicates",
            int(self.minimum_condition_replicates),
        )
        if self.protein_aware_method is not None:
            if (
                self.protein_aware_method
                not in SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS
            ):
                supported = ", ".join(
                    repr(value)
                    for value in SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS
                )
                raise WorkflowBoundaryError(
                    seam="differential.execution_config.protein_aware_method",
                    next_action=(
                        "resolve protein-aware differential execution from a "
                        "supported DifferentialProteinAwareModelConfig"
                    ),
                    details={"supported_methods": supported},
                    message_prefix=("differential workflow boundary validation failed"),
                )
            object.__setattr__(
                self,
                "protein_aware_method",
                str(self.protein_aware_method),
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
    block_ids: tuple[str, ...] | None
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
        block_ids = (
            None
            if self.block_ids is None
            else tuple(str(value) for value in self.block_ids)
        )
        if block_ids is not None and len(block_ids) != len(self.sample_order):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.block_id_alignment",
                next_action=(
                    "carry duplicate-correlation block IDs in the exact execution "
                    "sample order"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if (
            self.paired_design_policy == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
            and block_ids is None
        ):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.duplicate_correlation_blocks",
                next_action=(
                    "retain the validated block_id vector separately from the "
                    "fixed-effects design"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if (
            self.paired_design_policy == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
            and self.block_column_metadata is not None
        ):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.duplicate_correlation_design",
                next_action=(
                    "do not combine fixed block coefficients with a block "
                    "correlation structure"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if (
            self.paired_design_policy == PAIRED_DESIGN_POLICY_FIXED_BLOCK
            and block_ids is not None
        ):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.fixed_block_block_vector",
                next_action=(
                    "represent fixed_block solely through fixed-effect block columns"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        object.__setattr__(self, "block_ids", block_ids)
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
class ProteinAwareDifferentialResolvedInputs:
    """Execution-ready protein-aware inputs resolved before statistical fitting."""

    computation_request: ProteinAwareDifferentialComputationRequest
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs
    matched_pairs: pd.DataFrame
    candidate_matched_pairs: pd.DataFrame
    resolved_protein_covariates: pd.DataFrame
    site_eligibility_metadata: pd.DataFrame
    full_site_ids: tuple[str, ...]
    tested_site_ids: tuple[str, ...]
    sample_order: tuple[str, ...]
    base_design: DesignMatrix
    base_contrasts: ContrastMatrix
    method_id: DifferentialProteinAwareModelMethod
    preparation_policy: str
    protein_mapping_policy: str
    eligibility_counts: tuple[tuple[str, int], ...]
    status_counts: tuple[tuple[str, int], ...]
    reason_counts: tuple[tuple[str, int], ...]

    def __post_init__(self) -> None:
        if not isinstance(
            cast(object, self.computation_request),
            ProteinAwareDifferentialComputationRequest,
        ):
            raise WorkflowBoundaryError(
                seam="differential.protein_aware_inputs.computation_request",
                next_action=(
                    "resolve a ProteinAwareDifferentialComputationRequest before "
                    "protein-aware differential fitting"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if not isinstance(
            cast(object, self.feature_eligibility_inputs),
            DifferentialFeatureEligibilityInputs,
        ):
            raise WorkflowBoundaryError(
                seam="differential.protein_aware_inputs.feature_eligibility",
                next_action=(
                    "carry full-index protein-aware eligibility metadata with "
                    "resolved computation inputs"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if self.method_id not in SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS:
            raise WorkflowBoundaryError(
                seam="differential.protein_aware_inputs.method",
                next_action="select a supported protein-aware differential method",
                details={"method": str(self.method_id)},
                message_prefix="differential workflow boundary validation failed",
            )
        if not isinstance(cast(object, self.base_design), DesignMatrix):
            _raise_protein_aware_resolved_inputs_error(
                seam="base_design",
                next_action=(
                    "carry the validated ordinary design matrix with resolved "
                    "protein-aware inputs"
                ),
            )
        if not isinstance(cast(object, self.base_contrasts), ContrastMatrix):
            _raise_protein_aware_resolved_inputs_error(
                seam="base_contrasts",
                next_action=(
                    "carry the validated ordinary contrast matrix with resolved "
                    "protein-aware inputs"
                ),
            )
        full_site_ids = tuple(str(value) for value in self.full_site_ids)
        tested_site_ids = tuple(str(value) for value in self.tested_site_ids)
        sample_order = tuple(str(value) for value in self.sample_order)
        if not full_site_ids:
            _raise_protein_aware_resolved_inputs_error(
                seam="site_ids",
                next_action="carry the full phosphosite index with resolved inputs",
            )
        if len(set(full_site_ids)) != len(full_site_ids):
            _raise_protein_aware_resolved_inputs_error(
                seam="site_ids",
                next_action="carry unique full phosphosite identifiers",
                details={"duplicate_full_site_ids": _duplicates(full_site_ids)[:5]},
            )
        if len(set(tested_site_ids)) != len(tested_site_ids):
            _raise_protein_aware_resolved_inputs_error(
                seam="tested_site_ids",
                next_action="carry unique tested phosphosite identifiers",
                details={"duplicate_tested_site_ids": _duplicates(tested_site_ids)[:5]},
            )
        full_site_id_set = set(full_site_ids)
        unexpected_tested_ids = [
            site_id for site_id in tested_site_ids if site_id not in full_site_id_set
        ]
        if unexpected_tested_ids:
            _raise_protein_aware_resolved_inputs_error(
                seam="tested_site_ids",
                next_action="carry tested phosphosite identifiers from the full index",
                details={"unexpected_tested_site_ids": unexpected_tested_ids[:5]},
            )
        tested_site_id_set = set(tested_site_ids)
        expected_tested_order = tuple(
            site_id for site_id in full_site_ids if site_id in tested_site_id_set
        )
        if tested_site_ids != expected_tested_order:
            _raise_protein_aware_resolved_inputs_error(
                seam="tested_site_order",
                next_action=(
                    "preserve full-index phosphosite order when carrying tested "
                    "protein-aware identifiers"
                ),
                details={
                    "tested_site_ids": list(tested_site_ids),
                    "expected_tested_site_ids": list(expected_tested_order),
                },
            )
        if self.computation_request.sample_order != sample_order:
            _raise_protein_aware_resolved_inputs_error(
                seam="sample_order",
                next_action=(
                    "carry one validated analysis sample order through the resolved "
                    "artifact and computation request"
                ),
                details={
                    "sample_order": list(sample_order),
                    "computation_sample_order": list(
                        self.computation_request.sample_order
                    ),
                },
            )
        if self.method_id != self.computation_request.method_id:
            _raise_protein_aware_resolved_inputs_error(
                seam="method",
                next_action=(
                    "carry one supported protein-aware method through the resolved "
                    "artifact and computation request"
                ),
                details={
                    "method_id": str(self.method_id),
                    "computation_method_id": str(self.computation_request.method_id),
                },
            )
        preparation_policy = _require_non_empty_protein_aware_text(
            self.preparation_policy,
            seam="preparation_policy",
            label="protein-aware preparation policy",
        )
        protein_mapping_policy = _require_non_empty_protein_aware_text(
            self.protein_mapping_policy,
            seam="protein_mapping_policy",
            label="protein-aware mapping policy",
        )

        feature_metadata = pd.DataFrame(
            self.feature_eligibility_inputs.feature_metadata,
            copy=True,
        )
        result_status = pd.Series(
            self.feature_eligibility_inputs.result_status,
            copy=True,
        )
        feature_eligibility_inputs = DifferentialFeatureEligibilityInputs(
            feature_metadata=feature_metadata,
            result_status=result_status,
            testable_feature_ids=tested_site_ids,
            attach_to_result_tables=(
                self.feature_eligibility_inputs.attach_to_result_tables
            ),
        )
        matched_pairs = pd.DataFrame(self.matched_pairs, copy=True)
        candidate_matched_pairs = pd.DataFrame(
            self.candidate_matched_pairs,
            copy=True,
        )
        resolved_protein_covariates = pd.DataFrame(
            self.resolved_protein_covariates,
            copy=True,
        )
        site_eligibility_metadata = pd.DataFrame(
            self.site_eligibility_metadata,
            copy=True,
        )

        _require_index_labels(
            feature_eligibility_inputs.feature_metadata.index,
            expected=full_site_ids,
            seam="feature_eligibility_metadata",
            label="feature eligibility metadata",
        )
        _require_index_labels(
            feature_eligibility_inputs.result_status.index,
            expected=full_site_ids,
            seam="feature_eligibility_status",
            label="feature eligibility result status",
        )
        if self.feature_eligibility_inputs.testable_feature_ids != tested_site_ids:
            _raise_protein_aware_resolved_inputs_error(
                seam="feature_eligibility_tested_ids",
                next_action=(
                    "carry the same tested phosphosite identifiers in eligibility "
                    "and computation inputs"
                ),
                details={
                    "feature_testable_site_ids": list(
                        self.feature_eligibility_inputs.testable_feature_ids
                    ),
                    "tested_site_ids": list(tested_site_ids),
                },
            )
        _require_index_labels(
            site_eligibility_metadata.index,
            expected=full_site_ids,
            seam="site_eligibility_metadata",
            label="site eligibility metadata",
        )
        if not site_eligibility_metadata.equals(feature_metadata):
            _raise_protein_aware_resolved_inputs_error(
                seam="site_eligibility_metadata",
                next_action=(
                    "carry the full-index protein-aware eligibility metadata used "
                    "by feature eligibility"
                ),
            )
        _require_index_labels(
            self.computation_request.phosphosite_matrix.index,
            expected=tested_site_ids,
            seam="computation_phosphosite_matrix",
            label="computation phosphosite matrix",
        )
        _require_index_labels(
            self.computation_request.phosphosite_matrix.columns,
            expected=sample_order,
            seam="sample_order",
            label="computation phosphosite matrix columns",
        )
        _require_index_labels(
            self.computation_request.base_design.frame.index,
            expected=sample_order,
            seam="base_design",
            label="base design rows",
        )
        _require_index_labels(
            self.computation_request.resolved_protein_covariates.columns,
            expected=sample_order,
            seam="resolved_protein_covariates",
            label="resolved protein covariate columns",
        )
        computation_design_frame = pd.DataFrame(
            self.computation_request.base_design.frame,
            copy=False,
        )
        base_design_frame = pd.DataFrame(self.base_design.frame, copy=False)
        computation_contrast_frame = pd.DataFrame(
            self.computation_request.base_contrasts.frame,
            copy=False,
        )
        base_contrast_frame = pd.DataFrame(self.base_contrasts.frame, copy=False)
        if not computation_design_frame.equals(base_design_frame):
            _raise_protein_aware_resolved_inputs_error(
                seam="base_design",
                next_action=(
                    "carry the same ordinary design matrix into the resolved artifact "
                    "and computation request"
                ),
            )
        if not computation_contrast_frame.equals(base_contrast_frame):
            _raise_protein_aware_resolved_inputs_error(
                seam="base_contrasts",
                next_action=(
                    "carry the same ordinary contrast matrix into the resolved "
                    "artifact and computation request"
                ),
            )
        if not self.computation_request.matched_pairs.equals(matched_pairs):
            _raise_protein_aware_resolved_inputs_error(
                seam="matched_pairs",
                next_action=(
                    "carry the same tested site-to-protein matches into the resolved "
                    "artifact and computation request"
                ),
            )
        if not self.computation_request.resolved_protein_covariates.equals(
            resolved_protein_covariates
        ):
            _raise_protein_aware_resolved_inputs_error(
                seam="resolved_protein_covariates",
                next_action=(
                    "carry the same resolved protein covariates into the resolved "
                    "artifact and computation request"
                ),
            )
        matched_site_ids = _frame_column_strings(
            matched_pairs,
            "site_key",
            seam="matched_pairs",
        )
        if matched_site_ids != tested_site_ids:
            _raise_protein_aware_resolved_inputs_error(
                seam="matched_pairs",
                next_action=(
                    "align tested matched pairs one-to-one with tested phosphosite "
                    "identifiers"
                ),
                details={
                    "matched_pair_site_ids": list(matched_site_ids),
                    "tested_site_ids": list(tested_site_ids),
                },
            )
        candidate_site_ids = _frame_column_strings(
            candidate_matched_pairs,
            "site_key",
            seam="candidate_matched_pairs",
        )
        if len(set(candidate_site_ids)) != len(candidate_site_ids):
            _raise_protein_aware_resolved_inputs_error(
                seam="candidate_matched_pairs",
                next_action="carry candidate matched pairs with unique site_key values",
                details={
                    "duplicate_candidate_site_ids": _duplicates(candidate_site_ids)[:5]
                },
            )
        unexpected_candidate_ids = [
            site_id for site_id in candidate_site_ids if site_id not in full_site_id_set
        ]
        if unexpected_candidate_ids:
            _raise_protein_aware_resolved_inputs_error(
                seam="candidate_matched_pairs",
                next_action="carry candidate matched pairs from the full site index",
                details={"unexpected_candidate_site_ids": unexpected_candidate_ids[:5]},
            )
        candidate_site_id_set = set(candidate_site_ids)
        expected_candidate_order = tuple(
            site_id for site_id in full_site_ids if site_id in candidate_site_id_set
        )
        if candidate_site_ids != expected_candidate_order:
            _raise_protein_aware_resolved_inputs_error(
                seam="candidate_matched_pairs",
                next_action=(
                    "preserve full-index phosphosite order when carrying candidate "
                    "matched pairs"
                ),
                details={
                    "candidate_site_ids": list(candidate_site_ids),
                    "expected_candidate_site_ids": list(expected_candidate_order),
                },
            )
        tested_total_row_keys = _frame_column_strings(
            matched_pairs,
            "total_protein_row_key",
            seam="matched_pairs",
        )
        expected_covariate_index = tuple(dict.fromkeys(tested_total_row_keys))
        _require_index_labels(
            resolved_protein_covariates.index,
            expected=expected_covariate_index,
            seam="resolved_protein_covariates",
            label="resolved protein covariate rows",
        )
        object.__setattr__(self, "full_site_ids", full_site_ids)
        object.__setattr__(self, "tested_site_ids", tested_site_ids)
        object.__setattr__(self, "sample_order", sample_order)
        object.__setattr__(
            self,
            "feature_eligibility_inputs",
            feature_eligibility_inputs,
        )
        object.__setattr__(
            self,
            "matched_pairs",
            matched_pairs,
        )
        object.__setattr__(
            self,
            "candidate_matched_pairs",
            candidate_matched_pairs,
        )
        object.__setattr__(
            self,
            "resolved_protein_covariates",
            resolved_protein_covariates,
        )
        object.__setattr__(
            self,
            "site_eligibility_metadata",
            site_eligibility_metadata,
        )
        object.__setattr__(
            self,
            "eligibility_counts",
            tuple((str(name), int(count)) for name, count in self.eligibility_counts),
        )
        object.__setattr__(
            self,
            "status_counts",
            tuple((str(status), int(count)) for status, count in self.status_counts),
        )
        object.__setattr__(
            self,
            "reason_counts",
            tuple((str(reason), int(count)) for reason, count in self.reason_counts),
        )
        object.__setattr__(self, "preparation_policy", preparation_policy)
        object.__setattr__(self, "protein_mapping_policy", protein_mapping_policy)


def _require_index_labels(
    index: pd.Index,
    *,
    expected: tuple[str, ...],
    seam: str,
    label: str,
) -> None:
    actual = tuple(str(value) for value in index.tolist())
    if actual == expected:
        return
    _raise_protein_aware_resolved_inputs_error(
        seam=seam,
        next_action=f"align {label} with resolved protein-aware inputs",
        details={"actual": list(actual), "expected": list(expected)},
    )


def _frame_column_strings(
    frame: pd.DataFrame,
    column: str,
    *,
    seam: str,
) -> tuple[str, ...]:
    if column not in frame.columns:
        _raise_protein_aware_resolved_inputs_error(
            seam=seam,
            next_action=f"carry a {column!r} column in resolved protein-aware inputs",
        )
    return tuple(str(value) for value in frame.loc[:, column].tolist())


def _duplicates(values: tuple[str, ...]) -> list[str]:
    return [value for value in dict.fromkeys(values) if values.count(value) > 1]


def _require_non_empty_protein_aware_text(
    value: object,
    *,
    seam: str,
    label: str,
) -> str:
    if value is None:
        _raise_protein_aware_resolved_inputs_error(
            seam=seam,
            next_action=f"carry a non-empty {label} with resolved inputs",
        )
    text = str(value).strip()
    if not text:
        _raise_protein_aware_resolved_inputs_error(
            seam=seam,
            next_action=f"carry a non-empty {label} with resolved inputs",
        )
    return text


def _raise_protein_aware_resolved_inputs_error(
    *,
    seam: str,
    next_action: str,
    details: dict[str, object] | None = None,
) -> NoReturn:
    raise WorkflowBoundaryError(
        seam=f"differential.protein_aware_inputs.{seam}",
        next_action=next_action,
        details=details,
        message_prefix="differential workflow boundary validation failed",
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


class DifferentialDatasetEligibilityValidatorContract(Protocol):
    """Internal collaborator contract for differential dataset eligibility."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        imputed_value_policy: DifferentialImputedValuePolicy,
        allow_suspicious_declared_input_scale: bool,
        dataset_view: DatasetInternalView,
    ) -> None: ...


class DifferentialTechnicalReplicatePlannerContract(Protocol):
    """Internal collaborator contract for differential replicate planning."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        technical_replicate_policy: TechnicalReplicatePolicy,
        dataset_view: DatasetInternalView,
    ) -> TechnicalReplicateAggregationPlan: ...


class DifferentialDesignValidatorContract(Protocol):
    """Internal collaborator contract for differential design validation."""

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        contrasts: tuple[Contrast, ...],
        allow_design_subset: bool,
        minimum_condition_replicates: int,
        paired_design_policy: PairedDesignPolicy,
        dataset_view: DatasetInternalView,
    ) -> ValidatedExperimentalDesignContract: ...


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
    "DifferentialDatasetEligibilityValidatorContract",
    "DifferentialDesignValidatorContract",
    "DifferentialAnalysisInterpreterContract",
    "DifferentialTechnicalReplicatePlannerContract",
    "DifferentialAnalysisValidatorContract",
    "InterpretedDifferentialAnalysisRequest",
    "ProteinAwareDifferentialResolvedInputs",
    "ValidatedDifferentialAnalysisRequest",
]
