"""Public request models."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING

from phospy.api.configs import (
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    DatasetPreprocessingConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseSiteSequenceConflictPolicy,
    SignalomeConfig,
)
from phospy.datasets.builders.contracts import DatasetInput
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.design.models import Contrast, ExperimentalDesign, SampleDesignRecord
from phospy.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    EmpiricalBayesConfig,
)
from phospy.errors.validation import WorkflowValidationError
from phospy.evidence import dataset_resolution as _dataset_resolution
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.transformations.models import IntensityScaleKind, QuantitativeMeaning

if TYPE_CHECKING:
    from phospy.api.results import KinaseWorkflowResult

DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING = (
    _dataset_resolution.DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
)
DATASET_MULTI_SITE_POLICY_KEEP_JOINT = (
    _dataset_resolution.DATASET_MULTI_SITE_POLICY_KEEP_JOINT
)
DATASET_MULTI_SITE_POLICY_REJECT = _dataset_resolution.DATASET_MULTI_SITE_POLICY_REJECT
DATASET_MULTI_SITE_POLICY_SPLIT = _dataset_resolution.DATASET_MULTI_SITE_POLICY_SPLIT
DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE = (
    _dataset_resolution.DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE
)
DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED = (
    _dataset_resolution.DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED
)

MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG = "benjamini_hochberg"
SUPPORTED_MULTIPLE_TESTING_METHODS: tuple[str, ...] = (
    MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG,
)


class TechnicalReplicatePolicy(str, Enum):
    """Policy for handling repeated biological replicate IDs."""

    REJECT = "reject"
    MEAN = "mean"
    MEDIAN = "median"


@dataclass(frozen=True, slots=True)
class MultipleTestingConfig:
    """Public multiple-testing policy for differential analysis."""

    method: str = MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG

    def __post_init__(self) -> None:
        if self.method not in SUPPORTED_MULTIPLE_TESTING_METHODS:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_MULTIPLE_TESTING_METHODS
            )
            raise WorkflowValidationError(
                f"differential.multiple_testing.method must be one of: {supported}"
            )


@dataclass(frozen=True, slots=True)
class DatasetBuildRequest:
    """Request for building an ``AnalysisReadyPhosphoDataset``.

    Supported public inputs are pandas ``DataFrame`` values or file paths.
    ``site_resolution_mode`` selects one of two explicit lanes:

    - ``site_level_resolved``: provide ``phospho`` + ``site_metadata``.
    - ``peptide_evidence``: provide ``peptide_evidence`` plus
      ``peptide_evidence_sample_intensity_columns`` and ``multi_site_policy``.

    Preprocessing policy remains builder owned via ``preprocessing_config`` and
    must still converge on a strict, missing-value-free
    ``AnalysisReadyPhosphoDataset`` boundary.
    """

    phospho: DatasetInput | None = None
    site_metadata: DatasetInput | None = None
    sample_metadata: DatasetInput | None = None
    total: DatasetInput | None = None
    site_resolution_mode: str = DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED
    peptide_evidence: DatasetInput | None = None
    peptide_evidence_sample_intensity_columns: tuple[str, ...] | None = None
    peptide_site_mapping: DatasetInput | None = None
    multi_site_policy: str | None = None
    organism: Organism | None = None
    preprocessing_config: DatasetPreprocessingConfig = field(
        default_factory=DatasetPreprocessingConfig
    )
    input_intensity_scale: IntensityScaleKind | str | None = None
    quantitative_meaning: QuantitativeMeaning | str | None = None


@dataclass(frozen=True, slots=True)
class KinaseWorkflowRequest:
    """Request for the public kinase workflow.

    `site_sequence_conflict_policy` controls how dataset and reference
    site-sequence disagreements are handled during interpretation.
    """

    dataset: AnalysisReadyPhosphoDataset
    references: ReferencePreset | ReferenceBundle = ReferencePreset.AUTO
    scoring_config: KinaseScoringConfig = field(default_factory=KinaseScoringConfig)
    prediction_config: KinasePredictionConfig = field(
        default_factory=KinasePredictionConfig
    )
    activity_config: KinaseActivityConfig | None = field(
        default_factory=KinaseActivityConfig
    )
    site_sequence_conflict_policy: KinaseSiteSequenceConflictPolicy = (
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE
    )


@dataclass(frozen=True, slots=True)
class SignalomeWorkflowRequest:
    """Request for the public signalome workflow.

    Signalome execution requires explicit protein identity per interpreted site via
    ``dataset.site_metadata.protein_id``.
    """

    kinase_result: KinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)


@dataclass(frozen=True, slots=True)
class DifferentialAnalysisRequest:
    """Request for the public differential workflow."""

    dataset: AnalysisReadyPhosphoDataset
    design: ExperimentalDesign
    contrasts: tuple[Contrast, ...]
    technical_replicate_policy: TechnicalReplicatePolicy | str = (
        TechnicalReplicatePolicy.REJECT
    )
    allow_design_subset: bool = False
    minimum_condition_replicates: int = 2
    empirical_bayes: EmpiricalBayesConfig = field(default_factory=EmpiricalBayesConfig)
    multiple_testing: MultipleTestingConfig = field(
        default_factory=MultipleTestingConfig
    )

    def __post_init__(self) -> None:
        design = self.design
        if not isinstance(design, ExperimentalDesign):
            raise WorkflowValidationError(
                "differential workflow request design must be ExperimentalDesign"
            )
        contrasts = self.contrasts
        if not isinstance(contrasts, tuple):
            contrasts = tuple(contrasts)  # type: ignore[arg-type]
        if not contrasts:
            raise WorkflowValidationError(
                "differential workflow request contrasts must include at least one Contrast"
            )
        for contrast in contrasts:
            if not isinstance(contrast, Contrast):
                raise WorkflowValidationError(
                    "differential workflow request contrasts must contain Contrast values"
                )
        technical_replicate_policy = self.technical_replicate_policy
        if isinstance(technical_replicate_policy, str):
            try:
                technical_replicate_policy = TechnicalReplicatePolicy(
                    technical_replicate_policy
                )
            except ValueError as error:
                supported = ", ".join(
                    sorted(policy.value for policy in TechnicalReplicatePolicy)
                )
                raise WorkflowValidationError(
                    "differential workflow request technical_replicate_policy must be "
                    f"one of: {supported}"
                ) from error
        if not isinstance(technical_replicate_policy, TechnicalReplicatePolicy):
            raise WorkflowValidationError(
                "differential workflow request technical_replicate_policy must be a "
                "TechnicalReplicatePolicy"
            )
        if not isinstance(self.allow_design_subset, bool):
            raise WorkflowValidationError(
                "differential workflow request allow_design_subset must be a bool"
            )
        if not isinstance(self.minimum_condition_replicates, int):
            raise WorkflowValidationError(
                "differential workflow request minimum_condition_replicates must be an int"
            )
        if self.minimum_condition_replicates < 1:
            raise WorkflowValidationError(
                "differential workflow request minimum_condition_replicates must be >= 1"
            )
        if not isinstance(self.empirical_bayes, EmpiricalBayesConfig):
            raise WorkflowValidationError(
                "differential workflow request empirical_bayes must be EmpiricalBayesConfig"
            )
        if not isinstance(self.multiple_testing, MultipleTestingConfig):
            raise WorkflowValidationError(
                "differential workflow request multiple_testing must be MultipleTestingConfig"
            )
        object.__setattr__(self, "design", design)
        object.__setattr__(self, "contrasts", contrasts)
        object.__setattr__(
            self,
            "technical_replicate_policy",
            technical_replicate_policy,
        )


__all__ = [
    "ContrastMatrix",
    "DesignMatrix",
    "DatasetBuildRequest",
    "SampleDesignRecord",
    "ExperimentalDesign",
    "Contrast",
    "DifferentialAnalysisRequest",
    "EmpiricalBayesConfig",
    "KinaseWorkflowRequest",
    "MultipleTestingConfig",
    "TechnicalReplicatePolicy",
    "SignalomeWorkflowRequest",
]
