"""Public request models.

Request dataclasses are lightweight command payloads. They intentionally do not
perform scientific validation during construction; the relevant dataset builder
or workflow validator enforces request compatibility before interpretation and
execution. Config dataclasses may still self-validate local configuration
invariants.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from phospy.contracts.configs import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    DatasetPreprocessingConfig,
    DifferentialAnalysisConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseReferenceDisplayAmbiguityPolicy,
    KinaseScoringConfig,
    KinaseSiteSequenceConflictPolicy,
    SignalomeConfig,
)
from phospy.science.datasets.builders.contracts import DatasetInput
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.models import (
    Contrast,
    ExperimentalDesign,
    SampleDesignRecord,
)
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    EmpiricalBayesConfig,
)
from phospy.science.evidence import dataset_resolution as _dataset_resolution
from phospy.science.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)

if TYPE_CHECKING:
    from phospy.contracts.results import KinaseWorkflowResult

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


@dataclass(frozen=True, slots=True)
class DatasetBuildRequest:
    """Request for building an ``AnalysisReadyPhosphoDataset``.

    Construction stores the payload only. Dataset input shape, source types,
    preprocessing compatibility, site-resolution mode, and scientific
    readiness are validated by ``AnalysisReadyDatasetBuilder.run(...)`` before
    interpretation or execution.

    Supported public inputs are pandas ``DataFrame`` values or file paths.
    ``site_resolution_mode`` selects one of two explicit lanes:

    - ``site_level_resolved``: provide ``phospho`` + ``site_metadata``.
    - ``peptide_evidence``: provide ``peptide_evidence`` plus
      ``peptide_evidence_sample_intensity_columns`` and ``multi_site_policy``.

    Preprocessing policy remains builder owned via ``preprocessing_config`` and
    must still converge on a strict, missing-value-free
    ``AnalysisReadyPhosphoDataset`` boundary. Opaque non-STY site tokens remain
    disallowed by default and require explicit ``allow_opaque_site_values=True``.
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
    allow_opaque_site_values: bool = False
    organism: Organism | None = None
    preprocessing_config: DatasetPreprocessingConfig = field(
        default_factory=DatasetPreprocessingConfig
    )
    input_intensity_scale: IntensityScaleKind | str | None = None
    quantitative_meaning: QuantitativeMeaning | str | None = None


@dataclass(frozen=True, slots=True)
class KinaseWorkflowRequest:
    """Request for the public kinase workflow.

    Construction stores the payload only. Dataset compatibility, reference
    compatibility, scoring/prediction/activity config coherence, localisation
    requirements, and reference-projection policies are validated by
    ``KinaseWorkflow.run(...)`` before interpretation or execution.

    `site_sequence_conflict_policy` controls how dataset and reference
    site-sequence disagreements are handled during interpretation.

    `reference_display_ambiguity_policy` controls whether display-level
    kinase-substrate reference rows may project to more than one dataset
    `site_key`. The default is conservative and rejects ambiguous projection.
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
    reference_display_ambiguity_policy: KinaseReferenceDisplayAmbiguityPolicy = (
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR
    )


@dataclass(frozen=True, slots=True)
class SignalomeWorkflowRequest:
    """Request for the public signalome workflow.

    Construction stores the payload only. The upstream kinase result, signalome
    config, score/prediction matrix alignment, site identity, and protein
    grouping metadata are validated by ``SignalomeWorkflow.run(...)`` before
    interpretation or execution.

    Signalome execution requires explicit protein grouping metadata per
    interpreted site via ``dataset.site_metadata.protein_id``.
    """

    kinase_result: KinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)


@dataclass(frozen=True, slots=True)
class DifferentialAnalysisRequest:
    """Request for the public differential workflow.

    Construction stores the payload only. Dataset eligibility, design/sample
    alignment, contrast validity, replicate requirements, and config coherence
    are validated by ``DifferentialAnalysisWorkflow.run(...)`` before
    interpretation or statistical execution.
    """

    dataset: AnalysisReadyPhosphoDataset
    design: ExperimentalDesign
    contrasts: tuple[Contrast, ...]
    config: DifferentialAnalysisConfig = field(
        default_factory=DifferentialAnalysisConfig
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
    "SignalomeWorkflowRequest",
]
