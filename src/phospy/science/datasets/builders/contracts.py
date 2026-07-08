"""Dataset builder contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.science.datasets.builders.sequence_derivation import (
    SiteSequenceDerivationReport,
)
from phospy.science.datasets.builders.transformation_resolver import (
    ResolvedIntensityScale,
)
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.batch_correction_metadata import (
    ResolvedBatchCorrectionMetadata,
)
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.science.evidence.dataset_resolution import (
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    PeptideEvidenceResolutionSummary,
)
from phospy.science.references.models import Organism
from phospy.science.sites.identifiers import SiteIdentifierNormalisationReport
from phospy.science.transformations.models import (
    DeclaredIntensityScaleDiagnosticPolicy,
    IntensityScaleEstablishmentMode,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    IntensityTransformationEvent,
    QuantitativeMeaning,
)

if TYPE_CHECKING:
    from phospy.contracts.requests import DatasetBuildRequest
    from phospy.science.datasets.models import AnalysisReadyPhosphoDataset

DatasetInput = pd.DataFrame | str | Path | PathLike[str]


@dataclass(frozen=True, slots=True)
class InterpretedDatasetBuildRequest:
    """Resolved builder input after request interpretation."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    organism: Organism | None
    preprocessing_plan: PreprocessingPlan = field(
        default_factory=PreprocessingPlan.default
    )
    declared_input_intensity_scale_kind: IntensityScaleKind | None = None
    declared_input_intensity_scale_source: str | None = None
    site_identifier_normalisation: SiteIdentifierNormalisationReport | None = None
    site_sequence_derivation: SiteSequenceDerivationReport | None = None
    quantitative_meaning: QuantitativeMeaning | None = None
    site_resolution_mode: str = DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED
    multi_site_policy: str | None = None
    allow_opaque_site_values: bool = False
    allow_suspicious_declared_input_intensity_scale: bool = False
    peptide_evidence_resolution: PeptideEvidenceResolutionSummary | None = None
    corrected_preprocessing_output: CorrectedPreprocessingOutput | None = None


@dataclass(frozen=True, slots=True)
class PreprocessedDatasetBuildTables:
    """Tables after internal preprocessing and before scale-state establishment."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    comparisons: pd.DataFrame | None = None
    imputation_observation_mask: pd.DataFrame | None = None
    comparison_group_stats: pd.DataFrame | None = None
    comparison_pair_stats: pd.DataFrame | None = None
    preprocessing_row_counts: pd.DataFrame | None = None
    preprocessing_operations: pd.DataFrame | None = None
    row_audit: pd.DataFrame | None = None
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None = None
    duplicate_site_resolution: pd.DataFrame | None = None
    metadata_conflicts: pd.DataFrame | None = None
    batch_correction_metadata: ResolvedBatchCorrectionMetadata | None = None
    batch_correction_report: BatchCorrectionReport | None = None
    intensity_transformation_event: IntensityTransformationEvent | None = None


class DatasetBuildValidatorContract(Protocol):
    """Internal contract for dataset build request validation."""

    def run(self, request: DatasetBuildRequest) -> DatasetBuildRequest: ...


class DatasetBuildInterpreterContract(Protocol):
    """Internal contract for request interpretation into executable inputs."""

    def run(self, request: DatasetBuildRequest) -> InterpretedDatasetBuildRequest: ...


class DatasetBuildExecutorContract(Protocol):
    """Internal contract for constructing the final dataset."""

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset: ...


class DatasetPreprocessorContract(Protocol):
    """Internal contract for dataset preprocessing before scale-state setup."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        plan: PreprocessingPlan,
        corrected_preprocessing_output: CorrectedPreprocessingOutput | None = None,
    ) -> PreprocessedDatasetBuildTables: ...


class DatasetIntensityScaleResolverContract(Protocol):
    """Internal contract for dataset intensity-scale establishment."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None,
        expected_scale_kind: IntensityScaleKind | None = None,
        declared_input_scale_state: IntensityScaleState | None = None,
        declared_input_establishment_mode: IntensityScaleEstablishmentMode | None = (
            None
        ),
        input_declaration_source: str | None = None,
        scale_establishment_parameters: Mapping[str, object] | None = None,
        scale_establishment_evidence_level: IntensityScaleEvidenceLevel | None = None,
        establishment_transformer_name: str | None = None,
        establishment_trace_id: str | None = None,
        declared_scale_diagnostic_policy: DeclaredIntensityScaleDiagnosticPolicy
        | str = DeclaredIntensityScaleDiagnosticPolicy.WARN,
    ) -> ResolvedIntensityScale: ...


__all__ = [
    "DatasetBuildExecutorContract",
    "DatasetBuildInterpreterContract",
    "DatasetIntensityScaleResolverContract",
    "DatasetPreprocessorContract",
    "DatasetBuildValidatorContract",
    "DatasetInput",
    "InterpretedDatasetBuildRequest",
    "PreprocessedDatasetBuildTables",
]
