"""Dataset builder contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.science.evidence.dataset_resolution import (
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
)
from phospy.science.references.models import Organism
from phospy.science.sites.identifiers import SiteIdentifierNormalisationReport
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)

if TYPE_CHECKING:
    from phospy.api.requests import DatasetBuildRequest
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
    site_sequence_derivation: dict[str, object] | None = None
    quantitative_meaning: QuantitativeMeaning | None = None
    site_resolution_mode: str = DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED
    multi_site_policy: str | None = None
    allow_opaque_site_values: bool = False
    peptide_evidence_resolution: dict[str, object] | None = None


@dataclass(frozen=True, slots=True)
class PreprocessedDatasetBuildTables:
    """Tables after internal preprocessing and before scale-state establishment."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    comparisons: pd.DataFrame | None = None
    comparison_group_stats: pd.DataFrame | None = None
    comparison_pair_stats: pd.DataFrame | None = None
    preprocessing_row_counts: pd.DataFrame | None = None
    preprocessing_operations: pd.DataFrame | None = None
    row_audit: pd.DataFrame | None = None
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None = None
    duplicate_site_resolution: pd.DataFrame | None = None
    metadata_conflicts: pd.DataFrame | None = None


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
    ) -> PreprocessedDatasetBuildTables: ...


__all__ = [
    "DatasetBuildExecutorContract",
    "DatasetBuildInterpreterContract",
    "DatasetPreprocessorContract",
    "DatasetBuildValidatorContract",
    "DatasetInput",
    "InterpretedDatasetBuildRequest",
    "PreprocessedDatasetBuildTables",
]
