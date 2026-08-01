"""Preprocessing result and stage execution DTOs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.provenance.models import (
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    BatchCorrectionProvenance,
    DeterminismKind,
    ReproducibilityCaveat,
    TableFingerprint,
)
from phospy.science.datasets.preprocessing.report_schema import (
    ComparisonGroupStatsRow,
    ComparisonPairStatsRow,
    DuplicateSiteResolutionRow,
    MetadataConflictRow,
    PreprocessingRowAuditRow,
)
from phospy.science.transformations.models import IntensityTransformationEvent
from phospy.science.transformations.quantitative_contracts import (
    QuantitativeOperationContract,
    QuantitativeTransitionEvidence,
)

if TYPE_CHECKING:
    from phospy.science.datasets.preprocessing.trace import PreprocessingState

StageOwnedPreprocessingReportValue = (
    PreprocessingRowAuditRow
    | DuplicateSiteResolutionRow
    | MetadataConflictRow
    | ComparisonGroupStatsRow
    | ComparisonPairStatsRow
)


@dataclass(frozen=True, slots=True)
class ComparisonBuildResult:
    """Structured comparison-building output with provenance sidecars."""

    comparisons: pd.DataFrame
    comparison_group_stats: pd.DataFrame
    comparison_pair_stats: pd.DataFrame


@dataclass(frozen=True, slots=True)
class DuplicateSiteResolutionResult:
    """Duplicate-site policy output with structured provenance tables."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    dropped_row_count: int
    duplicate_site_resolution: pd.DataFrame
    metadata_conflicts: pd.DataFrame
    duplicate_aggregation_diagnostics: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class PreprocessingReportRow:
    """Structured stage-owned contribution to preprocessing report assembly."""

    table: str
    values: StageOwnedPreprocessingReportValue


@dataclass(frozen=True, slots=True)
class PreprocessingStageResult:
    """Structured output for a single preprocessing stage execution."""

    state: PreprocessingState
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    report_rows: Sequence[PreprocessingReportRow] = ()
    batch_correction_provenance: BatchCorrectionProvenance | None = None
    intensity_transformation_event: IntensityTransformationEvent | None = None
    quantitative_transition_evidence: QuantitativeTransitionEvidence | None = None


@dataclass(frozen=True, slots=True)
class PreprocessingStageExecution:
    """Executed preprocessing stage provenance trace."""

    stage: str
    operation: str
    parameters: dict[str, object]
    input_shape: tuple[int, int]
    output_shape: tuple[int, int]
    input_hash: str
    output_hash: str
    phospho_input_hash: str | None = None
    phospho_output_hash: str | None = None
    dropped_row_ids: tuple[str, ...] = ()
    dropped_row_count: int = 0
    schema_version: int = PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3
    consumed_input_tables: tuple[TableFingerprint, ...] = ()
    produced_output_tables: tuple[TableFingerprint, ...] = ()
    backend: str | None = None
    random_seed: int | None = None
    determinism: DeterminismKind | str = DeterminismKind.DETERMINISTIC
    reproducibility_caveats: tuple[ReproducibilityCaveat, ...] = ()
    is_deterministic: bool = True
    imputed_cell_count: int = 0
    imputed_row_ids: tuple[str, ...] = ()
    notes: str | None = None
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    batch_correction_provenance: BatchCorrectionProvenance | None = None
    intensity_transformation_event: IntensityTransformationEvent | None = None
    quantitative_transition_evidence: QuantitativeTransitionEvidence | None = None
    quantitative_contract: QuantitativeOperationContract | None = None

    @property
    def input_rows(self) -> int:
        return int(self.input_shape[0])

    @property
    def output_rows(self) -> int:
        return int(self.output_shape[0])


class PreprocessingStage(Protocol):
    """Single internal preprocessing stage contract."""

    stage_key: str

    def run(self, state: PreprocessingState) -> PreprocessingStageResult: ...
