"""Processing-state data models."""

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.datasets.preprocessing.policy_models import (
    MissingDataPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.errors.input import PhosPyInputError
from phospy.transformations.models import IntensityScaleState, QuantitativeMeaning

from .missing_data import MissingDataDiagnostics
from .total_protein import TotalProteinCorrectionDiagnostics


@dataclass(frozen=True, slots=True)
class MissingDataState:
    """Missing-data policy state at the analysis-ready dataset boundary."""

    policy: MissingDataPolicy
    min_observed_values: int | None
    complete_matrix: bool
    imputed: bool
    diagnostics: MissingDataDiagnostics | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "policy",
            MissingDataPolicy.parse(
                self.policy,
                field_name="dataset processing state missing_data.policy",
            ),
        )
        if self.diagnostics is None:
            return
        if isinstance(self.diagnostics, MissingDataDiagnostics):
            return
        normalized = MissingDataDiagnostics.from_payload(
            self.diagnostics,
            field_name="dataset processing state missing_data.diagnostics",
        )
        object.__setattr__(self, "diagnostics", normalized)


@dataclass(frozen=True, slots=True)
class NormalisationState:
    """Normalisation policy state at the analysis-ready dataset boundary."""

    policy: str


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionRowDiagnostic:
    """Durable per-row site-sequence resolution diagnostic record."""

    row_index: int
    row_id: str
    site_id: str | None
    status: str
    existing_site_sequence: str | None
    fasta_site_sequence: str | None
    resolved_site_sequence: str | None
    action: str
    reason: str | None
    conflict_policy: str | None
    resolver_version: str | None
    fasta_source_path: str | None
    fasta_sha256: str | None


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionState:
    """Dataset site-sequence FASTA-resolution state at preprocessing boundary."""

    configured: bool
    mode: str | None
    flank_size: int | None
    fasta_source_path: str | None
    fasta_source_label: str | None
    fasta_sha256: str | None
    resolver_version: str | None
    resolved_site_count: int
    unresolved_site_count: int
    unresolved_counts_by_reason: dict[str, int]
    filled_missing_count: int
    replaced_existing_count: int
    preserved_existing_count: int
    existing_sequence_conflict_count: int
    conflict_policy: str | None = None
    row_diagnostics: tuple[SiteSequenceResolutionRowDiagnostic, ...] = ()


@dataclass(frozen=True, slots=True)
class TotalProteinCorrectionState:
    """Total-protein correction state at the analysis-ready dataset boundary."""

    policy: TotalProteinCorrectionPolicy
    applied: bool
    formula: str | None = None
    requires_log_scale: bool | None = False
    input_scale: str | None = None
    output_scale: str | None = None
    quantitative_meaning: QuantitativeMeaning | None = None
    diagnostics: TotalProteinCorrectionDiagnostics | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "policy",
            TotalProteinCorrectionPolicy.parse(
                self.policy,
                field_name="dataset processing state total_protein_correction.policy",
            ),
        )
        quantitative_meaning = self.quantitative_meaning
        if quantitative_meaning is not None and not isinstance(
            quantitative_meaning, QuantitativeMeaning
        ):
            try:
                quantitative_meaning = QuantitativeMeaning(str(quantitative_meaning))
            except ValueError as exc:
                supported = ", ".join(member.value for member in QuantitativeMeaning)
                raise PhosPyInputError(
                    "dataset processing state total_protein_correction."
                    "quantitative_meaning must be one of: "
                    f"{supported}"
                ) from exc
            object.__setattr__(self, "quantitative_meaning", quantitative_meaning)
        if self.diagnostics is None:
            return
        if isinstance(self.diagnostics, TotalProteinCorrectionDiagnostics):
            return
        normalized = TotalProteinCorrectionDiagnostics.from_payload(
            self.diagnostics,
            field_name="dataset processing state total_protein_correction.diagnostics",
        )
        object.__setattr__(self, "diagnostics", normalized)


@dataclass(frozen=True, slots=True)
class SiteMatrixState:
    """Site-matrix construction state at the analysis-ready dataset boundary."""

    policy: str
    constructed: bool
    missing_data_policy: str
    minimum_observed_values: int | None
    duplicate_site_policy: str


@dataclass(frozen=True, slots=True)
class ComparisonState:
    """Comparison-building state at the analysis-ready dataset boundary."""

    policy: str
    sample_group_column: str
    pairs: tuple[tuple[str, str], ...] | None


@dataclass(frozen=True, slots=True)
class RuvReadinessState:
    """RUV-compatible preprocessing readiness reporting state."""

    enabled: bool
    ready: bool
    reasons: tuple[str, ...]
    control_feature_column: str
    replicate_group_column: str
    batch_column: str | None
    control_feature_count: int
    replicate_group_count: int
    batch_count: int | None
    requires_complete_matrix: bool
    matrix_complete: bool
    imputation_method_id: str | None
    missingness_mask_preserved: bool


def default_ruv_readiness_state() -> RuvReadinessState:
    """Return the default disabled readiness state."""

    return RuvReadinessState(
        enabled=False,
        ready=False,
        reasons=("not configured",),
        control_feature_column="is_control_feature",
        replicate_group_column="replicate_group",
        batch_column="batch",
        control_feature_count=0,
        replicate_group_count=0,
        batch_count=0,
        requires_complete_matrix=True,
        matrix_complete=False,
        imputation_method_id=None,
        missingness_mask_preserved=False,
    )


@dataclass(frozen=True, slots=True)
class DatasetProcessingState:
    """Compact summary of preprocessing state at the analysis-ready boundary."""

    intensity_scale: IntensityScaleState
    site_sequence_resolution: SiteSequenceResolutionState
    missing_data: MissingDataState
    normalisation: NormalisationState
    total_protein_correction: TotalProteinCorrectionState
    site_matrix: SiteMatrixState
    comparisons: ComparisonState
    ruv_readiness: RuvReadinessState = field(
        default_factory=default_ruv_readiness_state
    )
