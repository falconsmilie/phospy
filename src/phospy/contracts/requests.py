"""Public request models.

Request dataclasses are lightweight command payloads. They intentionally do not
perform scientific validation during construction; the relevant dataset builder
or workflow validator enforces request compatibility before interpretation and
execution. Config dataclasses may still self-validate local configuration
invariants.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

from phospy.contracts.configs import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    DatasetPreprocessingConfig,
    DifferentialAnalysisConfig,
    EnrichmentConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseReferenceDisplayAmbiguityPolicy,
    KinaseScoringConfig,
    KinaseSiteSequenceConflictPolicy,
    SignalomeConfig,
)
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.builders.contracts import DatasetInput
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.models import (
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    ExperimentalDesign,
    FixedEffectCovariate,
    SampleDesignRecord,
)
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    EmpiricalBayesConfig,
)
from phospy.science.enrichment.models import (
    EnrichmentIdentifierKind,
    EnrichmentSet,
    EnrichmentSetCollection,
    GeneSetCollection,
    PtmSetCollection,
    _normalise_identifier_sequence,
    _require_collection_matches_identifier_kind,
    _require_identifier_kind,
    _require_non_empty_string,
)
from phospy.science.evidence import dataset_resolution as _dataset_resolution
from phospy.science.references.kinase_library import KinaseLibraryResource
from phospy.science.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)

if TYPE_CHECKING:
    from phospy.contracts.results import KinaseWorkflowResult, PhosphositeImportResult

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
class PhosphositeImportRequest:
    """Request for translating an upstream phosphosite table into builder inputs.

    Construction stores the payload only. Import-source type checks, required
    column mapping validation, localisation-confidence parsing, and candidate
    table construction are owned by ``PhosphositeImporter.run(...)``. Dataset
    scientific readiness remains owned by ``AnalysisReadyDatasetBuilder``.

    ``sample_intensity_columns`` is explicit by design. It accepts either a
    sequence of source column names, which keeps sample IDs equal to source
    column names, or a mapping of ``source_column -> sample_id``. Importers do
    not infer sample groups, contrasts, batches, or differential designs from
    these names.
    """

    source: DatasetInput
    sample_intensity_columns: Mapping[str, str] | Sequence[str]
    gene_symbol_column: str = "gene_symbol"
    site_column: str = "site"
    row_id_column: str | None = None
    protein_id_column: str | None = None
    protein_accession_column: str | None = None
    protein_identifier_column: str | None = None
    protein_namespace_column: str | None = None
    organism_column: str | None = None
    isoform_id_column: str | None = None
    site_sequence_column: str | None = None
    display_id_column: str | None = None
    site_key_column: str | None = None
    localisation_confidence_column: str | None = None
    localisation_confidence_scale: str = "probability"
    peptide_row_id_column: str | None = None
    unique_feature_id_column: str | None = None
    peptide_sequence_column: str | None = None
    modified_peptide_sequence_column: str | None = None
    peptide_site_string_column: str | None = None
    peptide_site_id_column: str | None = None
    source_name: str = "phosphosite_import"


class PhosphositeImporter(Protocol):
    """Protocol for upstream phosphosite importers."""

    def run(self, request: PhosphositeImportRequest) -> PhosphositeImportResult: ...


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
class EnrichmentWorkflowRequest:
    """Request for native enrichment.

    Construction records user intent and validates only local contract
    invariants: one input source, explicit identifier semantics, explicit
    gene/PTM-set collection, explicit non-empty background universe, and config
    type. It does not calculate enrichment statistics, load online resources,
    or infer a background universe.
    """

    identifier_column: str
    identifier_kind: EnrichmentIdentifierKind
    set_collection: EnrichmentSetCollection
    background_universe: Sequence[str]
    config: EnrichmentConfig = field(default_factory=EnrichmentConfig)
    input_table: DatasetInput | None = None
    selected_identifiers: Sequence[str] | None = None

    def __post_init__(self) -> None:
        identifier_column = _require_non_empty_string(
            self.identifier_column,
            field_name="enrichment_request.identifier_column",
        )
        identifier_kind = _require_identifier_kind(
            self.identifier_kind,
            field_name="enrichment_request.identifier_kind",
        )
        set_collection = _require_collection_matches_identifier_kind(
            self.set_collection,
            identifier_kind=identifier_kind,
            field_name="enrichment_request.set_collection",
        )
        background_universe = _normalise_identifier_sequence(
            self.background_universe,
            field_name="enrichment_request.background_universe",
        )
        if not isinstance(self.config, EnrichmentConfig):
            raise WorkflowValidationError(
                "enrichment_request.config must be EnrichmentConfig"
            )
        has_input_table = self.input_table is not None
        has_selected_identifiers = self.selected_identifiers is not None
        if has_input_table == has_selected_identifiers:
            raise WorkflowValidationError(
                "enrichment_request requires exactly one of input_table or "
                "selected_identifiers"
            )
        selected_identifiers = (
            None
            if self.selected_identifiers is None
            else _normalise_identifier_sequence(
                self.selected_identifiers,
                field_name="enrichment_request.selected_identifiers",
            )
        )
        object.__setattr__(self, "identifier_column", identifier_column)
        object.__setattr__(self, "identifier_kind", identifier_kind)
        object.__setattr__(self, "set_collection", set_collection)
        object.__setattr__(self, "background_universe", background_universe)
        object.__setattr__(self, "selected_identifiers", selected_identifiers)


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
    kinase_library_resource: KinaseLibraryResource | None = None


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
    "BatchCovariate",
    "ContrastMatrix",
    "CategoricalCovariate",
    "DesignMatrix",
    "DatasetBuildRequest",
    "PhosphositeImporter",
    "PhosphositeImportRequest",
    "SampleDesignRecord",
    "ExperimentalDesign",
    "FixedEffectCovariate",
    "Contrast",
    "ContinuousCovariate",
    "DifferentialAnalysisRequest",
    "EmpiricalBayesConfig",
    "EnrichmentIdentifierKind",
    "EnrichmentSet",
    "EnrichmentSetCollection",
    "EnrichmentWorkflowRequest",
    "GeneSetCollection",
    "KinaseWorkflowRequest",
    "PtmSetCollection",
    "SignalomeWorkflowRequest",
]
