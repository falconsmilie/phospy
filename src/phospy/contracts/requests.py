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
from typing import Protocol

from phospy.contracts import dataset_build as _dataset_build
from phospy.contracts.configs import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    DifferentialAnalysisConfig,
    EnrichmentConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseReferenceDisplayAmbiguityPolicy,
    KinaseScoringConfig,
    KinaseSiteSequenceConflictPolicy,
    SignalomeConfig,
)
from phospy.contracts.dataset_build import DatasetBuildRequest, DatasetInput
from phospy.contracts.enrichment_identifier_sets import (
    EnrichmentDerivedQuantitativeSetProvenance,
    EnrichmentDerivedSetMissingValueRule,
    EnrichmentDerivedSetSourceResultKind,
    EnrichmentDerivedSetThresholdDirection,
    EnrichmentDerivedSetValueMeaning,
    EnrichmentDerivedSetValueScale,
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
)
from phospy.contracts.results.kinase import KinaseWorkflowResult
from phospy.provenance.models import InputIntensityScaleEvidence
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.contrast_helpers import (
    all_pairwise_contrasts,
    contrasts_vs_control,
)
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
)
from phospy.science.references.kinase_library_models import KinaseLibraryResource
from phospy.science.references.models import ReferenceBundle, ReferencePreset

DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING = (
    _dataset_build.DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
)
DATASET_MULTI_SITE_POLICY_REJECT = _dataset_build.DATASET_MULTI_SITE_POLICY_REJECT
DATASET_MULTI_SITE_POLICY_SPLIT = _dataset_build.DATASET_MULTI_SITE_POLICY_SPLIT
DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE = (
    _dataset_build.DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE
)
DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED = (
    _dataset_build.DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED
)


@dataclass(frozen=True, slots=True, eq=False)
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

    Python equality and hashing are identity-based. The importer owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

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

    def run(self, request: PhosphositeImportRequest) -> object: ...


@dataclass(frozen=True, slots=True, eq=False)
class EnrichmentWorkflowRequest:
    """Request for native enrichment.

    Construction stores the payload only. Identifier semantics, collection and
    background compatibility, input-source selection, input emptiness, and
    config type are validated by ``EnrichmentWorkflow.run(...)`` before
    interpretation or execution.

    Python equality and hashing are identity-based. The workflow owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

    identifier_column: str
    identifier_kind: EnrichmentIdentifierKind
    set_collection: EnrichmentSetCollection
    background_universe: Sequence[str]
    config: EnrichmentConfig = field(default_factory=EnrichmentConfig)
    input_table: DatasetInput | None = None
    selected_identifiers: Sequence[str] | None = None
    selected_identifier_provenance: EnrichmentIdentifierSetProvenance | None = None
    background_identifier_provenance: EnrichmentIdentifierSetProvenance | None = None


@dataclass(frozen=True, slots=True, eq=False)
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

    Python equality and hashing are identity-based. The workflow owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

    dataset: AnalysisReadyPhosphoDataset
    references: ReferencePreset | ReferenceBundle = ReferencePreset.AUTO
    scoring_config: KinaseScoringConfig | None = None
    prediction_config: KinasePredictionConfig = field(
        default_factory=KinasePredictionConfig
    )
    activity_config: KinaseActivityConfig | None = None
    site_sequence_conflict_policy: KinaseSiteSequenceConflictPolicy = (
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE
    )
    reference_display_ambiguity_policy: KinaseReferenceDisplayAmbiguityPolicy = (
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR
    )
    kinase_library_resource: KinaseLibraryResource | None = None


@dataclass(frozen=True, slots=True, eq=False)
class SignalomeWorkflowRequest:
    """Request for the public signalome workflow.

    Construction stores the payload only. The upstream kinase result, signalome
    config, score/prediction matrix alignment, site identity, and protein
    grouping metadata are validated by ``SignalomeWorkflow.run(...)`` before
    interpretation or execution.

    Signalome execution requires explicit protein grouping metadata per
    interpreted site via ``dataset.site_metadata.protein_group_id``. The legacy
    ``dataset.site_metadata.protein_id`` column is accepted only as a migration
    alias.

    Python equality and hashing are identity-based. The workflow owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

    kinase_result: KinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)


@dataclass(frozen=True, slots=True, eq=False)
class DifferentialAnalysisRequest:
    """Request for the public differential workflow.

    Construction stores the payload only. Dataset eligibility, design/sample
    alignment, contrast validity, replicate requirements, and config coherence
    are validated by ``DifferentialAnalysisWorkflow.run(...)`` before
    interpretation or statistical execution.

    Python equality and hashing are identity-based. The workflow owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

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
    "DatasetInput",
    "PhosphositeImporter",
    "PhosphositeImportRequest",
    "SampleDesignRecord",
    "all_pairwise_contrasts",
    "contrasts_vs_control",
    "ExperimentalDesign",
    "FixedEffectCovariate",
    "Contrast",
    "ContinuousCovariate",
    "DifferentialAnalysisRequest",
    "EmpiricalBayesConfig",
    "EnrichmentDerivedQuantitativeSetProvenance",
    "EnrichmentDerivedSetMissingValueRule",
    "EnrichmentDerivedSetSourceResultKind",
    "EnrichmentDerivedSetThresholdDirection",
    "EnrichmentDerivedSetValueMeaning",
    "EnrichmentDerivedSetValueScale",
    "EnrichmentIdentifierKind",
    "EnrichmentIdentifierSetProvenance",
    "EnrichmentIdentifierSetSourceType",
    "EnrichmentSet",
    "EnrichmentSetCollection",
    "EnrichmentWorkflowRequest",
    "GeneSetCollection",
    "InputIntensityScaleEvidence",
    "KinaseWorkflowRequest",
    "PtmSetCollection",
    "SignalomeWorkflowRequest",
]
