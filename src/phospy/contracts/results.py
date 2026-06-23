"""Public result models."""
# pyright: reportMissingTypeStubs=false, reportUnnecessaryIsInstance=false
# pandas has no bundled stubs here; runtime boundary guards are intentionally retained.

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, TypeAlias

import pandas as pd

from phospy.contracts.configs import EnrichmentConfig
from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import WorkflowValidationError
from phospy.frames.ownership import (
    export_dataframe,
    export_optional_dataframe,
    own_dataframe,
    own_optional_dataframe,
)
from phospy.provenance.models import RunProvenance
from phospy.science.activities.models import (
    ActivityMethodDiagnostics,
    KinaseActivityResult,
    KseaZScoreActivityDiagnostics,
    SsgseaSubstrateEnrichmentActivityDiagnostics,
    WeightedSubstrateActivityDiagnostics,
)
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwareMappingDiagnostics,
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)
from phospy.science.differential.models import (
    DifferentialAnalysisResult,
    DifferentialContrastDefinition,
    DifferentialDesignMatrixSummary,
    DifferentialEmpiricalBayesProvenance,
    DifferentialFixedEffectCovariateProvenance,
    DifferentialMissingValuePolicyProvenance,
    DifferentialPolicyProvenance,
    DifferentialReplicatePolicyProvenance,
    DifferentialStatisticalTestingProvenance,
    DifferentialTechnicalReplicateGroup,
    DifferentialUnsupportedDesignPolicyProvenance,
)
from phospy.science.enrichment.models import (
    EnrichmentIdentifierKind,
    EnrichmentResultRecord,
    EnrichmentSetCollection,
    _normalise_identifier_sequence,
    _require_collection_matches_identifier_kind,
    _require_identifier_kind,
)
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.references.models import ReferenceBundle
from phospy.science.signalomes.constants import (
    DISPLAY_ID_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_SITE,
    SITE_KEY_COLUMN,
)
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAlignmentDiagnostics,
    SignalomeAssignments,
    SignalomeModules,
    SignalomeModuleSelectionDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
    default_signalome_alignment_diagnostics,
    default_signalome_module_selection_diagnostics,
    default_signalome_score_preconditioning_diagnostics,
)
from phospy.science.sites.validation import require_site_key_series
from phospy.tables.kinase import KinaseSubstrateContributionTable
from phospy.tables.signalome import (
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
)

if TYPE_CHECKING:
    from phospy.contracts.requests import DatasetBuildRequest


ImporterQualityStatus: TypeAlias = Literal[
    "reported",
    "not_reported",
    "not_applicable",
]

IMPORTER_QUALITY_STATUS_REPORTED = "reported"
IMPORTER_QUALITY_STATUS_NOT_REPORTED = "not_reported"
IMPORTER_QUALITY_STATUS_NOT_APPLICABLE = "not_applicable"


@dataclass(frozen=True, slots=True)
class ImporterDetectedIntensityColumn:
    """Source intensity column mapped to one PhosPy sample ID."""

    source_column: str
    sample_id: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "source_column",
            _validate_non_empty_quality_text(
                self.source_column,
                field_name="importer_quality.detected_intensity_column.source_column",
            ),
        )
        object.__setattr__(
            self,
            "sample_id",
            _validate_non_empty_quality_text(
                self.sample_id,
                field_name="importer_quality.detected_intensity_column.sample_id",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible intensity-column payload."""

        return {
            "source_column": self.source_column,
            "sample_id": self.sample_id,
        }


@dataclass(frozen=True, slots=True)
class ImporterQualityCount:
    """Optional count with explicit reporting status."""

    status: ImporterQualityStatus = IMPORTER_QUALITY_STATUS_NOT_REPORTED
    count: int | None = None
    source_column: str | None = None
    policy: str | None = None
    reason: str | None = None

    def __post_init__(self) -> None:
        status = _validate_quality_status(
            self.status,
            field_name="importer_quality.count.status",
        )
        count = _validate_optional_non_negative_quality_int(
            self.count,
            field_name="importer_quality.count.count",
        )
        if status == IMPORTER_QUALITY_STATUS_REPORTED and count is None:
            raise PhosPyInputError(
                "importer_quality.count.count must be provided when status='reported'"
            )
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "count", count)
        object.__setattr__(
            self,
            "source_column",
            _validate_optional_quality_text(
                self.source_column,
                field_name="importer_quality.count.source_column",
            ),
        )
        object.__setattr__(
            self,
            "policy",
            _validate_optional_quality_text(
                self.policy,
                field_name="importer_quality.count.policy",
            ),
        )
        object.__setattr__(
            self,
            "reason",
            _validate_optional_quality_text(
                self.reason,
                field_name="importer_quality.count.reason",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible count payload."""

        return {
            "status": self.status,
            "count": self.count,
            "source_column": self.source_column,
            "policy": self.policy,
            "reason": self.reason,
        }


@dataclass(frozen=True, slots=True)
class ImporterMissingIntensitySummary:
    """Missing intensity counts observed during importer parsing."""

    status: ImporterQualityStatus = IMPORTER_QUALITY_STATUS_NOT_REPORTED
    total_missing_values: int | None = None
    rows_with_any_missing_intensity: int | None = None
    missing_values_by_sample_id: Mapping[str, int] = field(default_factory=dict)
    missing_values_by_source_column: Mapping[str, int] = field(default_factory=dict)
    reason: str | None = None

    def __post_init__(self) -> None:
        status = _validate_quality_status(
            self.status,
            field_name="importer_quality.missing_intensity.status",
        )
        total_missing_values = _validate_optional_non_negative_quality_int(
            self.total_missing_values,
            field_name="importer_quality.missing_intensity.total_missing_values",
        )
        rows_with_any_missing_intensity = _validate_optional_non_negative_quality_int(
            self.rows_with_any_missing_intensity,
            field_name=(
                "importer_quality.missing_intensity.rows_with_any_missing_intensity"
            ),
        )
        if status == IMPORTER_QUALITY_STATUS_REPORTED:
            if total_missing_values is None:
                raise PhosPyInputError(
                    "importer_quality.missing_intensity.total_missing_values must "
                    "be provided when status='reported'"
                )
            if rows_with_any_missing_intensity is None:
                raise PhosPyInputError(
                    "importer_quality.missing_intensity."
                    "rows_with_any_missing_intensity must be provided when "
                    "status='reported'"
                )
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "total_missing_values", total_missing_values)
        object.__setattr__(
            self,
            "rows_with_any_missing_intensity",
            rows_with_any_missing_intensity,
        )
        object.__setattr__(
            self,
            "missing_values_by_sample_id",
            _validate_quality_count_mapping(
                self.missing_values_by_sample_id,
                field_name="importer_quality.missing_intensity.by_sample_id",
            ),
        )
        object.__setattr__(
            self,
            "missing_values_by_source_column",
            _validate_quality_count_mapping(
                self.missing_values_by_source_column,
                field_name="importer_quality.missing_intensity.by_source_column",
            ),
        )
        object.__setattr__(
            self,
            "reason",
            _validate_optional_quality_text(
                self.reason,
                field_name="importer_quality.missing_intensity.reason",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible missing-intensity payload."""

        return {
            "status": self.status,
            "total_missing_values": self.total_missing_values,
            "rows_with_any_missing_intensity": self.rows_with_any_missing_intensity,
            "missing_values_by_sample_id": dict(self.missing_values_by_sample_id),
            "missing_values_by_source_column": dict(
                self.missing_values_by_source_column
            ),
            "reason": self.reason,
        }


@dataclass(frozen=True, slots=True)
class ImporterLocalisationConfidenceSummary:
    """Localisation-confidence parsing facts when an importer has them."""

    status: ImporterQualityStatus = IMPORTER_QUALITY_STATUS_NOT_REPORTED
    source_column: str | None = None
    output_column: str | None = None
    scale: str | None = None
    row_count: int | None = None
    missing_count: int | None = None
    invalid_count: int | None = None
    invalid_examples: tuple[str, ...] = ()
    reason: str | None = None

    def __post_init__(self) -> None:
        status = _validate_quality_status(
            self.status,
            field_name="importer_quality.localisation_confidence.status",
        )
        row_count = _validate_optional_non_negative_quality_int(
            self.row_count,
            field_name="importer_quality.localisation_confidence.row_count",
        )
        missing_count = _validate_optional_non_negative_quality_int(
            self.missing_count,
            field_name="importer_quality.localisation_confidence.missing_count",
        )
        invalid_count = _validate_optional_non_negative_quality_int(
            self.invalid_count,
            field_name="importer_quality.localisation_confidence.invalid_count",
        )
        if status == IMPORTER_QUALITY_STATUS_REPORTED:
            for field_name, value in (
                ("source_column", self.source_column),
                ("output_column", self.output_column),
                ("scale", self.scale),
            ):
                if not isinstance(value, str) or value.strip() == "":
                    raise PhosPyInputError(
                        "importer_quality.localisation_confidence."
                        f"{field_name} must be provided when status='reported'"
                    )
            for field_name, value in (
                ("row_count", row_count),
                ("missing_count", missing_count),
                ("invalid_count", invalid_count),
            ):
                if value is None:
                    raise PhosPyInputError(
                        "importer_quality.localisation_confidence."
                        f"{field_name} must be provided when status='reported'"
                    )
        object.__setattr__(self, "status", status)
        object.__setattr__(
            self,
            "source_column",
            _validate_optional_quality_text(
                self.source_column,
                field_name="importer_quality.localisation_confidence.source_column",
            ),
        )
        object.__setattr__(
            self,
            "output_column",
            _validate_optional_quality_text(
                self.output_column,
                field_name="importer_quality.localisation_confidence.output_column",
            ),
        )
        object.__setattr__(
            self,
            "scale",
            _validate_optional_quality_text(
                self.scale,
                field_name="importer_quality.localisation_confidence.scale",
            ),
        )
        object.__setattr__(self, "row_count", row_count)
        object.__setattr__(self, "missing_count", missing_count)
        object.__setattr__(self, "invalid_count", invalid_count)
        object.__setattr__(
            self,
            "invalid_examples",
            tuple(
                _validate_non_empty_quality_text(
                    value,
                    field_name=(
                        "importer_quality.localisation_confidence.invalid_examples"
                    ),
                )
                for value in self.invalid_examples
            ),
        )
        object.__setattr__(
            self,
            "reason",
            _validate_optional_quality_text(
                self.reason,
                field_name="importer_quality.localisation_confidence.reason",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible localisation payload."""

        return {
            "status": self.status,
            "source_column": self.source_column,
            "output_column": self.output_column,
            "scale": self.scale,
            "row_count": self.row_count,
            "missing_count": self.missing_count,
            "invalid_count": self.invalid_count,
            "invalid_examples": list(self.invalid_examples),
            "reason": self.reason,
        }


@dataclass(frozen=True, slots=True)
class ImporterFlaggedRowSummary:
    """Contaminant, reverse, and decoy counts when a format reports them."""

    contaminant: ImporterQualityCount = field(default_factory=ImporterQualityCount)
    reverse: ImporterQualityCount = field(default_factory=ImporterQualityCount)
    decoy: ImporterQualityCount = field(default_factory=ImporterQualityCount)

    def __post_init__(self) -> None:
        for field_name, value in (
            ("contaminant", self.contaminant),
            ("reverse", self.reverse),
            ("decoy", self.decoy),
        ):
            if not isinstance(value, ImporterQualityCount):
                raise PhosPyInputError(
                    f"importer_quality.flagged_rows.{field_name} must be "
                    "ImporterQualityCount"
                )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible flagged-row payload."""

        return {
            "contaminant": self.contaminant.to_payload(),
            "reverse": self.reverse.to_payload(),
            "decoy": self.decoy.to_payload(),
        }


@dataclass(frozen=True, slots=True)
class ImporterDuplicateKeySummary:
    """Duplicate site-key and display-key counts when importer input has them."""

    site_key: ImporterQualityCount = field(default_factory=ImporterQualityCount)
    display_key: ImporterQualityCount = field(default_factory=ImporterQualityCount)
    duplicate_site_candidate_rows: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.site_key, ImporterQualityCount):
            raise PhosPyInputError(
                "importer_quality.duplicate_keys.site_key must be ImporterQualityCount"
            )
        if not isinstance(self.display_key, ImporterQualityCount):
            raise PhosPyInputError(
                "importer_quality.duplicate_keys.display_key must be "
                "ImporterQualityCount"
            )
        object.__setattr__(
            self,
            "duplicate_site_candidate_rows",
            _validate_optional_non_negative_quality_int(
                self.duplicate_site_candidate_rows,
                field_name=(
                    "importer_quality.duplicate_keys.duplicate_site_candidate_rows"
                ),
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible duplicate-key payload."""

        return {
            "site_key": self.site_key.to_payload(),
            "display_key": self.display_key.to_payload(),
            "duplicate_site_candidate_rows": self.duplicate_site_candidate_rows,
        }


@dataclass(frozen=True, slots=True)
class ImporterQualityReport:
    """Input-quality facts reported by an importer before dataset construction.

    The report is informational. It does not validate analysis readiness, infer
    experimental design, assign conditions, or change imported values.
    """

    source_name: str = "phosphosite_import"
    row_count_status: ImporterQualityStatus = IMPORTER_QUALITY_STATUS_NOT_REPORTED
    rows_read: int | None = None
    rows_retained: int | None = None
    rows_dropped: int | None = None
    intensity_column_status: ImporterQualityStatus = (
        IMPORTER_QUALITY_STATUS_NOT_REPORTED
    )
    detected_intensity_columns: tuple[ImporterDetectedIntensityColumn, ...] = ()
    missing_intensity: ImporterMissingIntensitySummary = field(
        default_factory=ImporterMissingIntensitySummary
    )
    localisation_confidence: ImporterLocalisationConfidenceSummary = field(
        default_factory=ImporterLocalisationConfidenceSummary
    )
    flagged_rows: ImporterFlaggedRowSummary = field(
        default_factory=ImporterFlaggedRowSummary
    )
    duplicate_keys: ImporterDuplicateKeySummary = field(
        default_factory=ImporterDuplicateKeySummary
    )
    format_specific: Mapping[str, object] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        row_count_status = _validate_quality_status(
            self.row_count_status,
            field_name="importer_quality.row_count_status",
        )
        intensity_column_status = _validate_quality_status(
            self.intensity_column_status,
            field_name="importer_quality.intensity_column_status",
        )
        rows_read = _validate_optional_non_negative_quality_int(
            self.rows_read,
            field_name="importer_quality.rows_read",
        )
        rows_retained = _validate_optional_non_negative_quality_int(
            self.rows_retained,
            field_name="importer_quality.rows_retained",
        )
        rows_dropped = _validate_optional_non_negative_quality_int(
            self.rows_dropped,
            field_name="importer_quality.rows_dropped",
        )
        if row_count_status == IMPORTER_QUALITY_STATUS_REPORTED:
            for field_name, value in (
                ("rows_read", rows_read),
                ("rows_retained", rows_retained),
                ("rows_dropped", rows_dropped),
            ):
                if value is None:
                    raise PhosPyInputError(
                        f"importer_quality.{field_name} must be provided when "
                        "row_count_status='reported'"
                    )
        if not isinstance(self.missing_intensity, ImporterMissingIntensitySummary):
            raise PhosPyInputError(
                "importer_quality.missing_intensity must be "
                "ImporterMissingIntensitySummary"
            )
        if not isinstance(
            self.localisation_confidence,
            ImporterLocalisationConfidenceSummary,
        ):
            raise PhosPyInputError(
                "importer_quality.localisation_confidence must be "
                "ImporterLocalisationConfidenceSummary"
            )
        if not isinstance(self.flagged_rows, ImporterFlaggedRowSummary):
            raise PhosPyInputError(
                "importer_quality.flagged_rows must be ImporterFlaggedRowSummary"
            )
        if not isinstance(self.duplicate_keys, ImporterDuplicateKeySummary):
            raise PhosPyInputError(
                "importer_quality.duplicate_keys must be ImporterDuplicateKeySummary"
            )
        if not isinstance(self.format_specific, Mapping):
            raise PhosPyInputError("importer_quality.format_specific must be a mapping")
        detected_columns = tuple(self.detected_intensity_columns)
        for column in detected_columns:
            if not isinstance(column, ImporterDetectedIntensityColumn):
                raise PhosPyInputError(
                    "importer_quality.detected_intensity_columns must contain "
                    "ImporterDetectedIntensityColumn values"
                )
        object.__setattr__(
            self,
            "source_name",
            _validate_source_name(self.source_name),
        )
        object.__setattr__(self, "row_count_status", row_count_status)
        object.__setattr__(self, "rows_read", rows_read)
        object.__setattr__(self, "rows_retained", rows_retained)
        object.__setattr__(self, "rows_dropped", rows_dropped)
        object.__setattr__(self, "intensity_column_status", intensity_column_status)
        object.__setattr__(self, "detected_intensity_columns", detected_columns)
        object.__setattr__(self, "format_specific", dict(self.format_specific))
        object.__setattr__(
            self,
            "warnings",
            tuple(_validate_warning(value) for value in self.warnings),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible importer quality payload."""

        return {
            "source_name": self.source_name,
            "row_count_status": self.row_count_status,
            "rows_read": self.rows_read,
            "rows_retained": self.rows_retained,
            "rows_dropped": self.rows_dropped,
            "intensity_column_status": self.intensity_column_status,
            "detected_intensity_columns": [
                column.to_payload() for column in self.detected_intensity_columns
            ],
            "missing_intensity": self.missing_intensity.to_payload(),
            "localisation_confidence": self.localisation_confidence.to_payload(),
            "flagged_rows": self.flagged_rows.to_payload(),
            "duplicate_keys": self.duplicate_keys.to_payload(),
            "format_specific": dict(self.format_specific),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True, slots=True, init=False)
class PhosphositeImportResult:
    """Candidate tables produced by an upstream phosphosite importer.

    Import results are not analysis-ready datasets. They expose normalized
    PhosPy input candidates plus diagnostics so callers can pass the candidate
    tables into ``AnalysisReadyDatasetBuilder`` without bypassing the builder's
    validation, preprocessing, site-key derivation, or peptide-evidence
    resolution responsibilities.
    """

    _phospho_matrix_candidate: pd.DataFrame
    _site_metadata_candidate: pd.DataFrame
    _peptide_evidence: pd.DataFrame | None
    _sample_column_mapping: dict[str, str]
    localisation_confidence_column: str | None
    warnings: tuple[str, ...]
    diagnostics: dict[str, object]
    source_name: str
    quality_report: ImporterQualityReport

    def __init__(
        self,
        *,
        phospho_matrix_candidate: pd.DataFrame,
        site_metadata_candidate: pd.DataFrame,
        peptide_evidence: pd.DataFrame | None = None,
        sample_column_mapping: dict[str, str],
        localisation_confidence_column: str | None = None,
        warnings: tuple[str, ...] = (),
        diagnostics: dict[str, object] | None = None,
        source_name: str = "phosphosite_import",
        quality_report: ImporterQualityReport | None = None,
        _assume_owned: bool = False,
    ) -> None:
        phospho = own_dataframe(
            phospho_matrix_candidate,
            field_name="phosphosite_import_result.phospho_matrix_candidate",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        site_metadata = own_dataframe(
            site_metadata_candidate,
            field_name="phosphosite_import_result.site_metadata_candidate",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        evidence = own_optional_dataframe(
            peptide_evidence,
            field_name="phosphosite_import_result.peptide_evidence",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        mapping = _validate_sample_column_mapping(sample_column_mapping)
        if localisation_confidence_column is not None and not isinstance(
            localisation_confidence_column,
            str,
        ):
            raise PhosPyInputError(
                "phosphosite_import_result.localisation_confidence_column must be "
                "a string or None"
            )
        if (
            isinstance(localisation_confidence_column, str)
            and localisation_confidence_column.strip() == ""
        ):
            raise PhosPyInputError(
                "phosphosite_import_result.localisation_confidence_column must be "
                "non-empty when provided"
            )
        warning_values = tuple(_validate_warning(value) for value in warnings)
        if diagnostics is not None and not isinstance(diagnostics, dict):
            raise PhosPyInputError(
                "phosphosite_import_result.diagnostics must be a dict or None"
            )
        source_name_value = _validate_source_name(source_name)
        if quality_report is None:
            quality_report_value = ImporterQualityReport(
                source_name=source_name_value,
                warnings=warning_values,
            )
        elif isinstance(quality_report, ImporterQualityReport):
            quality_report_value = quality_report
        else:
            raise PhosPyInputError(
                "phosphosite_import_result.quality_report must be "
                "ImporterQualityReport or None"
            )

        object.__setattr__(self, "_phospho_matrix_candidate", phospho)
        object.__setattr__(self, "_site_metadata_candidate", site_metadata)
        object.__setattr__(self, "_peptide_evidence", evidence)
        object.__setattr__(self, "_sample_column_mapping", mapping)
        object.__setattr__(
            self,
            "localisation_confidence_column",
            localisation_confidence_column,
        )
        object.__setattr__(self, "warnings", warning_values)
        object.__setattr__(self, "diagnostics", dict(diagnostics or {}))
        object.__setattr__(self, "source_name", source_name_value)
        object.__setattr__(self, "quality_report", quality_report_value)

    @property
    def phospho_matrix_candidate(self) -> pd.DataFrame:
        """Return a defensive snapshot of the site-by-sample matrix candidate."""

        return export_dataframe(self._phospho_matrix_candidate)

    @property
    def site_metadata_candidate(self) -> pd.DataFrame:
        """Return a defensive snapshot of the site metadata candidate."""

        return export_dataframe(self._site_metadata_candidate)

    @property
    def peptide_evidence(self) -> pd.DataFrame | None:
        """Return a defensive snapshot of optional peptide evidence."""

        return export_optional_dataframe(self._peptide_evidence)

    @property
    def sample_column_mapping(self) -> dict[str, str]:
        """Return a defensive ``source_column -> sample_id`` mapping snapshot."""

        return dict(self._sample_column_mapping)

    @property
    def peptide_evidence_sample_intensity_columns(self) -> tuple[str, ...]:
        """Return PhosPy sample IDs used as peptide-evidence intensity columns."""

        return tuple(self._sample_column_mapping.values())

    def to_dataset_build_request(
        self,
        *,
        site_resolution_mode: str = "site_level_resolved",
        multi_site_policy: str | None = None,
        sample_metadata: object | None = None,
        total: object | None = None,
        organism: object | None = None,
        preprocessing_config: object | None = None,
        allow_opaque_site_values: bool = False,
        input_intensity_scale: object | None = None,
        quantitative_meaning: object | None = None,
    ) -> DatasetBuildRequest:
        """Create a ``DatasetBuildRequest`` from importer candidates.

        This method intentionally returns a builder request rather than a
        dataset. The builder still owns analysis-ready validation,
        preprocessing, site identity derivation, and peptide-evidence
        resolution.
        """

        from phospy.contracts.requests import DatasetBuildRequest
        from phospy.science.evidence.dataset_resolution import (
            DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
        )

        common_kwargs = {
            "sample_metadata": sample_metadata,
            "total": total,
            "organism": organism,
            "allow_opaque_site_values": allow_opaque_site_values,
            "input_intensity_scale": input_intensity_scale,
            "quantitative_meaning": quantitative_meaning,
        }
        if preprocessing_config is not None:
            common_kwargs["preprocessing_config"] = preprocessing_config

        if site_resolution_mode == DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED:
            if multi_site_policy is not None:
                raise PhosPyInputError(
                    "phosphosite import result multi_site_policy is only valid for "
                    "site_resolution_mode='peptide_evidence'"
                )
            return DatasetBuildRequest(
                phospho=self.phospho_matrix_candidate,
                site_metadata=self.site_metadata_candidate,
                site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
                **common_kwargs,
            )

        if site_resolution_mode == DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE:
            if self._peptide_evidence is None:
                raise PhosPyInputError(
                    "phosphosite import result has no peptide_evidence candidate"
                )
            if multi_site_policy is None:
                raise PhosPyInputError(
                    "phosphosite import result peptide_evidence handoff requires "
                    "multi_site_policy"
                )
            return DatasetBuildRequest(
                site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
                peptide_evidence=self.peptide_evidence,
                peptide_evidence_sample_intensity_columns=(
                    self.peptide_evidence_sample_intensity_columns
                ),
                multi_site_policy=multi_site_policy,
                **common_kwargs,
            )

        raise PhosPyInputError(
            "phosphosite import result site_resolution_mode must be one of: "
            "'site_level_resolved', 'peptide_evidence'"
        )


@dataclass(frozen=True, slots=True)
class EnrichmentWorkflowResult:
    """Top-level native enrichment result container.

    The result contract stores an explicit enrichment result shape. Direct
    construction validates only local container consistency; no enrichment
    statistics are calculated here.

    ``table``, ``result_table``, and ``to_dataframe()`` return in-memory
    defensive snapshots only. Exporting, formatting, plotting, and report
    generation belong to IO or presentation adapters.
    """

    identifier_kind: EnrichmentIdentifierKind
    set_collection: EnrichmentSetCollection
    config: EnrichmentConfig
    records: tuple[EnrichmentResultRecord, ...] = ()
    unmatched_identifiers: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    method_metadata: Mapping[str, object] = field(default_factory=dict)
    background_summary: Mapping[str, object] = field(default_factory=dict)
    set_collection_summary: Mapping[str, object] = field(default_factory=dict)
    provenance: RunProvenance | None = None
    _result_table: pd.DataFrame = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        identifier_kind = _require_identifier_kind(
            self.identifier_kind,
            field_name="enrichment_result.identifier_kind",
        )
        set_collection = _require_collection_matches_identifier_kind(
            self.set_collection,
            identifier_kind=identifier_kind,
            field_name="enrichment_result.set_collection",
        )
        if not isinstance(self.config, EnrichmentConfig):
            raise WorkflowValidationError(
                "enrichment_result.config must be EnrichmentConfig"
            )
        records = tuple(self.records)
        for record in records:
            if not isinstance(record, EnrichmentResultRecord):
                raise WorkflowValidationError(
                    "enrichment_result.records must contain "
                    "EnrichmentResultRecord values"
                )
            if record.identifier_kind != identifier_kind:
                raise WorkflowValidationError(
                    "enrichment_result.records identifier_kind values must match "
                    "enrichment_result.identifier_kind"
                )
            if record.collection_kind != set_collection.collection_kind:
                raise WorkflowValidationError(
                    "enrichment_result.records collection_kind values must match "
                    "enrichment_result.set_collection"
                )
        unmatched_identifiers = _normalise_identifier_sequence(
            self.unmatched_identifiers,
            field_name="enrichment_result.unmatched_identifiers",
            allow_empty=True,
        )
        warnings = tuple(_validate_enrichment_warning(value) for value in self.warnings)
        if not isinstance(self.diagnostics, Mapping):
            raise WorkflowValidationError(
                "enrichment_result.diagnostics must be a mapping"
            )
        if not isinstance(self.method_metadata, Mapping):
            raise WorkflowValidationError(
                "enrichment_result.method_metadata must be a mapping"
            )
        if not isinstance(self.background_summary, Mapping):
            raise WorkflowValidationError(
                "enrichment_result.background_summary must be a mapping"
            )
        if not isinstance(self.set_collection_summary, Mapping):
            raise WorkflowValidationError(
                "enrichment_result.set_collection_summary must be a mapping"
            )
        if self.provenance is not None and not isinstance(
            self.provenance, RunProvenance
        ):
            raise WorkflowValidationError(
                "enrichment_result.provenance must be RunProvenance or None"
            )
        result_table = own_dataframe(
            _enrichment_records_to_dataframe(records),
            field_name="enrichment_result.table",
            error_type=WorkflowValidationError,
            assume_owned=True,
        )
        object.__setattr__(self, "identifier_kind", identifier_kind)
        object.__setattr__(self, "set_collection", set_collection)
        object.__setattr__(self, "records", records)
        object.__setattr__(self, "unmatched_identifiers", unmatched_identifiers)
        object.__setattr__(self, "warnings", warnings)
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))
        object.__setattr__(self, "method_metadata", dict(self.method_metadata))
        object.__setattr__(self, "background_summary", dict(self.background_summary))
        object.__setattr__(
            self,
            "set_collection_summary",
            dict(self.set_collection_summary),
        )
        object.__setattr__(self, "_result_table", result_table)

    @property
    def table(self) -> pd.DataFrame:
        """Return a defensive snapshot of the enrichment result table."""

        return export_dataframe(self._result_table)

    @property
    def result_table(self) -> pd.DataFrame:
        """Return a defensive snapshot of the enrichment result table."""

        return export_dataframe(self._result_table)

    def to_dataframe(self) -> pd.DataFrame:
        """Return a defensive snapshot of the enrichment result table."""

        return export_dataframe(self._result_table)


def _enrichment_records_to_dataframe(
    records: tuple[EnrichmentResultRecord, ...],
) -> pd.DataFrame:
    columns = [
        "term_id",
        "term_name",
        "collection_kind",
        "identifier_kind",
        "input_overlap_count",
        "background_overlap_count",
        "set_size",
        "overlap_identifiers",
        "p_value",
        "adjusted_p_value",
        "correction_method",
        "enrichment_ratio",
    ]
    rows = [
        {
            "term_id": record.term_id,
            "term_name": record.term_name,
            "collection_kind": record.collection_kind,
            "identifier_kind": record.identifier_kind,
            "input_overlap_count": record.input_overlap_count,
            "background_overlap_count": record.background_overlap_count,
            "set_size": record.set_size,
            "overlap_identifiers": record.overlap_identifiers,
            "p_value": record.p_value,
            "adjusted_p_value": record.adjusted_p_value,
            "correction_method": record.correction_method,
            "enrichment_ratio": record.enrichment_ratio,
        }
        for record in records
    ]
    return pd.DataFrame(rows, columns=columns)


@dataclass(frozen=True, slots=True)
class KinaseWorkflowPreprocessingAttritionSummary:
    """Preprocessing-owned site attrition counters composed into kinase results."""

    input_rows: int
    rows_removed_during_preprocessing: int
    rows_removed_invalid_or_missing_site_identifiers: int
    duplicate_sites_merged_or_resolved: int
    output_rows: int
    sequence_complete_sites: int | None = None


@dataclass(frozen=True, slots=True)
class KinaseWorkflowScoringAttritionSummary:
    """Kinase workflow site-eligibility counters after preprocessing."""

    rows_removed_invalid_or_missing_site_identifiers: int
    final_quantitative_sites_entering_scoring: int
    sites_with_valid_site_sequence: int
    sites_without_usable_site_sequence: int
    sites_eligible_for_motif_scoring: int
    sites_with_kinase_substrate_reference_profile_evidence: int
    sites_contributing_to_final_fused_prediction_scoring_output: int
    sites_contributing_to_activity_scoring: int | None


@dataclass(frozen=True, slots=True)
class KinaseWorkflowSiteAttritionSummary:
    """Compact, user-facing kinase site attrition summary."""

    preprocessing: KinaseWorkflowPreprocessingAttritionSummary
    scoring: KinaseWorkflowScoringAttritionSummary


@dataclass(frozen=True, slots=True)
class KinaseEligibilityReport:
    """Compact, user-facing kinase workflow eligibility counters.

    `eligible_kinases` counts kinases whose projected, sequence-supported
    quantified substrates meet `KinaseScoringConfig.min_substrates`.
    `excluded_kinases_below_min_substrates` counts overlapping kinases with too
    few usable substrates. The report is count-based; it does not add per-kinase
    weak-support flags to scoring result tables.
    """

    total_dataset_sites: int
    sequence_complete_sites: int
    localisation_eligible_sites: int | None
    reference_overlap_sites: int
    excluded_no_reference_match: int
    excluded_low_localisation: int | None
    eligible_kinases: int
    excluded_kinases_below_min_substrates: int


@dataclass(frozen=True, slots=True, init=False)
class KinaseWorkflowResult:
    """Top-level public kinase workflow result.

    This is a workflow-owned container, not a direct user-construction
    validator. Workflow execution is the supported construction path for
    scientific coherence across `dataset`, `references`, scoring, prediction,
    optional activity, eligibility, attrition, and provenance. The nested stage
    result objects own their public table schemas; this container keeps the
    workflow-assembled objects together without re-running workflow validation.
    """

    dataset: AnalysisReadyPhosphoDataset
    references: ReferenceBundle
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult
    eligibility_report: KinaseEligibilityReport | None = None
    site_attrition_summary: KinaseWorkflowSiteAttritionSummary | None = None
    activity_result: KinaseActivityResult | None = None
    provenance: RunProvenance | None = None
    _substrate_contributions: pd.DataFrame | None = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __init__(
        self,
        dataset: AnalysisReadyPhosphoDataset,
        references: ReferenceBundle,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult,
        eligibility_report: KinaseEligibilityReport | None = None,
        site_attrition_summary: KinaseWorkflowSiteAttritionSummary | None = None,
        activity_result: KinaseActivityResult | None = None,
        provenance: RunProvenance | None = None,
        substrate_contributions: pd.DataFrame | None = None,
        *,
        _assume_owned: bool = False,
    ) -> None:
        object.__setattr__(self, "dataset", dataset)
        object.__setattr__(self, "references", references)
        object.__setattr__(self, "scoring_result", scoring_result)
        object.__setattr__(self, "prediction_result", prediction_result)
        object.__setattr__(self, "eligibility_report", eligibility_report)
        object.__setattr__(self, "site_attrition_summary", site_attrition_summary)
        object.__setattr__(self, "activity_result", activity_result)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(
            self,
            "_substrate_contributions",
            _own_optional_kinase_substrate_contributions(
                substrate_contributions,
                assume_owned=_assume_owned,
            ),
        )

    @property
    def input_dataset_preprocessing_report(self) -> DatasetPreprocessingReport | None:
        """Return preprocessing provenance of the input analysis-ready dataset."""

        return self.dataset.preprocessing_report

    @property
    def substrate_contributions(self) -> pd.DataFrame | None:
        """Return optional substrate-level contribution rows."""

        return export_optional_dataframe(self._substrate_contributions)

    def substrate_contributions_dataframe(self) -> pd.DataFrame | None:
        """Return optional substrate-level contribution rows."""

        return export_optional_dataframe(self._substrate_contributions)


def _own_optional_kinase_substrate_contributions(
    table: pd.DataFrame | None,
    *,
    assume_owned: bool,
) -> pd.DataFrame | None:
    if table is None:
        return None
    return KinaseSubstrateContributionTable(
        frame=table,
        _assume_owned=assume_owned,
    ).frame


@dataclass(frozen=True, slots=True, init=False)
class SignalomeWorkflowResult:
    """Top-level public signalome workflow result.

    This is a workflow-owned result object. Direct construction is supported for
    bundle reconstruction and tests, but the signalome workflow is the
    recommended public construction path for scientific coherence across module
    assignment, module summary, kinase-network, and context sidecars.

    `expanded_signalome` is a flattened optional DataFrame contract populated by
    the supported signalome executor lane. It includes focal-kinase rows with
    linked-kinase metadata, regulated module IDs, and selected site-membership
    rows with stable `site_order`. `score_preconditioning_diagnostics` reports
    downstream-score row preconditioning counts and the active
    `SignalomeConfig.validation.score_preconditioning_policy`.
    `site_membership` and
    `protein_site_context` provide optional signalome provenance sidecars for
    site-level and protein-level phosphosite context.

    Provenance in this object describes owned internal state at creation time.
    Public export helpers return defensive snapshots; mutating exported
    DataFrames does not mutate this owning result. Constructor validation is
    intentionally limited to ownership, provenance type, and public table
    identity contracts. It does not run clustering, mapping, scoring, reference
    resolution, dataset repair, file export, plotting, or report formatting.
    """

    dataset: AnalysisReadyPhosphoDataset
    kinase_result: KinaseWorkflowResult
    module_assignments: SignalomeAssignments
    signalome_modules: SignalomeModules
    kinase_network: KinaseNetwork
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics = field(
        default_factory=default_signalome_module_selection_diagnostics
    )
    score_preconditioning_diagnostics: SignalomeScorePreconditioningDiagnostics = field(
        default_factory=default_signalome_score_preconditioning_diagnostics
    )
    alignment_diagnostics: SignalomeAlignmentDiagnostics = field(
        default_factory=default_signalome_alignment_diagnostics
    )
    provenance: RunProvenance | None = None
    _expanded_signalome: pd.DataFrame | None = field(init=False, repr=False)
    _site_membership: pd.DataFrame | None = field(init=False, repr=False)
    _protein_site_context: pd.DataFrame | None = field(init=False, repr=False)
    _init_payload: (
        tuple[pd.DataFrame | None, pd.DataFrame | None, pd.DataFrame | None, bool]
        | None
    ) = field(init=False, repr=False, default=None)

    def __init__(
        self,
        dataset: AnalysisReadyPhosphoDataset,
        kinase_result: KinaseWorkflowResult,
        module_assignments: SignalomeAssignments,
        signalome_modules: SignalomeModules,
        kinase_network: KinaseNetwork,
        module_selection_diagnostics: SignalomeModuleSelectionDiagnostics | None = None,
        score_preconditioning_diagnostics: SignalomeScorePreconditioningDiagnostics
        | None = None,
        alignment_diagnostics: SignalomeAlignmentDiagnostics | None = None,
        expanded_signalome: pd.DataFrame | None = None,
        site_membership: pd.DataFrame | None = None,
        protein_site_context: pd.DataFrame | None = None,
        provenance: RunProvenance | None = None,
        _assume_owned: bool = False,
    ) -> None:
        object.__setattr__(self, "dataset", dataset)
        object.__setattr__(self, "kinase_result", kinase_result)
        object.__setattr__(self, "module_assignments", module_assignments)
        object.__setattr__(self, "signalome_modules", signalome_modules)
        object.__setattr__(self, "kinase_network", kinase_network)
        object.__setattr__(
            self,
            "module_selection_diagnostics",
            (
                default_signalome_module_selection_diagnostics()
                if module_selection_diagnostics is None
                else module_selection_diagnostics
            ),
        )
        object.__setattr__(
            self,
            "score_preconditioning_diagnostics",
            (
                default_signalome_score_preconditioning_diagnostics()
                if score_preconditioning_diagnostics is None
                else score_preconditioning_diagnostics
            ),
        )
        object.__setattr__(
            self,
            "alignment_diagnostics",
            (
                default_signalome_alignment_diagnostics()
                if alignment_diagnostics is None
                else alignment_diagnostics
            ),
        )
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(
            self,
            "_init_payload",
            (
                expanded_signalome,
                site_membership,
                protein_site_context,
                _assume_owned,
            ),
        )
        self.__post_init__()

    def __post_init__(self) -> None:
        payload = self._init_payload
        if payload is None:
            raise WorkflowValidationError(
                "signalome_result internal initialization payload missing"
            )
        expanded_signalome, site_membership, protein_site_context, _assume_owned = (
            payload
        )
        expanded_signalome = own_optional_dataframe(
            expanded_signalome,
            field_name="signalome_result.expanded_signalome",
            error_type=WorkflowValidationError,
            assume_owned=_assume_owned,
        )
        site_membership = own_optional_dataframe(
            site_membership,
            field_name="signalome_result.site_membership",
            error_type=WorkflowValidationError,
            assume_owned=_assume_owned,
        )
        protein_site_context = own_optional_dataframe(
            protein_site_context,
            field_name="signalome_result.protein_site_context",
            error_type=WorkflowValidationError,
            assume_owned=_assume_owned,
        )
        _validate_signalome_result_site_level_identity(
            dataset=self.dataset,
            module_assignments=self.module_assignments,
            expanded_signalome=expanded_signalome,
            site_membership=site_membership,
        )
        if site_membership is not None:
            site_membership = SignalomeSiteContext(
                frame=site_membership,
                _assume_owned=True,
            ).frame
        if protein_site_context is not None:
            protein_site_context = SignalomeProteinSiteContext(
                frame=protein_site_context,
                _assume_owned=True,
            ).frame
        if self.provenance is not None and not isinstance(
            self.provenance, RunProvenance
        ):
            raise WorkflowValidationError(
                "signalome_result.provenance must be RunProvenance or None"
            )
        object.__setattr__(self, "_expanded_signalome", expanded_signalome)
        object.__setattr__(self, "_site_membership", site_membership)
        object.__setattr__(self, "_protein_site_context", protein_site_context)
        object.__setattr__(self, "_init_payload", None)

    @property
    def expanded_signalome(self) -> pd.DataFrame | None:
        """Return an expanded-signalome snapshot when available."""

        return export_optional_dataframe(self._expanded_signalome)

    @property
    def site_membership(self) -> pd.DataFrame | None:
        """Return a site-membership snapshot when available."""

        return export_optional_dataframe(self._site_membership)

    @property
    def protein_site_context(self) -> pd.DataFrame | None:
        """Return a protein-site context snapshot when available."""

        return export_optional_dataframe(self._protein_site_context)

    @property
    def input_dataset_preprocessing_report(self) -> DatasetPreprocessingReport | None:
        """Return preprocessing provenance of the input analysis-ready dataset."""

        return self.dataset.preprocessing_report

    @classmethod
    def _from_owned(
        cls,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        kinase_result: KinaseWorkflowResult,
        module_assignments: SignalomeAssignments,
        signalome_modules: SignalomeModules,
        kinase_network: KinaseNetwork,
        module_selection_diagnostics: SignalomeModuleSelectionDiagnostics | None = None,
        score_preconditioning_diagnostics: SignalomeScorePreconditioningDiagnostics
        | None = None,
        alignment_diagnostics: SignalomeAlignmentDiagnostics | None = None,
        expanded_signalome: pd.DataFrame | None = None,
        site_membership: pd.DataFrame | None = None,
        protein_site_context: pd.DataFrame | None = None,
        provenance: RunProvenance | None = None,
    ) -> SignalomeWorkflowResult:
        return cls(
            dataset=dataset,
            kinase_result=kinase_result,
            module_assignments=module_assignments,
            signalome_modules=signalome_modules,
            kinase_network=kinase_network,
            module_selection_diagnostics=module_selection_diagnostics,
            score_preconditioning_diagnostics=score_preconditioning_diagnostics,
            alignment_diagnostics=alignment_diagnostics,
            expanded_signalome=expanded_signalome,
            site_membership=site_membership,
            protein_site_context=protein_site_context,
            provenance=provenance,
            _assume_owned=True,
        )

    def to_dataframe(self) -> pd.DataFrame | None:
        """Return an expanded-signalome snapshot, not an export."""

        return export_optional_dataframe(self._expanded_signalome)

    def site_membership_dataframe(self) -> pd.DataFrame | None:
        """Return a site-membership snapshot, not an export."""

        return export_optional_dataframe(self._site_membership)

    def protein_site_context_dataframe(self) -> pd.DataFrame | None:
        """Return a protein-site context snapshot, not an export."""

        return export_optional_dataframe(self._protein_site_context)


def _validate_signalome_result_site_level_identity(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    module_assignments: SignalomeAssignments,
    expanded_signalome: pd.DataFrame | None,
    site_membership: pd.DataFrame | None,
) -> None:
    dataset_identity = _signalome_dataset_identity_lookup(dataset)
    _validate_site_level_signalome_rows(
        table=module_assignments.table,
        field_name="signalome_result.module_assignments.table",
        dataset_identity=dataset_identity,
    )
    _validate_expanded_signalome_identity(
        expanded_signalome,
        dataset_identity=dataset_identity,
    )
    _validate_site_level_signalome_rows(
        table=site_membership,
        field_name="signalome_result.site_membership",
        dataset_identity=dataset_identity,
    )


def _validate_expanded_signalome_identity(
    expanded_signalome: pd.DataFrame | None,
    *,
    dataset_identity: dict[str, str],
) -> None:
    missing = [
        column
        for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN)
        if expanded_signalome is not None and column not in expanded_signalome.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowValidationError(
            f"signalome_result.expanded_signalome is missing required columns: {joined}"
        )
    if expanded_signalome is None or expanded_signalome.empty:
        return
    site_rows = expanded_signalome
    if EXPANDED_SIGNALOME_ROW_KIND_COLUMN in expanded_signalome.columns:
        site_rows = expanded_signalome.loc[
            expanded_signalome.loc[:, EXPANDED_SIGNALOME_ROW_KIND_COLUMN].astype(str)
            == EXPANDED_SIGNALOME_ROW_KIND_SITE,
            :,
        ]
    _validate_site_level_signalome_rows(
        table=site_rows,
        field_name="signalome_result.expanded_signalome",
        dataset_identity=dataset_identity,
    )


def _validate_site_level_signalome_rows(
    *,
    table: pd.DataFrame | None,
    field_name: str,
    dataset_identity: dict[str, str],
) -> None:
    if table is None:
        return
    missing = [
        column for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN) if column not in table
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowValidationError(
            f"{field_name} is missing required columns: {joined}"
        )
    if table.empty:
        return
    for column_name in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN):
        invalid_count = int(
            sum(
                1
                for value in table.loc[:, column_name].tolist()
                if not isinstance(value, str) or value.strip() == ""
            )
        )
        if invalid_count:
            raise WorkflowValidationError(
                f"{field_name} site rows require non-empty "
                f"{column_name} values; invalid_count={invalid_count}"
            )
    site_keys = table.loc[:, SITE_KEY_COLUMN].astype(str)
    require_site_key_series(
        site_keys,
        field_name=f"{field_name}.{SITE_KEY_COLUMN}",
        error_type=WorkflowValidationError,
    )
    display_ids = table.loc[:, DISPLAY_ID_COLUMN].astype(str)
    missing_site_keys = [
        site_key
        for site_key in dict.fromkeys(site_keys.tolist())
        if site_key not in dataset_identity
    ]
    if missing_site_keys:
        preview = ", ".join(repr(value) for value in missing_site_keys[:5])
        suffix = "" if len(missing_site_keys) <= 5 else " ..."
        raise WorkflowValidationError(
            f"{field_name}.{SITE_KEY_COLUMN} values must align to "
            "signalome_result.dataset; missing_site_keys="
            f"{preview}{suffix}"
        )
    mismatches = [
        f"{site_key!r}: observed={display_id!r}, expected={dataset_identity[site_key]!r}"
        for site_key, display_id in zip(
            site_keys.tolist(),
            display_ids.tolist(),
            strict=True,
        )
        if dataset_identity[site_key] != display_id
    ]
    if mismatches:
        preview = "; ".join(mismatches[:5])
        suffix = "" if len(mismatches) <= 5 else " ..."
        raise WorkflowValidationError(
            f"{field_name}.{DISPLAY_ID_COLUMN} values must match "
            "signalome_result.dataset.site_metadata.display_id for each site_key; "
            f"mismatches={preview}{suffix}"
        )


def _signalome_dataset_identity_lookup(
    dataset: AnalysisReadyPhosphoDataset,
) -> dict[str, str]:
    site_metadata = DatasetInternalView(dataset).site_metadata
    missing = [
        column
        for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN)
        if column not in site_metadata.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowValidationError(
            "signalome_result.dataset.site_metadata is missing required columns: "
            f"{joined}"
        )
    site_keys = site_metadata.loc[:, SITE_KEY_COLUMN].astype(str)
    require_site_key_series(
        site_keys,
        field_name=f"signalome_result.dataset.site_metadata.{SITE_KEY_COLUMN}",
        error_type=WorkflowValidationError,
    )
    display_ids = site_metadata.loc[:, DISPLAY_ID_COLUMN]
    invalid_display_id_count = int(
        sum(
            1
            for value in display_ids.tolist()
            if not isinstance(value, str) or value.strip() == ""
        )
    )
    if invalid_display_id_count:
        raise WorkflowValidationError(
            "signalome_result.dataset.site_metadata.display_id values must be "
            "non-empty strings; "
            f"invalid_count={invalid_display_id_count}"
        )
    return {
        site_key: display_id
        for site_key, display_id in zip(
            site_keys.tolist(),
            display_ids.astype(str).tolist(),
            strict=True,
        )
    }


def _validate_quality_status(
    value: object,
    *,
    field_name: str,
) -> ImporterQualityStatus:
    if value in {
        IMPORTER_QUALITY_STATUS_REPORTED,
        IMPORTER_QUALITY_STATUS_NOT_REPORTED,
        IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    }:
        return value  # type: ignore[return-value]
    raise PhosPyInputError(
        f"{field_name} must be one of: 'reported', 'not_reported', 'not_applicable'"
    )


def _validate_non_empty_quality_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value.strip()


def _validate_optional_quality_text(
    value: object,
    *,
    field_name: str,
) -> str | None:
    if value is None:
        return None
    return _validate_non_empty_quality_text(value, field_name=field_name)


def _validate_optional_non_negative_quality_int(
    value: object,
    *,
    field_name: str,
) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be a non-negative integer or None")
    if value < 0:
        raise PhosPyInputError(f"{field_name} must be non-negative")
    return int(value)


def _validate_quality_count_mapping(
    value: Mapping[str, int],
    *,
    field_name: str,
) -> dict[str, int]:
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    counts: dict[str, int] = {}
    for key, count in value.items():
        normalized_key = _validate_non_empty_quality_text(
            key,
            field_name=f"{field_name}.key",
        )
        normalized_count = _validate_optional_non_negative_quality_int(
            count,
            field_name=f"{field_name}[{normalized_key!r}]",
        )
        if normalized_count is None:  # pragma: no cover - helper returns None only
            raise PhosPyInputError(f"{field_name}[{normalized_key!r}] must be an int")
        counts[normalized_key] = normalized_count
    return counts


def _validate_sample_column_mapping(value: dict[str, str]) -> dict[str, str]:
    if not isinstance(value, dict):
        raise PhosPyInputError(
            "phosphosite_import_result.sample_column_mapping must be a dict"
        )
    if not value:
        raise PhosPyInputError(
            "phosphosite_import_result.sample_column_mapping must not be empty"
        )
    normalized: dict[str, str] = {}
    for source_column, sample_id in value.items():
        if not isinstance(source_column, str) or source_column.strip() == "":
            raise PhosPyInputError(
                "phosphosite_import_result.sample_column_mapping source columns "
                "must be non-empty strings"
            )
        if not isinstance(sample_id, str) or sample_id.strip() == "":
            raise PhosPyInputError(
                "phosphosite_import_result.sample_column_mapping sample IDs must "
                "be non-empty strings"
            )
        normalized[source_column.strip()] = sample_id.strip()
    if len(set(normalized.values())) != len(normalized):
        raise PhosPyInputError(
            "phosphosite_import_result.sample_column_mapping sample IDs must be unique"
        )
    return normalized


def _validate_warning(value: object) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(
            "phosphosite_import_result.warnings must contain non-empty strings"
        )
    return value.strip()


def _validate_source_name(value: object) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(
            "phosphosite_import_result.source_name must be a non-empty string"
        )
    return value.strip()


def _validate_enrichment_warning(value: object) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise WorkflowValidationError(
            "enrichment_result.warnings must contain non-empty strings"
        )
    return value.strip()


__all__ = [
    "ActivityMethodDiagnostics",
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "DifferentialAnalysisResult",
    "DifferentialContrastDefinition",
    "DifferentialDesignMatrixSummary",
    "DifferentialEmpiricalBayesProvenance",
    "DifferentialFixedEffectCovariateProvenance",
    "DifferentialMissingValuePolicyProvenance",
    "DifferentialPolicyProvenance",
    "DifferentialReplicatePolicyProvenance",
    "DifferentialStatisticalTestingProvenance",
    "DifferentialTechnicalReplicateGroup",
    "DifferentialUnsupportedDesignPolicyProvenance",
    "EnrichmentResultRecord",
    "EnrichmentWorkflowResult",
    "IMPORTER_QUALITY_STATUS_NOT_APPLICABLE",
    "IMPORTER_QUALITY_STATUS_NOT_REPORTED",
    "IMPORTER_QUALITY_STATUS_REPORTED",
    "ImporterDetectedIntensityColumn",
    "ImporterDuplicateKeySummary",
    "ImporterFlaggedRowSummary",
    "ImporterLocalisationConfidenceSummary",
    "ImporterMissingIntensitySummary",
    "ImporterQualityCount",
    "ImporterQualityReport",
    "ImporterQualityStatus",
    "KinaseActivityResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowResult",
    "KseaZScoreActivityDiagnostics",
    "PhosphositeImportResult",
    "ProteinAwareMappingDiagnostics",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
    "ProteinAwareSiteEligibility",
    "SsgseaSubstrateEnrichmentActivityDiagnostics",
    "SignalomeWorkflowResult",
    "WeightedSubstrateActivityDiagnostics",
]
