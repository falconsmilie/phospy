"""Base public result contracts for importer handoff."""
# pyright: reportMissingTypeStubs=false, reportUnnecessaryIsInstance=false

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal, TypeAlias

import pandas as pd

from phospy.contracts.dataset_build import (
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    DatasetBuildRequest,
)
from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import (
    export_dataframe,
    export_optional_dataframe,
    own_dataframe,
    own_optional_dataframe,
)

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


__all__ = [
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
    "PhosphositeImportResult",
]
