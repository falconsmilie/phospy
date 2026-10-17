"""Durable protein-aware preparation contracts owned by Datasets."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import cast

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.comparison import dataframe_equals
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.frames.validation import (
    require_columns,
    require_exact_index_match,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.preparation_records.protein_aware import (
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS,
    PROTEIN_AWARE_MATCHED_PAIR_COLUMNS,
    PROTEIN_AWARE_SITE_ELIGIBILITY_COLUMNS,
    ProteinAwareMappingDiagnostics,
    ProteinAwarePreparationEligibility,
    ProteinAwareSampleAlignmentDiagnostics,
    ProteinAwareSiteEligibility,
    ProteinAwareTransformationStateDiagnostics,
)
from phospy.provenance.immutability import freeze_json_mapping_with_error_type
from phospy.run_lineage.models import RunProvenance
from phospy.run_lineage.serialization import to_payload as provenance_to_payload

PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION = 1


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ProteinAwarePreparationReport:
    """Protein-aware preparation audit report."""

    __hash__ = object.__hash__

    site_eligibility: tuple[ProteinAwareSiteEligibility, ...]
    mapping_diagnostics: ProteinAwareMappingDiagnostics
    sample_alignment: ProteinAwareSampleAlignmentDiagnostics
    transformation_state: ProteinAwareTransformationStateDiagnostics | None
    preparation_policy: str
    protein_mapping_policy: str
    policy_parameters: Mapping[str, object]
    provenance: RunProvenance | None
    schema_version: int
    _site_eligibility_table: pd.DataFrame = field(init=False, repr=False)

    def __init__(
        self,
        *,
        site_eligibility: Sequence[ProteinAwareSiteEligibility],
        mapping_diagnostics: ProteinAwareMappingDiagnostics | None = None,
        sample_alignment: ProteinAwareSampleAlignmentDiagnostics,
        transformation_state: ProteinAwareTransformationStateDiagnostics | None = None,
        preparation_policy: str = (
            DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS
        ),
        protein_mapping_policy: str = (
            DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS
        ),
        policy_parameters: Mapping[str, object] | None = None,
        provenance: RunProvenance | None = None,
        schema_version: int = PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION,
    ) -> None:
        resolved_site_eligibility = _normalize_site_eligibility(site_eligibility)
        _require_instance(
            sample_alignment,
            expected_type=ProteinAwareSampleAlignmentDiagnostics,
            error_message=(
                "protein_aware_preparation_report.sample_alignment must be "
                "ProteinAwareSampleAlignmentDiagnostics"
            ),
        )
        if transformation_state is not None:
            _require_instance(
                transformation_state,
                expected_type=ProteinAwareTransformationStateDiagnostics,
                error_message=(
                    "protein_aware_preparation_report.transformation_state must be "
                    "ProteinAwareTransformationStateDiagnostics or None"
                ),
            )
        resolved_mapping_diagnostics = (
            ProteinAwareMappingDiagnostics()
            if mapping_diagnostics is None
            else mapping_diagnostics
        )
        _require_instance(
            resolved_mapping_diagnostics,
            expected_type=ProteinAwareMappingDiagnostics,
            error_message=(
                "protein_aware_preparation_report.mapping_diagnostics must be "
                "ProteinAwareMappingDiagnostics"
            ),
        )
        if provenance is not None:
            _require_instance(
                provenance,
                expected_type=RunProvenance,
                error_message=(
                    "protein_aware_preparation_report.provenance must be "
                    "RunProvenance or None"
                ),
            )
        schema_version_value = cast(object, schema_version)
        if not isinstance(schema_version_value, int) or isinstance(
            schema_version_value, bool
        ):
            raise DatasetValidationError(
                "protein_aware_preparation_report.schema_version must be an int"
            )
        object.__setattr__(self, "site_eligibility", resolved_site_eligibility)
        object.__setattr__(self, "mapping_diagnostics", resolved_mapping_diagnostics)
        object.__setattr__(self, "sample_alignment", sample_alignment)
        object.__setattr__(self, "transformation_state", transformation_state)
        object.__setattr__(
            self,
            "preparation_policy",
            _normalize_policy(
                preparation_policy,
                field_name="protein_aware_preparation_report.preparation_policy",
                supported=DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
            ),
        )
        object.__setattr__(
            self,
            "protein_mapping_policy",
            _normalize_policy(
                protein_mapping_policy,
                field_name="protein_aware_preparation_report.protein_mapping_policy",
                supported=DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
            ),
        )
        object.__setattr__(
            self,
            "policy_parameters",
            freeze_json_mapping_with_error_type(
                policy_parameters or {},
                field_name="protein_aware_preparation_report.policy_parameters",
                error_type=DatasetValidationError,
            ),
        )
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "schema_version", int(schema_version))
        object.__setattr__(
            self,
            "_site_eligibility_table",
            _site_eligibility_table(resolved_site_eligibility),
        )

    @property
    def site_eligibility_table(self) -> pd.DataFrame:
        """Return the per-site eligibility table as a defensive snapshot."""

        return export_dataframe(self._site_eligibility_table)

    @property
    def sample_alignment_diagnostics(self) -> ProteinAwareSampleAlignmentDiagnostics:
        """Return typed sample-alignment diagnostics."""

        return self.sample_alignment

    @property
    def missing_protein_abundance_diagnostics(self) -> pd.DataFrame:
        """Return per-site missing protein-abundance diagnostics."""

        return self.mapping_diagnostics.missing_protein_abundance

    @property
    def ambiguous_mapping_diagnostics(self) -> pd.DataFrame:
        """Return per-site ambiguous mapping diagnostics."""

        return self.mapping_diagnostics.ambiguous_mapping

    @property
    def eligible_site_keys(self) -> tuple[str, ...]:
        return tuple(
            row.site_key
            for row in self.site_eligibility
            if row.eligibility
            is ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION
        )

    @property
    def fallback_site_keys(self) -> tuple[str, ...]:
        return tuple(
            row.site_key
            for row in self.site_eligibility
            if row.eligibility
            is ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY
        )

    @property
    def excluded_site_keys(self) -> tuple[str, ...]:
        return tuple(
            row.site_key
            for row in self.site_eligibility
            if row.eligibility
            is ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION
        )

    def site_eligibility_dataframe(self) -> pd.DataFrame:
        """Return a defensive site-eligibility table snapshot."""

        return export_dataframe(self._site_eligibility_table)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible audit report payload."""

        return {
            "schema_version": int(self.schema_version),
            "preparation_policy": self.preparation_policy,
            "protein_mapping_policy": self.protein_mapping_policy,
            "policy_parameters": _mapping_payload(self.policy_parameters),
            "sample_alignment": self.sample_alignment.to_payload(),
            "transformation_state": (
                None
                if self.transformation_state is None
                else self.transformation_state.to_payload()
            ),
            "site_eligibility": [row.to_payload() for row in self.site_eligibility],
            "site_eligibility_table": _dataframe_records_payload(
                self._site_eligibility_table
            ),
            "eligible_site_keys": list(self.eligible_site_keys),
            "fallback_site_keys": list(self.fallback_site_keys),
            "excluded_site_keys": list(self.excluded_site_keys),
            "mapping_diagnostics": self.mapping_diagnostics.to_payload(),
            "provenance": (
                None
                if self.provenance is None
                else provenance_to_payload(self.provenance)
            ),
        }

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another report owns the same report content."""

        if not isinstance(other, ProteinAwarePreparationReport):
            return False
        return (
            self.site_eligibility == other.site_eligibility
            and self.mapping_diagnostics.scientifically_equals(
                other.mapping_diagnostics
            )
            and self.sample_alignment == other.sample_alignment
            and self.transformation_state == other.transformation_state
            and self.preparation_policy == other.preparation_policy
            and self.protein_mapping_policy == other.protein_mapping_policy
            and self.policy_parameters == other.policy_parameters
            and self.provenance == other.provenance
            and self.schema_version == other.schema_version
            and dataframe_equals(
                self._site_eligibility_table,
                other._site_eligibility_table,
            )
        )


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ProteinAwarePreparationResult:
    """Aligned protein-aware preparation inputs and audit report."""

    __hash__ = object.__hash__

    report: ProteinAwarePreparationReport
    _matched_pairs: pd.DataFrame = field(init=False, repr=False)
    _protein_covariate_matrix: pd.DataFrame = field(init=False, repr=False)

    def __init__(
        self,
        *,
        matched_pairs: pd.DataFrame,
        protein_covariate_matrix: pd.DataFrame,
        report: ProteinAwarePreparationReport,
        _assume_owned: bool = False,
    ) -> None:
        _require_instance(
            report,
            expected_type=ProteinAwarePreparationReport,
            error_message=(
                "protein_aware_preparation_result.report must be "
                "ProteinAwarePreparationReport"
            ),
        )
        pairs = own_dataframe(
            matched_pairs,
            field_name="protein_aware_preparation_result.matched_pairs",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        covariates = own_dataframe(
            protein_covariate_matrix,
            field_name="protein_aware_preparation_result.protein_covariate_matrix",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        _validate_matched_pairs(pairs)
        _validate_protein_covariate_matrix(
            covariates,
            sample_columns=report.sample_alignment.phospho_sample_columns,
        )
        _validate_result_alignment(
            matched_pairs=pairs,
            protein_covariate_matrix=covariates,
            report=report,
        )
        object.__setattr__(self, "_matched_pairs", pairs)
        object.__setattr__(self, "_protein_covariate_matrix", covariates)
        object.__setattr__(self, "report", report)

    @property
    def matched_pairs(self) -> pd.DataFrame:
        """Return matched phosphosite/protein-row pairs."""

        return export_dataframe(self._matched_pairs)

    @property
    def protein_covariate_matrix(self) -> pd.DataFrame:
        """Return the sample-aligned total-protein covariate matrix."""

        return export_dataframe(self._protein_covariate_matrix)

    @property
    def site_eligibility_table(self) -> pd.DataFrame:
        return self.report.site_eligibility_table

    @property
    def missing_protein_abundance_diagnostics(self) -> pd.DataFrame:
        return self.report.missing_protein_abundance_diagnostics

    @property
    def ambiguous_mapping_diagnostics(self) -> pd.DataFrame:
        return self.report.ambiguous_mapping_diagnostics

    @property
    def sample_alignment_diagnostics(self) -> ProteinAwareSampleAlignmentDiagnostics:
        return self.report.sample_alignment_diagnostics

    @property
    def preparation_policy(self) -> str:
        return self.report.preparation_policy

    @property
    def protein_mapping_policy(self) -> str:
        return self.report.protein_mapping_policy

    @property
    def provenance(self) -> RunProvenance | None:
        return self.report.provenance

    def matched_pairs_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self._matched_pairs)

    def protein_covariate_matrix_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self._protein_covariate_matrix)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible preparation result payload."""

        return {
            "matched_pairs": _dataframe_records_payload(self._matched_pairs),
            "protein_covariate_matrix": _matrix_payload(self._protein_covariate_matrix),
            "protein_covariate_matrix_shape": [
                int(self._protein_covariate_matrix.shape[0]),
                int(self._protein_covariate_matrix.shape[1]),
            ],
            "report": self.report.to_payload(),
            "preparation_policy": self.preparation_policy,
            "protein_mapping_policy": self.protein_mapping_policy,
            "provenance": (
                None
                if self.provenance is None
                else provenance_to_payload(self.provenance)
            ),
        }

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another result owns the same preparation content."""

        if not isinstance(other, ProteinAwarePreparationResult):
            return False
        return (
            self.report.scientifically_equals(other.report)
            and dataframe_equals(self._matched_pairs, other._matched_pairs)
            and dataframe_equals(
                self._protein_covariate_matrix,
                other._protein_covariate_matrix,
            )
        )

    @classmethod
    def _from_owned(
        cls,
        *,
        matched_pairs: pd.DataFrame,
        protein_covariate_matrix: pd.DataFrame,
        report: ProteinAwarePreparationReport,
    ) -> ProteinAwarePreparationResult:
        return cls(
            matched_pairs=matched_pairs,
            protein_covariate_matrix=protein_covariate_matrix,
            report=report,
            _assume_owned=True,
        )


def _normalize_site_eligibility(
    values: Sequence[ProteinAwareSiteEligibility],
) -> tuple[ProteinAwareSiteEligibility, ...]:
    raw_values = cast(object, values)
    if not isinstance(raw_values, Sequence) or isinstance(raw_values, (str, bytes)):
        raise DatasetValidationError(
            "protein_aware_preparation_report.site_eligibility must be a sequence "
            "of ProteinAwareSiteEligibility rows"
        )
    raw_rows = tuple(cast(Sequence[object], raw_values))
    rows: list[ProteinAwareSiteEligibility] = []
    for position, raw_row in enumerate(raw_rows):
        if not isinstance(raw_row, ProteinAwareSiteEligibility):
            raise DatasetValidationError(
                "protein_aware_preparation_report.site_eligibility entries must be "
                f"ProteinAwareSiteEligibility; invalid_position={position}"
            )
        rows.append(raw_row)
    resolved_rows = tuple(rows)
    duplicate_keys = _duplicate_values(tuple(row.site_key for row in resolved_rows))
    if duplicate_keys:
        preview = ", ".join(repr(value) for value in duplicate_keys[:5])
        suffix = "" if len(duplicate_keys) <= 5 else " ..."
        raise DatasetValidationError(
            "protein_aware_preparation_report.site_eligibility site_key values "
            f"must be unique; duplicate_site_keys={preview}{suffix}"
        )
    return resolved_rows


def _duplicate_values(values: tuple[str, ...]) -> list[str]:
    seen: set[str] = set()
    duplicates: list[str] = []
    duplicate_seen: set[str] = set()
    for value in values:
        if value not in seen:
            seen.add(value)
            continue
        if value not in duplicate_seen:
            duplicates.append(value)
            duplicate_seen.add(value)
    return duplicates


def _site_eligibility_table(
    rows: tuple[ProteinAwareSiteEligibility, ...],
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "site_key": row.site_key,
                "eligibility": row.eligibility.value,
                "mapping_status": row.mapping_status.value,
                "protein_identifier": row.protein_identifier,
                "total_protein_row_key": row.total_protein_row_key,
                "reasons": tuple(row.reasons),
            }
            for row in rows
        ],
        columns=list(PROTEIN_AWARE_SITE_ELIGIBILITY_COLUMNS),
    )


def _validate_matched_pairs(matched_pairs: pd.DataFrame) -> None:
    require_columns(
        matched_pairs,
        field_name="protein_aware_preparation_result.matched_pairs",
        required_columns=PROTEIN_AWARE_MATCHED_PAIR_COLUMNS,
        error_type=DatasetValidationError,
    )
    for column in PROTEIN_AWARE_MATCHED_PAIR_COLUMNS:
        _require_non_empty_string_column(
            matched_pairs,
            field_name="protein_aware_preparation_result.matched_pairs",
            column_name=column,
        )


def _validate_protein_covariate_matrix(
    protein_covariate_matrix: pd.DataFrame,
    *,
    sample_columns: tuple[str, ...],
) -> None:
    field_name = "protein_aware_preparation_result.protein_covariate_matrix"
    require_unique_index(
        protein_covariate_matrix,
        field_name=field_name,
        error_type=DatasetValidationError,
    )
    require_unique_columns(
        protein_covariate_matrix,
        field_name=field_name,
        error_type=DatasetValidationError,
    )
    require_numeric_dataframe(
        protein_covariate_matrix,
        field_name=field_name,
        error_type=DatasetValidationError,
    )
    if sample_columns:
        require_exact_index_match(
            left=protein_covariate_matrix.columns,
            right=pd.Index(sample_columns),
            left_name=f"{field_name}.columns",
            right_name=(
                "protein_aware_preparation_report.sample_alignment."
                "phospho_sample_columns"
            ),
            error_type=DatasetValidationError,
        )


def _validate_result_alignment(
    *,
    matched_pairs: pd.DataFrame,
    protein_covariate_matrix: pd.DataFrame,
    report: ProteinAwarePreparationReport,
) -> None:
    matched_site_keys = tuple(matched_pairs.loc[:, "site_key"].astype(str).tolist())
    if matched_site_keys != report.eligible_site_keys:
        raise DatasetValidationError(
            "protein_aware_preparation_result.matched_pairs.site_key values must "
            "match report eligible_site_keys in order"
        )
    missing_rows = [
        row_key
        for row_key in matched_pairs.loc[:, "total_protein_row_key"]
        .astype(str)
        .tolist()
        if row_key not in protein_covariate_matrix.index
    ]
    if missing_rows:
        preview = ", ".join(repr(value) for value in dict.fromkeys(missing_rows))
        raise DatasetValidationError(
            "protein_aware_preparation_result.protein_covariate_matrix.index must "
            f"contain matched total_protein_row_key values; missing={preview}"
        )


def _require_non_empty_string_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
) -> None:
    invalid_count = sum(
        1
        for value in frame.loc[:, column_name].tolist()
        if not isinstance(value, str) or value.strip() == ""
    )
    if invalid_count:
        raise DatasetValidationError(
            f"{field_name}.{column_name} values must be non-empty strings; "
            f"invalid_count={invalid_count}"
        )


def _normalize_policy(
    value: object,
    *,
    field_name: str,
    supported: frozenset[str],
) -> str:
    normalized = _require_non_empty_string(value, field_name=field_name)
    if normalized not in supported:
        supported_values = ", ".join(sorted(repr(item) for item in supported))
        raise DatasetValidationError(f"{field_name} must be one of: {supported_values}")
    return normalized


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise DatasetValidationError(f"{field_name} must be a non-empty string")
    return value.strip()


def _require_instance(
    value: object,
    *,
    expected_type: type[object],
    error_message: str,
) -> None:
    if not isinstance(value, expected_type):
        raise DatasetValidationError(error_message)


def _dataframe_records_payload(frame: pd.DataFrame) -> list[dict[str, object]]:
    records = cast(list[dict[object, object]], frame.to_dict(orient="records"))
    return [
        {str(key): _payload_value(value) for key, value in record.items()}
        for record in records
    ]


def _matrix_payload(frame: pd.DataFrame) -> dict[str, object]:
    index_values = cast(list[object], frame.index.tolist())
    column_values = cast(list[object], frame.columns.tolist())
    matrix_values = cast(list[list[object]], frame.to_numpy(dtype="object").tolist())
    return {
        "index": [_payload_value(value) for value in index_values],
        "columns": [_payload_value(value) for value in column_values],
        "data": [[_payload_value(value) for value in row] for row in matrix_values],
    }


def _mapping_payload(value: Mapping[str, object]) -> dict[str, object]:
    return {str(key): _payload_value(item) for key, item in value.items()}


def _payload_value(value: object) -> object:
    if isinstance(value, tuple):
        return [_payload_value(item) for item in cast(tuple[object, ...], value)]
    if isinstance(value, list):
        return [_payload_value(item) for item in cast(list[object], value)]
    if isinstance(value, Mapping):
        return _mapping_payload(cast(Mapping[str, object], value))
    if _is_missing_value(value):
        return None
    return value


def _is_missing_value(value: object) -> bool:
    if isinstance(value, (list, tuple, dict, set, frozenset)):
        return False
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


__all__ = [
    "PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
]
