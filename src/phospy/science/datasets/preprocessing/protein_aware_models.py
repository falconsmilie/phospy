"""Typed protein-aware preparation result and audit report models."""

from __future__ import annotations

__phospy_contracts_facade_role__ = "science_owned_public_model"

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

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
from phospy.provenance.immutability import freeze_json_mapping_with_error_type
from phospy.provenance.models import RunProvenance
from phospy.provenance.serialization import to_payload as provenance_to_payload
from phospy.science.configs.preprocessing.total_protein import (
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS,
)
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    ProteinAwarePreparationEligibility,
    ProteinAwareSampleAlignmentDiagnostics,
    ProteinAwareTransformationStateDiagnostics,
)
from phospy.science.datasets.preprocessing.protein_mapping import (
    ProteinMappingStatus,
)

PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION = 1

PROTEIN_AWARE_MATCHED_PAIR_COLUMNS = (
    "site_key",
    "protein_identifier",
    "total_protein_row_key",
)

PROTEIN_AWARE_SITE_ELIGIBILITY_COLUMNS = (
    "site_key",
    "eligibility",
    "mapping_status",
    "protein_identifier",
    "total_protein_row_key",
    "reasons",
)

PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS = (
    "site_key",
    "protein_identifier",
    "mapping_status",
    "reason",
)

PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS = (
    "site_key",
    "mapping_status",
    "protein_identifier",
    "candidate_protein_identifiers",
    "candidate_total_protein_row_keys",
    "reason",
)


@dataclass(frozen=True, slots=True)
class ProteinAwareSiteEligibility:
    """Typed per-site protein-aware preparation eligibility row."""

    site_key: str
    eligibility: ProteinAwarePreparationEligibility
    mapping_status: ProteinMappingStatus
    protein_identifier: str | None = None
    total_protein_row_key: str | None = None
    reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "site_key",
            _require_non_empty_string(
                self.site_key,
                field_name="protein_aware_site_eligibility.site_key",
            ),
        )
        object.__setattr__(
            self,
            "eligibility",
            ProteinAwarePreparationEligibility.parse(
                self.eligibility,
                field_name="protein_aware_site_eligibility.eligibility",
            ),
        )
        object.__setattr__(
            self,
            "mapping_status",
            ProteinMappingStatus.parse(
                self.mapping_status,
                field_name="protein_aware_site_eligibility.mapping_status",
            ),
        )
        object.__setattr__(
            self,
            "protein_identifier",
            _optional_non_empty_string(
                self.protein_identifier,
                field_name="protein_aware_site_eligibility.protein_identifier",
            ),
        )
        object.__setattr__(
            self,
            "total_protein_row_key",
            _optional_non_empty_string(
                self.total_protein_row_key,
                field_name="protein_aware_site_eligibility.total_protein_row_key",
            ),
        )
        object.__setattr__(
            self,
            "reasons",
            _normalize_string_tuple(
                self.reasons,
                field_name="protein_aware_site_eligibility.reasons",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible eligibility payload."""

        return {
            "site_key": self.site_key,
            "eligibility": self.eligibility.value,
            "mapping_status": self.mapping_status.value,
            "protein_identifier": self.protein_identifier,
            "total_protein_row_key": self.total_protein_row_key,
            "reasons": list(self.reasons),
        }


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ProteinAwareMappingDiagnostics:
    """Machine-readable mapping diagnostics for preparation audit."""

    __hash__ = object.__hash__

    _missing_protein_abundance: pd.DataFrame = field(init=False, repr=False)
    _ambiguous_mapping: pd.DataFrame = field(init=False, repr=False)

    def __init__(
        self,
        *,
        missing_protein_abundance: pd.DataFrame | None = None,
        ambiguous_mapping: pd.DataFrame | None = None,
        _assume_owned: bool = False,
    ) -> None:
        missing = _owned_diagnostic_table(
            missing_protein_abundance,
            columns=PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS,
            field_name=("protein_aware_mapping_diagnostics.missing_protein_abundance"),
            assume_owned=_assume_owned,
        )
        ambiguous = _owned_diagnostic_table(
            ambiguous_mapping,
            columns=PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS,
            field_name="protein_aware_mapping_diagnostics.ambiguous_mapping",
            assume_owned=_assume_owned,
        )
        object.__setattr__(self, "_missing_protein_abundance", missing)
        object.__setattr__(self, "_ambiguous_mapping", ambiguous)

    @property
    def missing_protein_abundance(self) -> pd.DataFrame:
        """Return per-site missing protein-abundance diagnostics."""

        return export_dataframe(self._missing_protein_abundance)

    @property
    def ambiguous_mapping(self) -> pd.DataFrame:
        """Return per-site ambiguous phosphosite/protein mapping diagnostics."""

        return export_dataframe(self._ambiguous_mapping)

    def missing_protein_abundance_dataframe(self) -> pd.DataFrame:
        """Return a defensive missing-protein diagnostics snapshot."""

        return export_dataframe(self._missing_protein_abundance)

    def ambiguous_mapping_dataframe(self) -> pd.DataFrame:
        """Return a defensive ambiguous-mapping diagnostics snapshot."""

        return export_dataframe(self._ambiguous_mapping)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible diagnostics payload."""

        return {
            "missing_protein_abundance": _dataframe_records_payload(
                self._missing_protein_abundance
            ),
            "ambiguous_mapping": _dataframe_records_payload(self._ambiguous_mapping),
        }

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another diagnostics object owns the same tables."""

        if not isinstance(other, ProteinAwareMappingDiagnostics):
            return False
        return dataframe_equals(
            self._missing_protein_abundance,
            other._missing_protein_abundance,
        ) and dataframe_equals(
            self._ambiguous_mapping,
            other._ambiguous_mapping,
        )

    @classmethod
    def _from_owned(
        cls,
        *,
        missing_protein_abundance: pd.DataFrame | None = None,
        ambiguous_mapping: pd.DataFrame | None = None,
    ) -> ProteinAwareMappingDiagnostics:
        return cls(
            missing_protein_abundance=missing_protein_abundance,
            ambiguous_mapping=ambiguous_mapping,
            _assume_owned=True,
        )


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ProteinAwarePreparationReport:
    """Protein-aware preparation audit report.

    This report records prepared input eligibility, mapping diagnostics, sample
    alignment, and policy provenance only. It does not resolve mappings, subtract
    total protein, fit a model, or run differential analysis.
    """

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
        if not isinstance(schema_version, int) or isinstance(schema_version, bool):
            raise DatasetValidationError(
                "protein_aware_preparation_report.schema_version must be an int"
            )
        object.__setattr__(self, "site_eligibility", resolved_site_eligibility)
        object.__setattr__(
            self,
            "mapping_diagnostics",
            resolved_mapping_diagnostics,
        )
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

        payload: dict[str, object] = {
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
        return payload

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
    """Aligned protein-aware preparation inputs and audit report.

    The protein covariate matrix is an input contract for future differential
    modelling. Constructing this result does not run modelling and does not
    change `AnalysisReadyPhosphoDataset`.
    """

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
        """Return the per-site eligibility table."""

        return self.report.site_eligibility_table

    @property
    def missing_protein_abundance_diagnostics(self) -> pd.DataFrame:
        """Return missing protein-abundance diagnostics."""

        return self.report.missing_protein_abundance_diagnostics

    @property
    def ambiguous_mapping_diagnostics(self) -> pd.DataFrame:
        """Return ambiguous mapping diagnostics."""

        return self.report.ambiguous_mapping_diagnostics

    @property
    def sample_alignment_diagnostics(self) -> ProteinAwareSampleAlignmentDiagnostics:
        """Return sample-alignment diagnostics."""

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
        """Return a defensive matched-pairs table snapshot."""

        return export_dataframe(self._matched_pairs)

    def protein_covariate_matrix_dataframe(self) -> pd.DataFrame:
        """Return a defensive protein-covariate matrix snapshot."""

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


def _owned_diagnostic_table(
    value: pd.DataFrame | None,
    *,
    columns: tuple[str, ...],
    field_name: str,
    assume_owned: bool,
) -> pd.DataFrame:
    if value is None:
        return pd.DataFrame(columns=list(columns))
    table = own_dataframe(
        value,
        field_name=field_name,
        error_type=DatasetValidationError,
        assume_owned=assume_owned,
    )
    require_columns(
        table,
        field_name=field_name,
        required_columns=columns,
        error_type=DatasetValidationError,
    )
    return table


def _normalize_site_eligibility(
    values: Sequence[ProteinAwareSiteEligibility],
) -> tuple[ProteinAwareSiteEligibility, ...]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise DatasetValidationError(
            "protein_aware_preparation_report.site_eligibility must be a sequence "
            "of ProteinAwareSiteEligibility rows"
        )
    rows = tuple(values)
    for position, row in enumerate(rows):
        if not isinstance(row, ProteinAwareSiteEligibility):
            raise DatasetValidationError(
                "protein_aware_preparation_report.site_eligibility entries must be "
                f"ProteinAwareSiteEligibility; invalid_position={position}"
            )
    site_keys = [row.site_key for row in rows]
    duplicate_keys = [
        site_key
        for site_key in dict.fromkeys(site_keys)
        if site_keys.count(site_key) > 1
    ]
    if duplicate_keys:
        preview = ", ".join(repr(value) for value in duplicate_keys[:5])
        suffix = "" if len(duplicate_keys) <= 5 else " ..."
        raise DatasetValidationError(
            "protein_aware_preparation_report.site_eligibility site_key values "
            f"must be unique; duplicate_site_keys={preview}{suffix}"
        )
    return rows


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
    require_unique_index(
        protein_covariate_matrix,
        field_name="protein_aware_preparation_result.protein_covariate_matrix",
        error_type=DatasetValidationError,
    )
    require_unique_columns(
        protein_covariate_matrix,
        field_name="protein_aware_preparation_result.protein_covariate_matrix",
        error_type=DatasetValidationError,
    )
    require_numeric_dataframe(
        protein_covariate_matrix,
        field_name="protein_aware_preparation_result.protein_covariate_matrix",
        error_type=DatasetValidationError,
    )
    if sample_columns:
        require_exact_index_match(
            left=protein_covariate_matrix.columns,
            right=pd.Index(sample_columns),
            left_name=(
                "protein_aware_preparation_result.protein_covariate_matrix.columns"
            ),
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
    invalid_count = 0
    for value in frame.loc[:, column_name].tolist():
        if not isinstance(value, str) or value.strip() == "":
            invalid_count += 1
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


def _normalize_string_tuple(
    value: Sequence[str],
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(value, str):
        value = (value,)
    if not isinstance(value, Sequence) or isinstance(value, bytes):
        raise DatasetValidationError(f"{field_name} must be a sequence of strings")
    normalized: list[str] = []
    for item in value:
        normalized.append(_require_non_empty_string(item, field_name=field_name))
    return tuple(normalized)


def _optional_non_empty_string(value: object | None, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_non_empty_string(value, field_name=field_name)


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
    return [
        {str(key): _payload_value(value) for key, value in record.items()}
        for record in frame.to_dict(orient="records")
    ]


def _matrix_payload(frame: pd.DataFrame) -> dict[str, object]:
    return {
        "index": [_payload_value(value) for value in frame.index.tolist()],
        "columns": [_payload_value(value) for value in frame.columns.tolist()],
        "data": [
            [_payload_value(value) for value in row]
            for row in frame.to_numpy(dtype="object").tolist()
        ],
    }


def _mapping_payload(value: Mapping[str, object]) -> dict[str, object]:
    return {str(key): _payload_value(item) for key, item in value.items()}


def _payload_value(value: object) -> object:
    if isinstance(value, tuple):
        return [_payload_value(item) for item in value]
    if isinstance(value, list):
        return [_payload_value(item) for item in value]
    if isinstance(value, Mapping):
        return _mapping_payload(value)
    if _is_missing_value(value):
        return None
    return value


def _is_missing_value(value: object) -> bool:
    if isinstance(value, (list, tuple, dict, set, frozenset)):
        return False
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


__all__ = [
    "PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS",
    "PROTEIN_AWARE_MATCHED_PAIR_COLUMNS",
    "PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS",
    "PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION",
    "PROTEIN_AWARE_SITE_ELIGIBILITY_COLUMNS",
    "ProteinAwareMappingDiagnostics",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
    "ProteinAwareSiteEligibility",
]
