"""Passive protein-aware preparation records and persisted vocabulary."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.comparison import dataframe_equals
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.frames.validation import require_columns
from phospy.policies import PolicyEnum

DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED = "disabled"
DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS = "prepare_model_inputs"
DatasetProteinAwarePreparationPolicy = Literal["disabled", "prepare_model_inputs"]
DATASET_PROTEIN_AWARE_PREPARATION_POLICIES = frozenset(
    {
        DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED,
        DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS,
    }
)

DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS = (
    "require_unambiguous"
)
DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT = (
    "allow_missing_with_report"
)
DatasetProteinAwarePreparationMappingPolicy = Literal[
    "require_unambiguous",
    "allow_missing_with_report",
]
DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES = frozenset(
    {
        DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
        DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT,
    }
)

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


class ProteinMappingStatus(PolicyEnum):
    """Per-site status for protein-aware mapping resolution."""

    MATCHED = "matched"
    MISSING_SITE_PROTEIN_IDENTIFIER = "missing_site_protein_identifier"
    MISSING_TOTAL_PROTEIN_ROW = "missing_total_protein_row"
    AMBIGUOUS_SITE_PROTEIN_MAPPING = "ambiguous_site_protein_mapping"
    AMBIGUOUS_TOTAL_PROTEIN_MAPPING = "ambiguous_total_protein_mapping"


class ProteinAwarePreparationEligibility(PolicyEnum):
    """Per-site protein-aware preparation eligibility status."""

    ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION = "eligible_for_protein_aware_preparation"
    FALLBACK_TO_PHOSPHO_ONLY = "fallback_to_phospho_only"
    EXCLUDED_FROM_PREPARATION = "excluded_from_preparation"


@dataclass(frozen=True, slots=True)
class ProteinAwareSampleAlignmentDiagnostics:
    """Diagnostics describing phospho/total sample-column compatibility."""

    phospho_sample_columns: tuple[str, ...]
    total_protein_sample_columns: tuple[str, ...]
    exact_sample_order_match: bool
    sample_order_compatible: bool
    reordered_sample_columns: bool
    allow_reordered_samples: bool
    missing_total_protein_samples: tuple[str, ...]
    extra_total_protein_samples: tuple[str, ...]

    def to_payload(self) -> dict[str, object]:
        return {
            "phospho_sample_columns": list(self.phospho_sample_columns),
            "total_protein_sample_columns": list(self.total_protein_sample_columns),
            "exact_sample_order_match": self.exact_sample_order_match,
            "sample_order_compatible": self.sample_order_compatible,
            "reordered_sample_columns": self.reordered_sample_columns,
            "allow_reordered_samples": self.allow_reordered_samples,
            "missing_total_protein_samples": list(self.missing_total_protein_samples),
            "extra_total_protein_samples": list(self.extra_total_protein_samples),
        }


@dataclass(frozen=True, slots=True)
class ProteinAwareTransformationStateDiagnostics:
    """Diagnostics describing phospho/total transformation-state compatibility."""

    compatible: bool
    phospho_transformation_state: dict[str, object] | None
    total_protein_transformation_state: dict[str, object] | None

    def to_payload(self) -> dict[str, object]:
        return {
            "transformation_state_compatible": self.compatible,
            "phospho_transformation_state": self.phospho_transformation_state,
            "total_protein_transformation_state": (
                self.total_protein_transformation_state
            ),
        }


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
                self.site_key, field_name="protein_aware_site_eligibility.site_key"
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
                self.reasons, field_name="protein_aware_site_eligibility.reasons"
            ),
        )

    def to_payload(self) -> dict[str, object]:
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
            field_name="protein_aware_mapping_diagnostics.missing_protein_abundance",
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
        return export_dataframe(self._missing_protein_abundance)

    @property
    def ambiguous_mapping(self) -> pd.DataFrame:
        return export_dataframe(self._ambiguous_mapping)

    def missing_protein_abundance_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self._missing_protein_abundance)

    def ambiguous_mapping_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self._ambiguous_mapping)

    def to_payload(self) -> dict[str, object]:
        return {
            "missing_protein_abundance": _dataframe_records_payload(
                self._missing_protein_abundance
            ),
            "ambiguous_mapping": _dataframe_records_payload(self._ambiguous_mapping),
        }

    def scientifically_equals(self, other: object) -> bool:
        if not isinstance(other, ProteinAwareMappingDiagnostics):
            return False
        return dataframe_equals(
            self._missing_protein_abundance, other._missing_protein_abundance
        ) and dataframe_equals(self._ambiguous_mapping, other._ambiguous_mapping)

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


def _normalize_string_tuple(
    value: Sequence[str], *, field_name: str
) -> tuple[str, ...]:
    if isinstance(value, str):
        value = (value,)
    if not isinstance(value, Sequence) or isinstance(value, bytes):
        raise DatasetValidationError(f"{field_name} must be a sequence of strings")
    return tuple(
        _require_non_empty_string(item, field_name=field_name) for item in value
    )


def _optional_non_empty_string(value: object | None, *, field_name: str) -> str | None:
    return (
        None
        if value is None
        else _require_non_empty_string(value, field_name=field_name)
    )


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise DatasetValidationError(f"{field_name} must be a non-empty string")
    return value.strip()


def _dataframe_records_payload(frame: pd.DataFrame) -> list[dict[str, object]]:
    return [
        {str(key): _payload_value(value) for key, value in record.items()}
        for record in frame.to_dict(orient="records")
    ]


def _payload_value(value: object) -> object:
    if isinstance(value, tuple):
        return [_payload_value(item) for item in value]
    if isinstance(value, list):
        return [_payload_value(item) for item in value]
    if isinstance(value, Mapping):
        return {str(key): _payload_value(item) for key, item in value.items()}
    if not isinstance(value, (list, tuple, dict, set, frozenset)) and bool(
        pd.Series((value,), dtype="object").isna().iat[0]
    ):
        return None
    return value


__all__ = [
    "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES",
    "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT",
    "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS",
    "DATASET_PROTEIN_AWARE_PREPARATION_POLICIES",
    "DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED",
    "DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS",
    "DatasetProteinAwarePreparationMappingPolicy",
    "DatasetProteinAwarePreparationPolicy",
    "PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS",
    "PROTEIN_AWARE_MATCHED_PAIR_COLUMNS",
    "PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS",
    "PROTEIN_AWARE_SITE_ELIGIBILITY_COLUMNS",
    "ProteinAwareMappingDiagnostics",
    "ProteinAwarePreparationEligibility",
    "ProteinAwareSampleAlignmentDiagnostics",
    "ProteinAwareSiteEligibility",
    "ProteinAwareTransformationStateDiagnostics",
    "ProteinMappingStatus",
]
