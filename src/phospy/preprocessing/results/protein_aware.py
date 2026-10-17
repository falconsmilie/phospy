"""Passive results produced by protein-aware preprocessing."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

import phospy.preparation_records.protein_aware as _protein_records
from phospy.errors.validation import DatasetValidationError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.frames.validation import (
    require_columns,
    require_exact_index_match,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.provenance.immutability import freeze_json_mapping_with_error_type


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ProteinAwarePreparationStageOutput:
    """Owned, provenance-free output from protein-aware preprocessing."""

    __hash__ = object.__hash__
    site_eligibility: tuple[_protein_records.ProteinAwareSiteEligibility, ...]
    mapping_diagnostics: _protein_records.ProteinAwareMappingDiagnostics
    sample_alignment: _protein_records.ProteinAwareSampleAlignmentDiagnostics
    transformation_state: (
        _protein_records.ProteinAwareTransformationStateDiagnostics | None
    )
    preparation_policy: str
    protein_mapping_policy: str
    policy_parameters: Mapping[str, object]
    _matched_pairs: pd.DataFrame = field(init=False, repr=False)
    _protein_covariate_matrix: pd.DataFrame = field(init=False, repr=False)

    def __init__(
        self,
        *,
        matched_pairs: pd.DataFrame,
        protein_covariate_matrix: pd.DataFrame,
        site_eligibility: Sequence[_protein_records.ProteinAwareSiteEligibility],
        mapping_diagnostics: _protein_records.ProteinAwareMappingDiagnostics,
        sample_alignment: _protein_records.ProteinAwareSampleAlignmentDiagnostics,
        transformation_state: _protein_records.ProteinAwareTransformationStateDiagnostics
        | None,
        preparation_policy: str,
        protein_mapping_policy: str,
        policy_parameters: Mapping[str, object],
        _assume_owned: bool = False,
    ) -> None:
        rows = _normalize_site_eligibility(site_eligibility)
        _require_instance(
            mapping_diagnostics,
            expected_type=_protein_records.ProteinAwareMappingDiagnostics,
            error_message="protein_aware_preparation_stage_output.mapping_diagnostics must be ProteinAwareMappingDiagnostics",
        )
        _require_instance(
            sample_alignment,
            expected_type=_protein_records.ProteinAwareSampleAlignmentDiagnostics,
            error_message="protein_aware_preparation_stage_output.sample_alignment must be ProteinAwareSampleAlignmentDiagnostics",
        )
        if transformation_state is not None:
            _require_instance(
                transformation_state,
                expected_type=_protein_records.ProteinAwareTransformationStateDiagnostics,
                error_message="protein_aware_preparation_stage_output.transformation_state must be ProteinAwareTransformationStateDiagnostics or None",
            )
        resolved_policy = _normalize_policy(
            preparation_policy,
            field_name="protein_aware_preparation_stage_output.preparation_policy",
            supported=_protein_records.DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
        )
        resolved_mapping_policy = _normalize_policy(
            protein_mapping_policy,
            field_name="protein_aware_preparation_stage_output.protein_mapping_policy",
            supported=_protein_records.DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
        )
        parameters = freeze_json_mapping_with_error_type(
            policy_parameters,
            field_name="protein_aware_preparation_stage_output.policy_parameters",
            error_type=DatasetValidationError,
        )
        pairs = own_dataframe(
            matched_pairs,
            field_name="protein_aware_preparation_stage_output.matched_pairs",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        covariates = own_dataframe(
            protein_covariate_matrix,
            field_name="protein_aware_preparation_stage_output.protein_covariate_matrix",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        _validate_matched_pairs(pairs)
        _validate_protein_covariate_matrix(
            covariates, sample_columns=sample_alignment.phospho_sample_columns
        )
        eligible = tuple(
            row.site_key
            for row in rows
            if row.eligibility
            is _protein_records.ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION
        )
        _validate_stage_alignment(
            matched_pairs=pairs,
            protein_covariate_matrix=covariates,
            eligible_site_keys=eligible,
        )
        object.__setattr__(self, "_matched_pairs", pairs)
        object.__setattr__(self, "_protein_covariate_matrix", covariates)
        object.__setattr__(self, "site_eligibility", rows)
        object.__setattr__(self, "mapping_diagnostics", mapping_diagnostics)
        object.__setattr__(self, "sample_alignment", sample_alignment)
        object.__setattr__(self, "transformation_state", transformation_state)
        object.__setattr__(self, "preparation_policy", resolved_policy)
        object.__setattr__(self, "protein_mapping_policy", resolved_mapping_policy)
        object.__setattr__(self, "policy_parameters", parameters)

    @property
    def matched_pairs(self) -> pd.DataFrame:
        return export_dataframe(self._matched_pairs)

    @property
    def protein_covariate_matrix(self) -> pd.DataFrame:
        return export_dataframe(self._protein_covariate_matrix)

    @classmethod
    def _from_owned(
        cls,
        *,
        matched_pairs: pd.DataFrame,
        protein_covariate_matrix: pd.DataFrame,
        site_eligibility: Sequence[_protein_records.ProteinAwareSiteEligibility],
        mapping_diagnostics: _protein_records.ProteinAwareMappingDiagnostics,
        sample_alignment: _protein_records.ProteinAwareSampleAlignmentDiagnostics,
        transformation_state: _protein_records.ProteinAwareTransformationStateDiagnostics
        | None,
        preparation_policy: str,
        protein_mapping_policy: str,
        policy_parameters: Mapping[str, object],
    ) -> ProteinAwarePreparationStageOutput:
        return cls(
            matched_pairs=matched_pairs,
            protein_covariate_matrix=protein_covariate_matrix,
            site_eligibility=site_eligibility,
            mapping_diagnostics=mapping_diagnostics,
            sample_alignment=sample_alignment,
            transformation_state=transformation_state,
            preparation_policy=preparation_policy,
            protein_mapping_policy=protein_mapping_policy,
            policy_parameters=policy_parameters,
            _assume_owned=True,
        )


def _normalize_site_eligibility(
    values: Sequence[_protein_records.ProteinAwareSiteEligibility],
) -> tuple[_protein_records.ProteinAwareSiteEligibility, ...]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise DatasetValidationError(
            "protein_aware_preparation_stage_output.site_eligibility must be a sequence of ProteinAwareSiteEligibility rows"
        )
    rows = tuple(values)
    for position, row in enumerate(rows):
        if not isinstance(row, _protein_records.ProteinAwareSiteEligibility):
            raise DatasetValidationError(
                "protein_aware_preparation_stage_output.site_eligibility entries must be ProteinAwareSiteEligibility; "
                f"invalid_position={position}"
            )
    duplicate_keys = _duplicate_values(tuple(row.site_key for row in rows))
    if duplicate_keys:
        preview = ", ".join(repr(value) for value in duplicate_keys[:5])
        suffix = "" if len(duplicate_keys) <= 5 else " ..."
        raise DatasetValidationError(
            "protein_aware_preparation_stage_output.site_eligibility site_key values must be unique; "
            f"duplicate_site_keys={preview}{suffix}"
        )
    return rows


def _duplicate_values(values: tuple[str, ...]) -> list[str]:
    seen: set[str] = set()
    duplicate_seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
        elif value not in duplicate_seen:
            duplicate_seen.add(value)
            duplicates.append(value)
    return duplicates


def _validate_matched_pairs(frame: pd.DataFrame) -> None:
    field_name = "protein_aware_preparation_stage_output.matched_pairs"
    require_columns(
        frame,
        field_name=field_name,
        required_columns=_protein_records.PROTEIN_AWARE_MATCHED_PAIR_COLUMNS,
        error_type=DatasetValidationError,
    )
    for column in _protein_records.PROTEIN_AWARE_MATCHED_PAIR_COLUMNS:
        invalid_count = sum(
            1
            for value in frame.loc[:, column].tolist()
            if not isinstance(value, str) or value.strip() == ""
        )
        if invalid_count:
            raise DatasetValidationError(
                f"{field_name}.{column} values must be non-empty strings; invalid_count={invalid_count}"
            )


def _validate_protein_covariate_matrix(
    frame: pd.DataFrame, *, sample_columns: tuple[str, ...]
) -> None:
    field_name = "protein_aware_preparation_stage_output.protein_covariate_matrix"
    require_unique_index(
        frame, field_name=field_name, error_type=DatasetValidationError
    )
    require_unique_columns(
        frame, field_name=field_name, error_type=DatasetValidationError
    )
    require_numeric_dataframe(
        frame, field_name=field_name, error_type=DatasetValidationError
    )
    if sample_columns:
        require_exact_index_match(
            left=frame.columns,
            right=pd.Index(sample_columns),
            left_name=f"{field_name}.columns",
            right_name="protein_aware_preparation_stage_output.sample_alignment.phospho_sample_columns",
            error_type=DatasetValidationError,
        )


def _validate_stage_alignment(
    *,
    matched_pairs: pd.DataFrame,
    protein_covariate_matrix: pd.DataFrame,
    eligible_site_keys: tuple[str, ...],
) -> None:
    matched_site_keys = tuple(matched_pairs.loc[:, "site_key"].astype(str).tolist())
    if matched_site_keys != eligible_site_keys:
        raise DatasetValidationError(
            "protein_aware_preparation_stage_output.matched_pairs.site_key values must match eligible site keys in order"
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
            "protein_aware_preparation_stage_output.protein_covariate_matrix.index must contain matched total_protein_row_key values; "
            f"missing={preview}"
        )


def _normalize_policy(
    value: object, *, field_name: str, supported: frozenset[str]
) -> str:
    normalized = _protein_records._require_non_empty_string(
        value,
        field_name=field_name,
    )
    if normalized not in supported:
        supported_values = ", ".join(sorted(repr(item) for item in supported))
        raise DatasetValidationError(f"{field_name} must be one of: {supported_values}")
    return normalized


def _require_instance(
    value: object, *, expected_type: type[object], error_message: str
) -> None:
    if not isinstance(value, expected_type):
        raise DatasetValidationError(error_message)


__all__ = [
    "ProteinAwarePreparationStageOutput",
]
