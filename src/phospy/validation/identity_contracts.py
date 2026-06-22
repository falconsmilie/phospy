"""Shared phosphosite identity contracts for dataset, workflow, and result tables."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TypeVar, cast

import pandas as pd

from phospy.science.sites.identifiers import (
    ParsedSiteToken,
    canonicalize_site_components,
    canonicalize_site_identifier,
    canonicalize_site_series,
    try_parse_site_token,
)
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    decode_site_key,
    encode_site_key,
    require_positive_integer_position,
)
from phospy.science.sites.validation import (
    require_no_mixed_site_key_isoform_scope,
    require_site_key_index,
)
from phospy.validation.common.dataframes import (
    ValidationErrorType,
    require_columns,
    require_exact_index_match,
    require_non_empty_string_column,
)

ErrorType = TypeVar("ErrorType", bound=Exception)
_SITE_POSITION_CANDIDATE_COLUMNS = ("position", "site_position")
_SITE_KEY_INDEX_NAME = "site_key"
SITE_KEY_COLUMN = "site_key"
DISPLAY_ID_COLUMN = "display_id"
SITE_SEQUENCE_COLUMN = "site_sequence"
BASE_IDENTITY_COLUMNS = (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN)
PROTEIN_CONTEXT_COLUMNS = (
    "organism",
    "protein_namespace",
    "protein_identifier",
    "site",
)
DISPLAY_CONTEXT_COLUMNS = ("gene_symbol", "site")
ANALYSIS_READY_IDENTITY_COLUMNS = (
    SITE_KEY_COLUMN,
    DISPLAY_ID_COLUMN,
    "organism",
    "protein_namespace",
    "protein_identifier",
    "gene_symbol",
    "site",
    SITE_SEQUENCE_COLUMN,
)
RESULT_IDENTITY_COLUMNS = (
    SITE_KEY_COLUMN,
    DISPLAY_ID_COLUMN,
    "organism",
    "protein_namespace",
    "protein_identifier",
    "gene_symbol",
    "site",
)


class SequenceContextRequirement(str, Enum):
    """Sequence-context strictness for one identity boundary."""

    NONE = "none"
    PRESENT = "present"
    CENTRED = "centred"


@dataclass(frozen=True, slots=True)
class PhosphositeIdentityContract:
    """Reusable phosphosite identity requirements for one public boundary."""

    contract_id: str
    required_columns: tuple[str, ...] = BASE_IDENTITY_COLUMNS
    require_site_key_index: bool = True
    prefer_analysis_ready_index_diagnostics: bool = True
    require_site_key_column_index_coherence: bool = True
    check_site_key_column_index_before_uniqueness: bool = False
    require_unique_site_key: bool = True
    require_display_id: bool = True
    require_protein_context: bool = False
    require_site_key_metadata_coherence: bool = False
    sequence_context: SequenceContextRequirement = SequenceContextRequirement.NONE
    require_non_empty_sequence_context: bool = False
    allow_duplicate_display_id: bool = True


ANALYSIS_READY_DATASET_BASE_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="analysis_ready_dataset_base_identity",
    required_columns=ANALYSIS_READY_IDENTITY_COLUMNS,
    require_protein_context=True,
    sequence_context=SequenceContextRequirement.PRESENT,
)
ANALYSIS_READY_DATASET_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="analysis_ready_dataset_identity",
    required_columns=ANALYSIS_READY_IDENTITY_COLUMNS,
    require_protein_context=True,
    require_site_key_metadata_coherence=True,
    sequence_context=SequenceContextRequirement.PRESENT,
    require_non_empty_sequence_context=True,
)
WORKFLOW_INPUT_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="workflow_input_identity",
    required_columns=BASE_IDENTITY_COLUMNS,
    check_site_key_column_index_before_uniqueness=True,
)
WORKFLOW_PROTEIN_CONTEXT_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="workflow_protein_context_identity",
    required_columns=BASE_IDENTITY_COLUMNS + PROTEIN_CONTEXT_COLUMNS,
    check_site_key_column_index_before_uniqueness=True,
    require_protein_context=True,
    require_site_key_metadata_coherence=True,
)
WORKFLOW_SEQUENCE_CONTEXT_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="workflow_sequence_context_identity",
    required_columns=BASE_IDENTITY_COLUMNS,
    check_site_key_column_index_before_uniqueness=True,
    sequence_context=SequenceContextRequirement.CENTRED,
)
WORKFLOW_PROTEIN_SEQUENCE_CONTEXT_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="workflow_protein_sequence_context_identity",
    required_columns=BASE_IDENTITY_COLUMNS + PROTEIN_CONTEXT_COLUMNS,
    check_site_key_column_index_before_uniqueness=True,
    require_protein_context=True,
    require_site_key_metadata_coherence=True,
    sequence_context=SequenceContextRequirement.CENTRED,
)
RESULT_TABLE_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="result_table_identity",
    required_columns=RESULT_IDENTITY_COLUMNS,
    require_protein_context=True,
    prefer_analysis_ready_index_diagnostics=False,
)


def enforce_phosphosite_identity_contract(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    contract: PhosphositeIdentityContract,
    error_type: type[ErrorType],
    expected_index: pd.Index | None = None,
    expected_index_field_name: str | None = None,
    workflow_name: str | None = None,
    allow_opaque_site_values: bool = False,
    allow_gapped_sequence_context: bool = False,
    compare_raw_site_key_column_before_decode: bool = False,
) -> None:
    """Enforce a reusable identity contract without workflow-specific policy."""

    dataframe_error_type = cast(ValidationErrorType, error_type)
    require_columns(
        site_metadata,
        field_name=field_name,
        required_columns=contract.required_columns,
        error_type=dataframe_error_type,
    )
    if contract.require_site_key_index:
        if contract.prefer_analysis_ready_index_diagnostics:
            enforce_analysis_ready_site_key_index(
                site_metadata.index,
                field_name=f"{field_name}.index",
                error_type=error_type,
            )
        else:
            require_site_key_index(
                site_metadata.index,
                field_name=f"{field_name}.index",
                error_type=error_type,
            )
        if expected_index is not None:
            expected_name = (
                expected_index_field_name
                if expected_index_field_name is not None
                else "expected site_key index"
            )
            if contract.prefer_analysis_ready_index_diagnostics:
                enforce_analysis_ready_site_key_index(
                    expected_index,
                    field_name=expected_name,
                    error_type=error_type,
                )
            else:
                require_site_key_index(
                    expected_index,
                    field_name=expected_name,
                    error_type=error_type,
                )
    if compare_raw_site_key_column_before_decode:
        if contract.require_site_key_column_index_coherence:
            enforce_site_key_column_raw_matches_index(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
            )
            enforce_site_key_column(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
            )
        if contract.require_unique_site_key:
            enforce_unique_site_key_identity(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
            )
    else:
        if (
            contract.check_site_key_column_index_before_uniqueness
            and contract.require_site_key_column_index_coherence
        ):
            enforce_site_key_column_matches_index(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
            )
        if contract.require_unique_site_key:
            enforce_unique_site_key_identity(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
            )
        if (
            not contract.check_site_key_column_index_before_uniqueness
            and contract.require_site_key_column_index_coherence
        ):
            enforce_site_key_column_matches_index(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
            )
    if contract.require_display_id:
        enforce_display_id_column(
            site_metadata=site_metadata,
            field_name=field_name,
            error_type=error_type,
        )
    if expected_index is not None:
        require_exact_index_match(
            left=site_metadata.index,
            right=expected_index,
            left_name=f"{field_name}.index",
            right_name=(
                expected_index_field_name
                if expected_index_field_name is not None
                else "expected site_key index"
            ),
            error_type=dataframe_error_type,
        )
    if contract.require_protein_context:
        require_columns(
            site_metadata,
            field_name=field_name,
            required_columns=PROTEIN_CONTEXT_COLUMNS,
            error_type=dataframe_error_type,
        )
        if contract.require_site_key_metadata_coherence:
            enforce_site_key_matches_metadata(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
            )
    if contract.sequence_context is SequenceContextRequirement.PRESENT:
        require_columns(
            site_metadata,
            field_name=field_name,
            required_columns=(SITE_SEQUENCE_COLUMN,),
            error_type=dataframe_error_type,
        )
        if contract.require_non_empty_sequence_context:
            require_non_empty_string_column(
                site_metadata,
                field_name=field_name,
                column_name=SITE_SEQUENCE_COLUMN,
                error_type=dataframe_error_type,
            )
    if contract.sequence_context is SequenceContextRequirement.CENTRED:
        from phospy.validation.datasets.site_metadata import (
            enforce_centred_site_sequence_context,
        )

        enforce_centred_site_sequence_context(
            site_metadata=site_metadata,
            field_name=field_name,
            workflow_name=workflow_name or contract.contract_id,
            error_type=error_type,
            allow_gapped_sequence_context=allow_gapped_sequence_context,
            allow_unknown_site_residue=allow_opaque_site_values,
        )


def enforce_required_identity_text_columns(
    *,
    table: pd.DataFrame,
    field_name: str,
    columns: tuple[str, ...],
    error_type: type[ErrorType],
) -> None:
    """Require identity columns to contain non-empty strings."""

    dataframe_error_type = cast(ValidationErrorType, error_type)
    for column_name in columns:
        require_non_empty_string_column(
            table,
            field_name=field_name,
            column_name=column_name,
            error_type=dataframe_error_type,
        )


def enforce_analysis_ready_site_key_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[ErrorType],
    preview_limit: int = 5,
) -> pd.Index:
    """Require an analysis-ready public row index to be encoded site_key identity."""

    if index.name != _SITE_KEY_INDEX_NAME:
        if _looks_like_display_site_index(index):
            raise error_type(
                f"{field_name} is display-indexed direct construction; "
                "analysis-ready phosphosite row identity must use encoded "
                "site_key values with index.name='site_key', not GENE;SITE; labels"
            )
        raise error_type(
            f"{field_name} must be named 'site_key' for analysis-ready "
            "phosphosite row identity"
        )
    if _looks_like_display_site_index(index):
        raise error_type(
            f"{field_name} is display-indexed direct construction; analysis-ready "
            "phosphosite row identity must use encoded site_key values, not "
            "GENE;SITE; labels"
        )
    try:
        require_site_key_index(
            index,
            field_name=field_name,
            error_type=error_type,
            require_unique=False,
        )
    except error_type as exc:
        raise error_type(
            f"{field_name} must contain valid PhosPy site_key values; {exc}"
        ) from exc
    if index.is_unique:
        return index
    duplicate_values = _duplicate_label_values(index)
    preview = ", ".join(repr(value) for value in duplicate_values[:preview_limit])
    suffix = "" if len(duplicate_values) <= preview_limit else " ..."
    raise error_type(
        f"{field_name} must contain unique site_key values; "
        f"duplicate_site_key_values=[{preview}{suffix}]"
    )


def enforce_site_key_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "site_key",
) -> pd.Series:
    """Require one present site_key column with decodable encoded key values."""

    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    values = pd.Series(site_metadata[column_name], dtype="object")
    normalized: list[str] = []
    for row_id, raw_value in values.items():
        key = decode_site_key(
            raw_value,
            field_name=f"{field_name}.{column_name}[{row_id!r}]",
            error_type=error_type,
        )
        normalized.append(encode_site_key(key))
    return pd.Series(
        normalized,
        index=pd.Index(site_metadata.index),
        name=column_name,
        dtype="object",
    )


def enforce_display_id_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "display_id",
) -> pd.Series:
    """Require one present display_id column with recommended site identifiers."""

    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    return canonicalize_site_series(
        pd.Series(site_metadata[column_name], dtype="object"),
        field_name=f"{field_name}.{column_name}",
        error_type=error_type,
    )


def enforce_unique_site_key_identity(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str = "site_key",
    preview_limit: int = 5,
) -> pd.Series:
    """Require unique site_key identity and strict mixed-isoform consistency."""

    site_keys = enforce_site_key_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=site_key_column,
    )
    duplicate_mask = site_keys.duplicated(keep=False)
    if bool(duplicate_mask.any()):
        duplicate_values = list(
            dict.fromkeys(site_keys.loc[duplicate_mask].astype(str).tolist())
        )
        preview = ", ".join(repr(value) for value in duplicate_values[:preview_limit])
        suffix = "" if len(duplicate_values) <= preview_limit else " ..."
        raise error_type(
            f"{field_name}.{site_key_column} must be unique; "
            f"duplicate_values=[{preview}{suffix}]"
        )

    require_no_mixed_site_key_isoform_scope(
        site_keys=site_keys,
        field_name=f"{field_name}.{site_key_column}",
        error_type=error_type,
        preview_limit=preview_limit,
    )
    return site_keys


def enforce_site_key_column_matches_index(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str = "site_key",
    preview_limit: int = 5,
) -> pd.Series:
    """Require site_metadata.site_key values to exactly equal the row index."""

    site_keys = enforce_site_key_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=site_key_column,
    )
    raw_site_key_values = [
        str(value)
        for value in pd.Series(site_metadata[site_key_column], dtype="object").tolist()
    ]
    index_values = [str(value) for value in site_metadata.index.tolist()]
    mismatches: list[str] = []
    for position, (index_value, column_value) in enumerate(
        zip(index_values, raw_site_key_values, strict=True)
    ):
        if index_value == column_value:
            continue
        mismatches.append(
            f"position {position}: index={index_value!r}:"
            f"{site_key_column}={column_value!r}"
        )
    if not mismatches:
        return site_keys
    preview = ", ".join(mismatches[:preview_limit])
    suffix = "" if len(mismatches) <= preview_limit else " ..."
    raise error_type(
        f"{field_name}.{site_key_column} must exactly match {field_name}.index "
        f"values; mismatches=[{preview}{suffix}]"
    )


def enforce_site_key_column_raw_matches_index(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str = SITE_KEY_COLUMN,
    preview_limit: int = 5,
) -> None:
    """Require raw site_key column strings to exactly equal the row index."""

    if site_key_column not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {site_key_column}")
    site_key_values = [
        str(value)
        for value in pd.Series(site_metadata[site_key_column], dtype="object").tolist()
    ]
    index_values = [str(value) for value in site_metadata.index.tolist()]
    mismatches = [
        index_value
        for index_value, site_key in zip(index_values, site_key_values, strict=True)
        if index_value != site_key
    ]
    if not mismatches:
        return
    preview = ", ".join(repr(value) for value in mismatches[:preview_limit])
    suffix = "" if len(mismatches) <= preview_limit else " ..."
    raise error_type(
        f"{field_name}.{site_key_column} must exactly match {field_name}.index; "
        f"mismatched_labels={preview}{suffix}"
    )


def enforce_site_key_matches_metadata(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str = "site_key",
    isoform_column: str = "isoform_id",
    preview_limit: int = 5,
) -> pd.Series:
    """Require encoded site_key values to match metadata-derived protein keys."""

    site_keys = enforce_site_key_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=site_key_column,
    )

    mismatches: list[str] = []
    for row_id, encoded_site_key in site_keys.items():
        decoded_key = decode_site_key(
            encoded_site_key,
            field_name=f"{field_name}.{site_key_column}[{row_id!r}]",
            error_type=error_type,
        )
        metadata_key = build_protein_scoped_site_key(
            organism=_required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="organism",
                field_name=field_name,
                error_type=error_type,
            ),
            protein_namespace=_required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="protein_namespace",
                field_name=field_name,
                error_type=error_type,
            ),
            protein_identifier=_required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="protein_identifier",
                field_name=field_name,
                error_type=error_type,
            ),
            residue=_resolve_row_residue(
                site_metadata=site_metadata,
                row_id=row_id,
                field_name=field_name,
                error_type=error_type,
            ),
            position=_resolve_row_position(
                site_metadata=site_metadata,
                row_id=row_id,
                field_name=field_name,
                error_type=error_type,
            ),
            isoform_id=_optional_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name=isoform_column,
                field_name=field_name,
                error_type=error_type,
            ),
            field_name=f"{field_name}[{row_id!r}]",
            error_type=error_type,
        )
        if decoded_key != metadata_key:
            mismatches.append(
                f"{row_id!r}:observed={encoded_site_key!r}:"
                f"expected={encode_site_key(metadata_key)!r}"
            )

    if mismatches:
        preview = ", ".join(mismatches[:preview_limit])
        suffix = "" if len(mismatches) <= preview_limit else " ..."
        raise error_type(
            f"{field_name}.{site_key_column} must match metadata-derived "
            f"ProteinScopedPhosphositeKey values; mismatches=[{preview}{suffix}]"
        )
    return site_keys


def enforce_site_key_index(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str = "site_key",
    preview_limit: int = 5,
) -> None:
    """Require site_metadata index labels to match site_key values."""

    site_keys = enforce_site_key_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=site_key_column,
    )
    mismatches: list[str] = []
    for row_id, encoded_site_key in site_keys.items():
        index_value = row_id
        if not isinstance(index_value, str):
            mismatches.append(f"{row_id!r}:index={index_value!r}")
            continue
        if index_value != encoded_site_key:
            mismatches.append(
                f"{row_id!r}:index={index_value!r}:site_key={encoded_site_key!r}"
            )
    if not mismatches:
        return
    preview = ", ".join(mismatches[:preview_limit])
    suffix = "" if len(mismatches) <= preview_limit else " ..."
    raise error_type(
        f"{field_name}.index must match {field_name}.{site_key_column} when "
        f"enforced; mismatches=[{preview}{suffix}]"
    )


def enforce_result_identity_metadata_coherence(
    *,
    table: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
) -> None:
    """Require public result identity metadata to agree with encoded site_key."""

    optional_key_columns = {
        "organism": "organism",
        "protein_namespace": "protein_namespace",
        "protein_identifier": "protein_identifier",
    }
    for row_position, row_label in enumerate(table.index.tolist()):
        encoded_site_key = table.at[row_label, SITE_KEY_COLUMN]
        decoded_key = decode_site_key(
            encoded_site_key,
            field_name=f"{field_name}.{SITE_KEY_COLUMN}[{row_label!r}]",
            error_type=error_type,
        )
        encoded_site = f"{decoded_key.residue}{decoded_key.position}"
        row_site_value = table.at[row_label, "site"]
        parsed_site = try_parse_site_token(row_site_value)
        if parsed_site is None:
            raise error_type(
                _identity_incoherence_message(
                    field_name=field_name,
                    row_position=row_position,
                    row_label=row_label,
                    detail=(
                        f"row metadata site is {row_site_value!r}; "
                        "site must use strict 'S/T/Y<position>' tokens"
                    ),
                )
            )
        row_site = f"{parsed_site.residue}{parsed_site.position}"
        if row_site != encoded_site:
            raise error_type(
                _identity_incoherence_message(
                    field_name=field_name,
                    row_position=row_position,
                    row_label=row_label,
                    detail=(
                        f"site_key encodes {encoded_site} but row metadata "
                        f"site is {row_site_value!r}"
                    ),
                )
            )

        expected_display_id = canonicalize_site_components(
            table.at[row_label, "gene_symbol"],
            table.at[row_label, "site"],
            field_name=f"{field_name}[{row_label!r}].gene_symbol_site",
            error_type=error_type,
        )
        observed_display_id = canonicalize_site_identifier(
            table.at[row_label, DISPLAY_ID_COLUMN],
            field_name=f"{field_name}[{row_label!r}].display_id",
            error_type=error_type,
        )
        if observed_display_id != expected_display_id:
            raise error_type(
                _identity_incoherence_message(
                    field_name=field_name,
                    row_position=row_position,
                    row_label=row_label,
                    detail=(
                        "display_id does not match gene_symbol + site; "
                        f"expected {expected_display_id!r} but row metadata "
                        f"display_id is {table.at[row_label, DISPLAY_ID_COLUMN]!r}"
                    ),
                )
            )

        for column_name, key_attribute in optional_key_columns.items():
            if column_name not in table.columns:
                continue
            encoded_value = str(getattr(decoded_key, key_attribute))
            row_value = _required_present_identity_text(
                table.at[row_label, column_name],
                field_name=field_name,
                row_position=row_position,
                row_label=row_label,
                column_name=column_name,
                error_type=error_type,
            )
            if row_value == encoded_value:
                continue
            raise error_type(
                _identity_incoherence_message(
                    field_name=field_name,
                    row_position=row_position,
                    row_label=row_label,
                    detail=(
                        f"{column_name} is incoherent: site_key encodes "
                        f"{encoded_value!r} but row metadata {column_name} "
                        f"is {row_value!r}"
                    ),
                )
            )


def _required_text_value(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    column_name: str,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    value = site_metadata.at[row_id, column_name]
    if not isinstance(value, str):
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a non-empty string"
        )
    token = value.strip()
    if token == "":
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a non-empty string"
        )
    return token


def _optional_text_value(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    column_name: str,
    field_name: str,
    error_type: type[ErrorType],
) -> str | None:
    if column_name not in site_metadata.columns:
        return None
    value = site_metadata.at[row_id, column_name]
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a string when provided"
        )
    token = value.strip()
    if token == "":
        return None
    return token


def _resolve_row_residue(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    explicit_residue = (
        _optional_text_value(
            site_metadata=site_metadata,
            row_id=row_id,
            column_name="residue",
            field_name=field_name,
            error_type=error_type,
        )
        if "residue" in site_metadata.columns
        else None
    )
    parsed_site = _parse_row_site_token(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    if explicit_residue is None:
        if parsed_site is None:
            raise error_type(
                f"{field_name}[{row_id!r}] requires residue metadata or strict site "
                "token parsing to derive ProteinScopedPhosphositeKey"
            )
        return parsed_site.residue
    token = explicit_residue.upper()
    if len(token) != 1 or token not in {"S", "T", "Y"}:
        raise error_type(
            f"{field_name}[{row_id!r}].residue must be one of 'S', 'T', or 'Y'"
        )
    if parsed_site is not None and parsed_site.residue != token:
        raise error_type(
            f"{field_name}[{row_id!r}] has inconsistent residue metadata; "
            f"site_token={parsed_site.residue!r}, residue_column={token!r}"
        )
    return token


def _resolve_row_position(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> int:
    explicit_position = _resolve_explicit_position(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    parsed_site = _parse_row_site_token(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    if explicit_position is None:
        if parsed_site is None:
            raise error_type(
                f"{field_name}[{row_id!r}] requires position metadata or strict site "
                "token parsing to derive ProteinScopedPhosphositeKey"
            )
        return int(parsed_site.position)
    if parsed_site is not None and explicit_position != int(parsed_site.position):
        raise error_type(
            f"{field_name}[{row_id!r}] has inconsistent position metadata; "
            f"site_token={parsed_site.position}, position_column={explicit_position}"
        )
    return explicit_position


def _resolve_explicit_position(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> int | None:
    for column_name in _SITE_POSITION_CANDIDATE_COLUMNS:
        if column_name not in site_metadata.columns:
            continue
        raw_value = site_metadata.at[row_id, column_name]
        return require_positive_integer_position(
            raw_value,
            field_name=f"{field_name}[{row_id!r}].{column_name}",
            error_type=error_type,
        )
    return None


def _parse_row_site_token(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> ParsedSiteToken | None:
    if "site" not in site_metadata.columns:
        return None
    parsed = try_parse_site_token(site_metadata.at[row_id, "site"])
    if parsed is not None:
        return parsed
    raw_value = site_metadata.at[row_id, "site"]
    if _is_missing(raw_value):
        return None
    if isinstance(raw_value, str) and raw_value.strip() == "":
        return None
    raise error_type(
        f"{field_name}[{row_id!r}].site must use strict 'S/T/Y<position>' tokens"
    )


def _is_missing(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _looks_like_display_site_index(index: pd.Index) -> bool:
    values = index.tolist()
    if not values:
        return False
    if not all(isinstance(value, str) for value in values):
        return False
    if not any(";" in value for value in values):
        return False
    try:
        canonicalize_site_series(
            pd.Series(values, dtype="object"),
            field_name="analysis-ready row identity",
            error_type=ValueError,
        )
    except ValueError:
        return False
    return True


def _duplicate_label_values(index: pd.Index) -> list[object]:
    values = cast(list[object], index.tolist())
    duplicate_values: list[object] = []
    for value in values:
        if values.count(value) <= 1:
            continue
        if any(existing == value for existing in duplicate_values):
            continue
        duplicate_values.append(value)
    return duplicate_values


def _required_present_identity_text(
    value: object,
    *,
    field_name: str,
    row_position: int,
    row_label: object,
    column_name: str,
    error_type: type[ErrorType],
) -> str:
    if _is_missing(value):
        raise error_type(
            _identity_incoherence_message(
                field_name=field_name,
                row_position=row_position,
                row_label=row_label,
                detail=f"{column_name} must be a non-empty string when present",
            )
        )
    if not isinstance(value, str) or value.strip() == "":
        raise error_type(
            _identity_incoherence_message(
                field_name=field_name,
                row_position=row_position,
                row_label=row_label,
                detail=f"{column_name} must be a non-empty string when present",
            )
        )
    return value.strip()


def _identity_incoherence_message(
    *,
    field_name: str,
    row_position: int,
    row_label: object,
    detail: str,
) -> str:
    return (
        "Differential result identity metadata is inconsistent with site_key "
        f"at row {row_position} ({row_label!r}) in {field_name}: {detail}"
    )


__all__ = [
    "ANALYSIS_READY_DATASET_BASE_IDENTITY_CONTRACT",
    "ANALYSIS_READY_DATASET_IDENTITY_CONTRACT",
    "ANALYSIS_READY_IDENTITY_COLUMNS",
    "BASE_IDENTITY_COLUMNS",
    "DISPLAY_CONTEXT_COLUMNS",
    "DISPLAY_ID_COLUMN",
    "PhosphositeIdentityContract",
    "PROTEIN_CONTEXT_COLUMNS",
    "RESULT_IDENTITY_COLUMNS",
    "RESULT_TABLE_IDENTITY_CONTRACT",
    "SITE_KEY_COLUMN",
    "SITE_SEQUENCE_COLUMN",
    "SequenceContextRequirement",
    "WORKFLOW_INPUT_IDENTITY_CONTRACT",
    "WORKFLOW_PROTEIN_CONTEXT_IDENTITY_CONTRACT",
    "WORKFLOW_PROTEIN_SEQUENCE_CONTEXT_IDENTITY_CONTRACT",
    "WORKFLOW_SEQUENCE_CONTEXT_IDENTITY_CONTRACT",
    "enforce_analysis_ready_site_key_index",
    "enforce_display_id_column",
    "enforce_phosphosite_identity_contract",
    "enforce_required_identity_text_columns",
    "enforce_result_identity_metadata_coherence",
    "enforce_site_key_column",
    "enforce_site_key_column_matches_index",
    "enforce_site_key_column_raw_matches_index",
    "enforce_site_key_index",
    "enforce_site_key_matches_metadata",
    "enforce_unique_site_key_identity",
]
