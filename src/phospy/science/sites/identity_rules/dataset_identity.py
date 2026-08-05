"""Dataset and workflow phosphosite identity validation rules."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, TypeVar, cast

import pandas as pd

from phospy.frames.validation import (
    ValidationErrorType,
    require_columns,
    require_exact_index_match,
    require_non_empty_string_column,
)
from phospy.science.sites.identifiers import canonicalize_site_series
from phospy.science.sites.identity_rules.contracts import (
    PROTEIN_CONTEXT_COLUMNS,
    SITE_KEY_COLUMN,
    SITE_SEQUENCE_COLUMN,
    PhosphositeIdentityContract,
    SequenceContextRequirement,
)
from phospy.science.sites.identity_rules.parsing import (
    looks_like_display_site_index,
    optional_text_value,
    required_text_value,
    resolve_row_position,
    resolve_row_residue,
)
from phospy.science.sites.sequence_context import (
    SequenceContextContract,
    enforce_centred_site_sequence_context,
    enforce_site_sequence_context_contract,
)
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    decode_site_key,
    encode_site_key,
)
from phospy.science.sites.validation import (
    require_no_mixed_site_key_isoform_scope,
    require_site_key_index,
)

ErrorType = TypeVar("ErrorType", bound=Exception)
_SITE_KEY_INDEX_NAME = "site_key"


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
    sequence_context_contract: SequenceContextContract | None = None,
    scoring_mode: str | None = None,
    sequence_source_by_site: Mapping[Any, object] | None = None,
    allow_unknown_site_residue: bool | None = None,
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
        resolved_sequence_contract = (
            sequence_context_contract or contract.sequence_context_contract
        )
        if resolved_sequence_contract is None:
            enforce_centred_site_sequence_context(
                site_metadata=site_metadata,
                field_name=field_name,
                workflow_name=workflow_name or contract.contract_id,
                error_type=error_type,
                allow_gapped_sequence_context=allow_gapped_sequence_context,
                allow_unknown_site_residue=allow_opaque_site_values,
            )
            return
        enforce_site_sequence_context_contract(
            site_metadata=site_metadata,
            field_name=field_name,
            workflow_name=workflow_name or contract.contract_id,
            error_type=error_type,
            scoring_mode=scoring_mode,
            contract=resolved_sequence_contract,
            sequence_source_by_site=sequence_source_by_site,
            allow_unknown_site_residue=(
                allow_opaque_site_values
                if allow_unknown_site_residue is None
                else allow_unknown_site_residue
            ),
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
        if looks_like_display_site_index(index):
            raise error_type(
                f"{field_name} is display-indexed direct construction; "
                "analysis-ready phosphosite row identity must use encoded "
                "site_key values with index.name='site_key', not GENE;SITE; labels"
            )
        raise error_type(
            f"{field_name} must be named 'site_key' for analysis-ready "
            "phosphosite row identity"
        )
    if looks_like_display_site_index(index):
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
            organism=required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="organism",
                field_name=field_name,
                error_type=error_type,
            ),
            protein_namespace=required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="protein_namespace",
                field_name=field_name,
                error_type=error_type,
            ),
            protein_identifier=required_text_value(
                site_metadata=site_metadata,
                row_id=row_id,
                column_name="protein_identifier",
                field_name=field_name,
                error_type=error_type,
            ),
            residue=resolve_row_residue(
                site_metadata=site_metadata,
                row_id=row_id,
                field_name=field_name,
                error_type=error_type,
            ),
            position=resolve_row_position(
                site_metadata=site_metadata,
                row_id=row_id,
                field_name=field_name,
                error_type=error_type,
            ),
            isoform_id=optional_text_value(
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


__all__ = [
    "enforce_analysis_ready_site_key_index",
    "enforce_display_id_column",
    "enforce_phosphosite_identity_contract",
    "enforce_required_identity_text_columns",
    "enforce_site_key_column",
    "enforce_site_key_column_matches_index",
    "enforce_site_key_column_raw_matches_index",
    "enforce_site_key_index",
    "enforce_site_key_matches_metadata",
    "enforce_unique_site_key_identity",
]
