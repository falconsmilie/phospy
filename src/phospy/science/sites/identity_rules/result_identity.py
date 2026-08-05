"""Result-table phosphosite identity validation rules."""

from __future__ import annotations

from typing import TypeVar

import pandas as pd

from phospy.science.references.models import Organism
from phospy.science.sites.identifiers import (
    canonicalize_site_components,
    canonicalize_site_identifier,
)
from phospy.science.sites.identity_rules.contracts import (
    DISPLAY_ID_COLUMN,
    RESULT_IDENTITY_COLUMNS,
    RESULT_TABLE_IDENTITY_CONTRACT,
    SITE_KEY_COLUMN,
)
from phospy.science.sites.identity_rules.dataset_identity import (
    enforce_phosphosite_identity_contract,
    enforce_required_identity_text_columns,
)
from phospy.science.sites.identity_rules.messages import (
    identity_incoherence_message,
)
from phospy.science.sites.identity_rules.parsing import (
    is_missing,
    parse_site_token,
)
from phospy.science.sites.site_keys import decode_site_key

ErrorType = TypeVar("ErrorType", bound=Exception)


def enforce_result_table_identity_contract(
    *,
    table: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    context_label: str = "Identity metadata",
    identity_columns: tuple[str, ...] = RESULT_IDENTITY_COLUMNS,
) -> None:
    """Require a public result table to carry coherent site-key identity."""

    enforce_phosphosite_identity_contract(
        site_metadata=table,
        field_name=field_name,
        contract=RESULT_TABLE_IDENTITY_CONTRACT,
        error_type=error_type,
        compare_raw_site_key_column_before_decode=True,
    )
    enforce_required_identity_text_columns(
        table=table,
        field_name=field_name,
        columns=identity_columns,
        error_type=error_type,
    )
    enforce_result_identity_metadata_coherence(
        table=table,
        field_name=field_name,
        error_type=error_type,
        context_label=context_label,
    )


def enforce_result_identity_metadata_coherence(
    *,
    table: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    context_label: str = "Identity metadata",
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
        try:
            parsed_site = parse_site_token(
                row_site_value,
                field_name=f"{field_name}[{row_label!r}].site",
                error_type=error_type,
            )
        except error_type:
            raise error_type(
                identity_incoherence_message(
                    context_label=context_label,
                    field_name=field_name,
                    row_position=row_position,
                    row_label=row_label,
                    detail=(
                        f"row metadata site is {row_site_value!r}; "
                        "site must use strict 'S/T/Y<position>' tokens"
                    ),
                )
            ) from None
        row_site = f"{parsed_site.residue}{parsed_site.position}"
        if row_site != encoded_site:
            raise error_type(
                identity_incoherence_message(
                    context_label=context_label,
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
                identity_incoherence_message(
                    context_label=context_label,
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
            encoded_attribute = getattr(decoded_key, key_attribute)
            encoded_value = (
                encoded_attribute.value
                if isinstance(encoded_attribute, Organism)
                else str(encoded_attribute)
            )
            row_value = _required_present_identity_text(
                table.at[row_label, column_name],
                field_name=field_name,
                row_position=row_position,
                row_label=row_label,
                column_name=column_name,
                context_label=context_label,
                error_type=error_type,
            )
            if row_value == encoded_value:
                continue
            raise error_type(
                identity_incoherence_message(
                    context_label=context_label,
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


def _required_present_identity_text(
    value: object,
    *,
    field_name: str,
    row_position: int,
    row_label: object,
    column_name: str,
    error_type: type[ErrorType],
    context_label: str = "Identity metadata",
) -> str:
    if is_missing(value):
        raise error_type(
            identity_incoherence_message(
                context_label=context_label,
                field_name=field_name,
                row_position=row_position,
                row_label=row_label,
                detail=f"{column_name} must be a non-empty string when present",
            )
        )
    if not isinstance(value, str) or value.strip() == "":
        raise error_type(
            identity_incoherence_message(
                context_label=context_label,
                field_name=field_name,
                row_position=row_position,
                row_label=row_label,
                detail=f"{column_name} must be a non-empty string when present",
            )
        )
    return value.strip()


__all__ = [
    "enforce_result_identity_metadata_coherence",
    "enforce_result_table_identity_contract",
]
