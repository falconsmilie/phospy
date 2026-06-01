"""Shared phosphosite-identity contracts for workflow validators."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

import pandas as pd

from phospy.science.sites.identity import validate_no_conflicting_identity_collisions
from phospy.science.sites.validation import require_canonical_site_series
from phospy.validation.datasets.protein_scoped_site_identity import (
    enforce_display_id_column,
    enforce_site_key_column,
)
from phospy.validation.datasets.site_metadata import (
    enforce_centred_site_sequence_context,
    enforce_required_non_empty_string_column,
    enforce_site_identity_rows,
)

ErrorType = TypeVar("ErrorType", bound=Exception)


@dataclass(frozen=True, slots=True)
class WorkflowIdentityContract:
    """Workflow-specific phosphosite identity requirements."""

    workflow_name: str
    contract_id: str
    require_protein_identity: bool = False
    require_centred_sequence_context: bool = False
    allow_gapped_sequence_context: bool = False


DIFFERENTIAL_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="differential workflow request",
    contract_id="site_key_identity_minimum",
)
KINASE_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="kinase workflow request",
    contract_id="sty_site_identity_plus_sequence_context",
    require_centred_sequence_context=True,
    allow_gapped_sequence_context=True,
)
SIGNALOME_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="signalome workflow request",
    contract_id="protein_scoped_site_identity",
    require_protein_identity=True,
    require_centred_sequence_context=True,
    allow_gapped_sequence_context=True,
)


def enforce_workflow_site_identity_contract(
    *,
    site_metadata: pd.DataFrame,
    expected_index: pd.Index | None,
    field_name: str,
    contract: WorkflowIdentityContract,
    error_type: type[ErrorType],
    allow_opaque_site_values: bool,
) -> None:
    """Enforce one workflow's declared phosphosite identity contract."""

    try:
        has_site_key_column = "site_key" in site_metadata.columns
        if has_site_key_column:
            display_ids = enforce_display_id_column(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
                column_name="display_id",
            )
            site_keys_are_encoded = False
            try:
                enforce_site_key_column(
                    site_metadata=site_metadata,
                    field_name=field_name,
                    error_type=error_type,
                    column_name="site_key",
                )
                site_keys_are_encoded = True
            except error_type:
                require_canonical_site_series(
                    site_metadata.loc[:, "site_key"],
                    field_name=f"{field_name}.site_key",
                    error_type=error_type,
                )
            expected_index_is_encoded = _index_is_encoded_site_keys(expected_index)
            should_enforce_site_key_alignment = (
                not site_keys_are_encoded or expected_index_is_encoded
            )
            if should_enforce_site_key_alignment:
                _enforce_site_key_index_alignment(
                    site_metadata=site_metadata,
                    field_name=field_name,
                    error_type=error_type,
                    site_key_column="site_key",
                )
        elif "display_id" in site_metadata.columns:
            display_ids = enforce_display_id_column(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
                column_name="display_id",
            )
        else:
            display_ids = pd.Series(
                site_metadata.index.astype(str).tolist(),
                index=site_metadata.index.copy(),
                name="display_id",
            )
        site_metadata_for_display_checks = site_metadata.copy(deep=True)
        site_metadata_for_display_checks.index = pd.Index(
            display_ids.tolist(),
            name="display_id",
        )
        enforce_site_identity_rows(
            site_metadata=site_metadata_for_display_checks,
            field_name=field_name,
            error_type=error_type,
            allow_opaque_site_values=allow_opaque_site_values,
        )
        validate_no_conflicting_identity_collisions(
            site_metadata=site_metadata_for_display_checks,
            display_ids=display_ids,
            field_name=f"{field_name}.identity_collision_policy",
            error_type=error_type,
            allow_opaque_site_values=allow_opaque_site_values,
        )
        if contract.require_protein_identity:
            enforce_required_non_empty_string_column(
                site_metadata=site_metadata,
                field_name=field_name,
                workflow_name=contract.workflow_name,
                column_name="protein_id",
                error_type=error_type,
            )
        if contract.require_centred_sequence_context:
            enforce_centred_site_sequence_context(
                site_metadata=site_metadata,
                field_name=field_name,
                workflow_name=contract.workflow_name,
                error_type=error_type,
                allow_gapped_sequence_context=contract.allow_gapped_sequence_context,
                allow_unknown_site_residue=allow_opaque_site_values,
            )
    except error_type as exc:
        raise error_type(
            f"{contract.workflow_name} identity requirement failed "
            f"(contract={contract.contract_id}): {exc}"
        ) from exc


__all__ = [
    "DIFFERENTIAL_IDENTITY_CONTRACT",
    "KINASE_IDENTITY_CONTRACT",
    "SIGNALOME_IDENTITY_CONTRACT",
    "WorkflowIdentityContract",
    "enforce_workflow_site_identity_contract",
]


def _enforce_site_key_index_alignment(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_key_column: str,
    preview_limit: int = 5,
) -> None:
    site_keys = (
        site_metadata.loc[:, site_key_column].astype("string").str.strip().astype(str)
    )
    index_values = site_metadata.index.astype(str)
    mismatch_mask = site_keys.to_numpy() != index_values.to_numpy()
    if not bool(mismatch_mask.any()):
        return
    mismatches: list[str] = []
    for row_id, key_value, mismatched in zip(
        index_values.tolist(),
        site_keys.tolist(),
        mismatch_mask.tolist(),
        strict=False,
    ):
        if not bool(mismatched):
            continue
        mismatches.append(f"{row_id!r}:index={row_id!r}:site_key={key_value!r}")
    preview = ", ".join(mismatches[:preview_limit])
    suffix = "" if len(mismatches) <= preview_limit else " ..."
    raise error_type(
        f"{field_name}.index must match {field_name}.{site_key_column} when "
        f"enforced; mismatches=[{preview}{suffix}]"
    )


def _index_is_encoded_site_keys(index: pd.Index | None) -> bool:
    if index is None:
        return False
    values = index.astype(str).tolist()
    if not values:
        return False
    return all(value.startswith("phospy:v1|") for value in values)
