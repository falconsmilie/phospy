"""Shared phosphosite-identity contracts for workflow validators."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

import pandas as pd

from phospy.science.sites.identity import validate_no_conflicting_identity_collisions
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
    contract_id="display_site_identity_minimum",
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
    field_name: str,
    contract: WorkflowIdentityContract,
    error_type: type[ErrorType],
    allow_opaque_site_values: bool,
) -> None:
    """Enforce one workflow's declared phosphosite identity contract."""

    try:
        enforce_site_identity_rows(
            site_metadata=site_metadata,
            field_name=field_name,
            error_type=error_type,
            allow_opaque_site_values=allow_opaque_site_values,
        )
        validate_no_conflicting_identity_collisions(
            site_metadata=site_metadata,
            display_ids=pd.Series(
                site_metadata.index.tolist(),
                index=site_metadata.index.copy(),
                name="display_site_id",
                dtype="object",
            ),
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
