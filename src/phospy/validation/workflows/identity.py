"""Shared phosphosite-identity contracts for workflow validators."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

import pandas as pd

from phospy.errors.base import PhosPyError
from phospy.validation.identity_contracts import (
    WORKFLOW_PROTEIN_CONTEXT_IDENTITY_CONTRACT,
    WORKFLOW_PROTEIN_SEQUENCE_CONTEXT_IDENTITY_CONTRACT,
    WORKFLOW_SEQUENCE_CONTEXT_IDENTITY_CONTRACT,
    PhosphositeIdentityContract,
    enforce_phosphosite_identity_contract,
)

ErrorType = TypeVar("ErrorType", bound=PhosPyError)


@dataclass(frozen=True, slots=True)
class WorkflowIdentityContract:
    """Workflow-specific phosphosite identity requirements."""

    workflow_name: str
    contract_id: str
    identity_contract: PhosphositeIdentityContract
    allow_gapped_sequence_context: bool = False


DIFFERENTIAL_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="differential workflow request",
    contract_id="protein_scoped_site_identity",
    identity_contract=WORKFLOW_PROTEIN_CONTEXT_IDENTITY_CONTRACT,
)
KINASE_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="kinase workflow request",
    contract_id="sty_site_identity_plus_sequence_context",
    identity_contract=WORKFLOW_SEQUENCE_CONTEXT_IDENTITY_CONTRACT,
    allow_gapped_sequence_context=True,
)
SIGNALOME_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="signalome workflow request",
    contract_id="protein_scoped_site_identity",
    identity_contract=WORKFLOW_PROTEIN_SEQUENCE_CONTEXT_IDENTITY_CONTRACT,
    allow_gapped_sequence_context=True,
)


def enforce_workflow_site_identity_contract(
    *,
    site_metadata: pd.DataFrame,
    expected_index: pd.Index,
    expected_index_field_name: str,
    field_name: str,
    contract: WorkflowIdentityContract,
    error_type: type[ErrorType],
    allow_opaque_site_values: bool,
) -> None:
    """Enforce one workflow's declared phosphosite identity contract."""

    try:
        # Keep workflow ownership explicit while delegating implementation to
        # the shared contract layer. The composed helpers are:
        # enforce_analysis_ready_site_key_index(...)
        # enforce_site_key_column_matches_index(...)
        # enforce_site_key_matches_metadata(...)
        # enforce_display_id_column(...)
        # require_exact_index_match(...)
        # enforce_centred_site_sequence_context(...)
        enforce_phosphosite_identity_contract(
            site_metadata=site_metadata,
            expected_index=expected_index,
            expected_index_field_name=expected_index_field_name,
            field_name=field_name,
            contract=contract.identity_contract,
            error_type=error_type,
            workflow_name=contract.workflow_name,
            allow_opaque_site_values=allow_opaque_site_values,
            allow_gapped_sequence_context=contract.allow_gapped_sequence_context,
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
