"""Shared phosphosite-identity contracts for workflow validators."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

import pandas as pd

from phospy.errors.base import PhosPyError
from phospy.validation.common.dataframes import (
    require_columns,
    require_exact_index_match,
)
from phospy.validation.datasets.protein_scoped_site_identity import (
    enforce_analysis_ready_site_key_index,
    enforce_display_id_column,
    enforce_site_key_column_matches_index,
    enforce_site_key_matches_metadata,
)
from phospy.validation.datasets.site_metadata import (
    enforce_centred_site_sequence_context,
)

ErrorType = TypeVar("ErrorType", bound=PhosPyError)


@dataclass(frozen=True, slots=True)
class WorkflowIdentityContract:
    """Workflow-specific phosphosite identity requirements."""

    workflow_name: str
    contract_id: str
    # Requires dataset-level protein-scoped row identity metadata, not the
    # signalome-only site_metadata.protein_id grouping column.
    require_protein_identity: bool = False
    require_centred_sequence_context: bool = False
    allow_gapped_sequence_context: bool = False


DIFFERENTIAL_IDENTITY_CONTRACT = WorkflowIdentityContract(
    workflow_name="differential workflow request",
    contract_id="protein_scoped_site_identity",
    require_protein_identity=True,
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
    expected_index: pd.Index,
    expected_index_field_name: str,
    field_name: str,
    contract: WorkflowIdentityContract,
    error_type: type[ErrorType],
    allow_opaque_site_values: bool,
) -> None:
    """Enforce one workflow's declared phosphosite identity contract."""

    try:
        enforce_analysis_ready_site_key_index(
            expected_index,
            field_name=expected_index_field_name,
            error_type=error_type,
        )
        enforce_analysis_ready_site_key_index(
            site_metadata.index,
            field_name=f"{field_name}.index",
            error_type=error_type,
        )
        enforce_site_key_column_matches_index(
            site_metadata=site_metadata,
            field_name=field_name,
            error_type=error_type,
            site_key_column="site_key",
        )
        enforce_display_id_column(
            site_metadata=site_metadata,
            field_name=field_name,
            error_type=error_type,
            column_name="display_id",
        )
        require_exact_index_match(
            left=site_metadata.index,
            right=expected_index,
            left_name=f"{field_name}.index",
            right_name=expected_index_field_name,
            error_type=error_type,
        )
        if contract.require_protein_identity:
            require_columns(
                site_metadata,
                field_name=field_name,
                required_columns=(
                    "organism",
                    "protein_namespace",
                    "protein_identifier",
                    "site",
                ),
                error_type=error_type,
            )
            enforce_site_key_matches_metadata(
                site_metadata=site_metadata,
                field_name=field_name,
                error_type=error_type,
                site_key_column="site_key",
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
