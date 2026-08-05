"""Immutable phosphosite identity contract records."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from phospy.science.sites.sequence_context import (
    WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT,
    SequenceContextContract,
)

REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE = "reference_context_unknown"
SITE_KEY_COLUMN = "site_key"
DISPLAY_ID_COLUMN = "display_id"
SITE_SEQUENCE_COLUMN = "site_sequence"
REFERENCE_CONTEXT_IDENTITY_FIELDS = (
    "organism",
    "protein_namespace",
    "source_name",
    "source_version",
    "proteome_version",
    "reference_table_sha256",
)
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
    sequence_context_contract: SequenceContextContract | None = None
    allow_duplicate_display_id: bool = True


@dataclass(frozen=True, slots=True)
class ReferenceContextCompatibilityWarning:
    """Typed warning returned when unknown reference context is explicitly allowed."""

    operation: str
    missing_contexts: tuple[str, ...]
    left_reference_context_id: str | None
    right_reference_context_id: str | None
    code: str = REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE
    severity: str = "warning"

    @property
    def message(self) -> str:
        missing = ", ".join(self.missing_contexts)
        noun = "contexts are" if len(self.missing_contexts) > 1 else "context is"
        return (
            "Reference-context compatibility could not be proven because "
            f"{missing} {noun} unknown for operation={self.operation!r}."
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "operation": self.operation,
            "missing_contexts": list(self.missing_contexts),
            "left_reference_context_id": self.left_reference_context_id,
            "right_reference_context_id": self.right_reference_context_id,
        }


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
    sequence_context_contract=WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT,
)
WORKFLOW_PROTEIN_SEQUENCE_CONTEXT_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="workflow_protein_sequence_context_identity",
    required_columns=BASE_IDENTITY_COLUMNS + PROTEIN_CONTEXT_COLUMNS,
    check_site_key_column_index_before_uniqueness=True,
    require_protein_context=True,
    require_site_key_metadata_coherence=True,
    sequence_context=SequenceContextRequirement.CENTRED,
    sequence_context_contract=WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT,
)
RESULT_TABLE_IDENTITY_CONTRACT = PhosphositeIdentityContract(
    contract_id="result_table_identity",
    required_columns=RESULT_IDENTITY_COLUMNS,
    require_protein_context=True,
    prefer_analysis_ready_index_diagnostics=False,
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
    "REFERENCE_CONTEXT_IDENTITY_FIELDS",
    "REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE",
    "RESULT_IDENTITY_COLUMNS",
    "RESULT_TABLE_IDENTITY_CONTRACT",
    "ReferenceContextCompatibilityWarning",
    "SITE_KEY_COLUMN",
    "SITE_SEQUENCE_COLUMN",
    "SequenceContextContract",
    "SequenceContextRequirement",
    "WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT",
    "WORKFLOW_INPUT_IDENTITY_CONTRACT",
    "WORKFLOW_PROTEIN_CONTEXT_IDENTITY_CONTRACT",
    "WORKFLOW_PROTEIN_SEQUENCE_CONTEXT_IDENTITY_CONTRACT",
    "WORKFLOW_SEQUENCE_CONTEXT_IDENTITY_CONTRACT",
]
