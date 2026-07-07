"""Explicit input contracts for kinase workflow scoring modes."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.contracts.configs import (
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    KinaseScoringMode,
    normalize_kinase_scoring_mode,
)
from phospy.errors.workflows import PhosPyWorkflowError

KINASE_ANALYSIS_READY_SITE_METADATA_COLUMNS = (
    "site_key",
    "display_id",
    "organism",
    "protein_namespace",
    "protein_identifier",
    "gene_symbol",
    "site",
    "site_sequence",
)


@dataclass(frozen=True, slots=True)
class KinaseScoringModeInputContract:
    """Required workflow inputs for one canonical kinase scoring mode."""

    scoring_mode: KinaseScoringMode
    required_dataset_columns: tuple[str, ...]
    requires_site_sequence: bool
    requires_centered_sequence_context: bool
    requires_localisation_probability: bool
    requires_substrate_reference_overlap: bool
    requires_kinase_library_resource: bool
    requires_profile_construction: bool


def _contract(
    *,
    scoring_mode: KinaseScoringMode,
    requires_substrate_reference_overlap: bool,
    requires_kinase_library_resource: bool,
    requires_profile_construction: bool,
) -> KinaseScoringModeInputContract:
    return KinaseScoringModeInputContract(
        scoring_mode=scoring_mode,
        required_dataset_columns=KINASE_ANALYSIS_READY_SITE_METADATA_COLUMNS,
        requires_site_sequence=True,
        requires_centered_sequence_context=True,
        requires_localisation_probability=False,
        requires_substrate_reference_overlap=requires_substrate_reference_overlap,
        requires_kinase_library_resource=requires_kinase_library_resource,
        requires_profile_construction=requires_profile_construction,
    )


KINASE_SCORING_MODE_INPUT_CONTRACTS = {
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED: _contract(
        scoring_mode=KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=False,
        requires_profile_construction=True,
    ),
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF: _contract(
        scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=True,
        requires_profile_construction=True,
    ),
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY: _contract(
        scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
        requires_substrate_reference_overlap=False,
        requires_kinase_library_resource=True,
        requires_profile_construction=False,
    ),
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF: _contract(
        scoring_mode=KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=True,
        requires_profile_construction=True,
    ),
}


def kinase_scoring_mode_input_contract(
    scoring_mode: object,
) -> KinaseScoringModeInputContract:
    """Return the explicit input contract for a canonical scoring mode."""

    normalized = normalize_kinase_scoring_mode(scoring_mode)
    contract = KINASE_SCORING_MODE_INPUT_CONTRACTS.get(normalized)
    if contract is None:
        known_modes = ", ".join(sorted(KINASE_SCORING_MODE_INPUT_CONTRACTS))
        raise PhosPyWorkflowError(
            "kinase scoring mode input contract is not defined for "
            f"scoring_mode={normalized!r}; add an entry to "
            "KINASE_SCORING_MODE_INPUT_CONTRACTS before enabling this mode; "
            f"known_modes={known_modes}"
        )
    return contract


__all__ = [
    "KINASE_ANALYSIS_READY_SITE_METADATA_COLUMNS",
    "KINASE_SCORING_MODE_INPUT_CONTRACTS",
    "KinaseScoringModeInputContract",
    "kinase_scoring_mode_input_contract",
]
