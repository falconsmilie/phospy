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
from phospy.science.quantitative_method_contracts import (
    MethodQuantitativeInputContract,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)

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
    quantitative_contract: MethodQuantitativeInputContract
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
    quantitative_contract: MethodQuantitativeInputContract,
) -> KinaseScoringModeInputContract:
    return KinaseScoringModeInputContract(
        scoring_mode=scoring_mode,
        required_dataset_columns=KINASE_ANALYSIS_READY_SITE_METADATA_COLUMNS,
        quantitative_contract=quantitative_contract,
        requires_site_sequence=True,
        requires_centered_sequence_context=True,
        requires_localisation_probability=False,
        requires_substrate_reference_overlap=requires_substrate_reference_overlap,
        requires_kinase_library_resource=requires_kinase_library_resource,
        requires_profile_construction=requires_profile_construction,
    )


def _profile_quantitative_contract(
    *,
    scoring_mode: KinaseScoringMode,
    allow_mixed_total_protein_quantitative_meaning: bool = False,
) -> MethodQuantitativeInputContract:
    accepted_meanings = [
        QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
    ]
    if allow_mixed_total_protein_quantitative_meaning:
        accepted_meanings.append(
            QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE
        )
    return MethodQuantitativeInputContract(
        method_id=f"kinase_scoring.{scoring_mode}",
        accepted_scales=(IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2),
        accepted_meanings=tuple(accepted_meanings),
        required_centring=(
            "Requires centered phosphosite sequence context; does not center "
            "quantitative values during scoring."
        ),
        required_standardisation=(
            "No automatic quantitative standardisation; values are consumed on "
            "the declared input scale."
        ),
        missing_value_treatment=(
            "Profile construction follows scoring_config.profile_missing_value_strategy; "
            "missing values are never imputed by the scoring method and no "
            "method-level imputation is performed."
        ),
        profile_axis_requirements=(
            "Rows are protein-scoped site_key phosphosites; columns are aligned "
            "sample/profile abundance or total-corrected quantitative observations "
            "used for profile support."
        ),
        statistical_interpretation=(
            "Profile-derived relative support scores are within-run evidence "
            "summaries over the declared abundance/profile axis; linear and "
            "log2 inputs are scale-sensitive and are not numerically "
            "interchangeable."
        ),
        p_value_interpretation=None,
    )


def _motif_only_quantitative_contract(
    *,
    scoring_mode: KinaseScoringMode,
) -> MethodQuantitativeInputContract:
    return MethodQuantitativeInputContract(
        method_id=f"kinase_scoring.{scoring_mode}",
        accepted_scales=(IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2),
        accepted_meanings=(
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
            QuantitativeMeaning.UNKNOWN,
        ),
        required_centring=(
            "Requires centered phosphosite sequence context; quantitative "
            "centring is not applicable because phospho values are not consumed "
            "by motif-only scoring."
        ),
        required_standardisation=(
            "No quantitative standardisation is required or performed for "
            "motif-only scoring."
        ),
        missing_value_treatment=(
            "Phospho missing values are not read by motif-only scoring; no "
            "missing-value transformation or imputation is performed."
        ),
        profile_axis_requirements=(
            "Rows are protein-scoped site_key phosphosites with centered "
            "sequence context; quantitative columns are not used for motif-only "
            "score calculation."
        ),
        statistical_interpretation=(
            "Scores are sequence-motif support scores from the supplied Kinase "
            "Library-style resource, not abundance-profile statistics."
        ),
        p_value_interpretation=None,
        quantitative_input_required=False,
        scale_sensitivity=(
            "The quantitative matrix is not consumed by motif-only scoring; the "
            "resolved dataset scale is recorded for audit but does not affect "
            "motif scores."
        ),
    )


KINASE_SCORING_MODE_INPUT_CONTRACTS = {
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED: _contract(
        scoring_mode=KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=False,
        requires_profile_construction=True,
        quantitative_contract=_profile_quantitative_contract(
            scoring_mode=KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
        ),
    ),
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF: _contract(
        scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=True,
        requires_profile_construction=True,
        quantitative_contract=_profile_quantitative_contract(
            scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF
        ),
    ),
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY: _contract(
        scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
        requires_substrate_reference_overlap=False,
        requires_kinase_library_resource=True,
        requires_profile_construction=False,
        quantitative_contract=_motif_only_quantitative_contract(
            scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY
        ),
    ),
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF: _contract(
        scoring_mode=KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
        requires_substrate_reference_overlap=True,
        requires_kinase_library_resource=True,
        requires_profile_construction=True,
        quantitative_contract=_profile_quantitative_contract(
            scoring_mode=KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF
        ),
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


def kinase_scoring_method_quantitative_input_contract(
    scoring_mode: object,
    *,
    allow_mixed_total_protein_quantitative_meaning: bool = False,
) -> MethodQuantitativeInputContract:
    """Return the method-owned quantitative contract for one scoring mode."""

    normalized = normalize_kinase_scoring_mode(scoring_mode)
    if normalized == KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY:
        return _motif_only_quantitative_contract(scoring_mode=normalized)
    # Ensure the scoring mode itself is known before returning a profile contract.
    kinase_scoring_mode_input_contract(normalized)
    return _profile_quantitative_contract(
        scoring_mode=normalized,
        allow_mixed_total_protein_quantitative_meaning=(
            allow_mixed_total_protein_quantitative_meaning
        ),
    )


def all_kinase_scoring_method_quantitative_contracts() -> tuple[
    MethodQuantitativeInputContract,
    ...,
]:
    """Return documentation-ready quantitative contracts for scoring modes."""

    return tuple(
        kinase_scoring_method_quantitative_input_contract(scoring_mode)
        for scoring_mode in sorted(KINASE_SCORING_MODE_INPUT_CONTRACTS)
    )


__all__ = [
    "KINASE_ANALYSIS_READY_SITE_METADATA_COLUMNS",
    "KINASE_SCORING_MODE_INPUT_CONTRACTS",
    "KinaseScoringModeInputContract",
    "all_kinase_scoring_method_quantitative_contracts",
    "kinase_scoring_mode_input_contract",
    "kinase_scoring_method_quantitative_input_contract",
]
