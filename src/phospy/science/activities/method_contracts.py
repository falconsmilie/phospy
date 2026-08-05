"""Method-owned quantitative contracts for kinase activity methods."""

from __future__ import annotations

from phospy.science.activities.models import (
    KSEA_ZSCORE_ACTIVITY_METHOD,
    SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD,
    SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD,
)
from phospy.science.activities.semantics import (
    ActivityProfileAxis,
    ActivityQuantitativeSemantics,
)
from phospy.science.configs.kinase import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
    KINASE_ACTIVITY_METHODS,
    KinaseActivityMethod,
)
from phospy.science.quantitative_method_contracts import (
    MethodQuantitativeInputContract,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)


def simplified_weighted_substrate_activity_input_contract() -> (
    MethodQuantitativeInputContract
):
    """Return the quantitative contract for weighted substrate activity."""

    return MethodQuantitativeInputContract(
        method_id=SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD.activity_method_id,
        accepted_scales=(IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2),
        accepted_meanings=(
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
        ),
        accepted_activity_profile_axes=(
            ActivityProfileAxis.SAMPLE,
            ActivityProfileAxis.CONDITION_SUMMARY,
        ),
        accepted_activity_quantitative_semantics=(
            ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE,
            ActivityQuantitativeSemantics.CONDITION_SUMMARY_ABUNDANCE,
        ),
        required_centring=(
            "No method-level centring; activity values are weighted means on "
            "the declared input scale."
        ),
        required_standardisation=(
            "No automatic standardisation; linear and log2 abundance summaries "
            "have different meanings."
        ),
        missing_value_treatment=(
            "Missing substrate values are ignored per profile when computing "
            "weighted and thresholded means; no imputation is performed."
        ),
        profile_axis_requirements=(
            "Columns must represent sample-level abundance or explicit "
            "condition-summary abundance profiles."
        ),
        statistical_interpretation=(
            "Heuristic substrate-supported weighted mean; not a statistical "
            "enrichment test and not causal kinase activity proof."
        ),
        p_value_interpretation=None,
    )


def ksea_zscore_activity_input_contract() -> MethodQuantitativeInputContract:
    """Return the quantitative contract for KSEA-style z-score activity."""

    return MethodQuantitativeInputContract(
        method_id=KSEA_ZSCORE_ACTIVITY_METHOD.activity_method_id,
        accepted_scales=(IntensityScaleKind.LOG2,),
        accepted_meanings=(
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
            QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
            QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        ),
        accepted_activity_profile_axes=(
            ActivityProfileAxis.SAMPLE,
            ActivityProfileAxis.CONTRAST,
            ActivityProfileAxis.EFFECT,
        ),
        accepted_activity_quantitative_semantics=(
            ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE,
            ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE,
            ActivityQuantitativeSemantics.STANDARDISED_EFFECT,
        ),
        required_centring=(
            "Uses log2 sample, total-corrected, contrast, or effect profiles "
            "as declared by the dataset; no centring is performed in the method."
        ),
        required_standardisation=(
            "Requires log2 abundance, log2 total-corrected ratio, log2 contrast "
            "fold-change, or pre-standardised effect semantics; linear raw "
            "abundance is rejected."
        ),
        missing_value_treatment=(
            "Finite values define per-profile substrate and background sets; "
            "missing values are omitted from those calculations without imputation."
        ),
        profile_axis_requirements=(
            "Columns must represent log-scale sample profiles, contrasts, or "
            "standardised effect profiles; linear raw samples are rejected."
        ),
        statistical_interpretation=(
            "Unweighted substrate-set z-score enrichment over declared log-scale "
            "sample, contrast, or effect values with background variance checks."
        ),
        p_value_interpretation=(
            "Two-sided normal-approximation p-values for computed z-scores; "
            "available only when typed substrate-membership provenance declares "
            "the membership independent of the tested quantitative matrix. "
            "Eligible p-values use Benjamini-Hochberg q-value adjustment per "
            "profile when enabled; adaptive membership reports descriptive "
            "z-scores with p/q unavailable."
        ),
    )


def ssgsea_substrate_enrichment_activity_input_contract() -> (
    MethodQuantitativeInputContract
):
    """Return the quantitative contract for ssGSEA-style activity."""

    return MethodQuantitativeInputContract(
        method_id=SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD.activity_method_id,
        accepted_scales=(IntensityScaleKind.LOG2,),
        accepted_meanings=(
            QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
            QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        ),
        accepted_activity_profile_axes=(
            ActivityProfileAxis.CONTRAST,
            ActivityProfileAxis.EFFECT,
        ),
        accepted_activity_quantitative_semantics=(
            ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE,
            ActivityQuantitativeSemantics.STANDARDISED_EFFECT,
        ),
        required_centring=(
            "Uses ranked contrast/effect values supplied by the caller; no "
            "centring is performed inside the method."
        ),
        required_standardisation=(
            "Requires log2 contrast fold-change or pre-standardised effect "
            "semantics; raw abundance is rejected."
        ),
        missing_value_treatment=(
            "Only finite effect values enter the ranked background; missing "
            "values are omitted without imputation."
        ),
        profile_axis_requirements=(
            "Columns must represent contrasts or standardised effect profiles, "
            "not raw samples."
        ),
        statistical_interpretation=(
            "Rank-walk substrate-set enrichment summary over ordered effect "
            "values. Equal-valued finite sites are handled inside the method as "
            "tie blocks using the documented block-expectation policy, not row "
            "order or lexical site labels. Not PTM-SEA parity and not causal "
            "kinase activity proof."
        ),
        p_value_interpretation=(
            "No p-values are produced unless seeded permutations are requested; "
            "permutation p-values are two-sided empirical substrate-label "
            "permutation p-values, with Benjamini-Hochberg q-values per profile "
            "when enabled."
        ),
    )


def kinase_activity_method_quantitative_input_contract(
    method: KinaseActivityMethod | str,
) -> MethodQuantitativeInputContract:
    """Return the method-owned quantitative contract for an activity method."""

    method_text = str(method)
    if method_text == KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY:
        return simplified_weighted_substrate_activity_input_contract()
    if method_text == KINASE_ACTIVITY_METHOD_KSEA_ZSCORE:
        return ksea_zscore_activity_input_contract()
    if method_text == KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT:
        return ssgsea_substrate_enrichment_activity_input_contract()
    allowed = ", ".join(sorted(KINASE_ACTIVITY_METHODS))
    raise ValueError(f"kinase activity method must be one of: {allowed}")


def all_kinase_activity_method_quantitative_contracts() -> tuple[
    MethodQuantitativeInputContract,
    ...,
]:
    """Return documentation-ready quantitative contracts for activity methods."""

    return tuple(
        kinase_activity_method_quantitative_input_contract(method)
        for method in sorted(KINASE_ACTIVITY_METHODS)
    )


__all__ = [
    "all_kinase_activity_method_quantitative_contracts",
    "kinase_activity_method_quantitative_input_contract",
    "ksea_zscore_activity_input_contract",
    "simplified_weighted_substrate_activity_input_contract",
    "ssgsea_substrate_enrichment_activity_input_contract",
]
