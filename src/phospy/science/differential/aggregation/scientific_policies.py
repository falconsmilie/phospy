"""Scientific policy records for typed peptide-to-site estimate combination."""

from __future__ import annotations

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)


def build_peptide_to_site_aggregation_policy(
    *,
    uncertainty_method: str,
    min_estimates_per_site: int,
    dependence_policy: str,
    multiple_testing_method: str,
    input_mapping_policies: tuple[str, ...],
) -> ScientificPolicyRecord:
    """Build deterministic provenance for peptide-to-site estimate combination."""

    return ScientificPolicyRecord(
        id=ScientificPolicyId.PEPTIDE_TO_SITE_AGGREGATION,
        name=f"peptide_to_site_{uncertainty_method}_v1",
        version="1",
        description=(
            "Records a typed peptide-to-site post-hoc differential estimate "
            "combination run. The preferred PhosPy-origin lane remains "
            "sample-intensity peptide evidence resolution before differential "
            "model fitting."
        ),
        parameters={
            "support_status": "supported_typed_estimate_combination_v1",
            "uncertainty_method": str(uncertainty_method),
            "min_estimates_per_site": int(min_estimates_per_site),
            "dependence_policy": str(dependence_policy),
            "multiple_testing_method": str(multiple_testing_method),
            "input_mapping_policy_count": int(len(input_mapping_policies)),
            "input_mapping_policies": "|".join(input_mapping_policies),
            "single_estimate_policy": "pass_through_original_t_and_p_value",
            "finite_df_t_to_z_policy": (
                "signed_two_sided_p_value_conversion_when_z_combination_is_used"
            ),
        },
        assumptions=(
            "Single-estimate site outputs are pass-through summaries and are "
            "not labelled as meta-analysis.",
            "Multi-estimate post-hoc combination is supported only for estimates "
            "from independent source experiments or runs.",
            "Same-experiment peptide estimates from the same samples are rejected "
            "because their dependence is not modelled by this lane.",
            "Stouffer-style z combination converts finite-degree-of-freedom "
            "t evidence through signed two-sided p-values rather than z=t.",
            "Multiple-testing correction is delegated to the configured shared "
            "correction method.",
        ),
        output_scale=(
            "Site-level post-hoc differential estimate summary with explicit "
            "uncertainty and dependence provenance."
        ),
        quantitative_meaning=(
            "typed_posthoc_peptide_differential_estimate_combination"
        ),
    )


__all__ = ["build_peptide_to_site_aggregation_policy"]
