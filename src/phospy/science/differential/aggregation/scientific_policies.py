"""Experimental policy records for unsupported peptide-to-site aggregation."""

from __future__ import annotations

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)


def build_peptide_to_site_aggregation_policy(
    *,
    strategy: str,
    min_peptides_per_site: int,
    missing_variance_policy: str,
    stouffer_weighting: str,
    random_effect_tau2_floor: float,
    compatibility_mode_warning: bool,
) -> ScientificPolicyRecord:
    assumptions = [
        "This record is emitted only by an internal/experimental compatibility "
        "route, not by a supported PhosPy public differential API.",
        "Aggregation consumes same-experiment peptide-level differential model "
        "outputs without refitting a site-level statistical model.",
        "The current fixed-effect, random-effect, inverse-variance, Stouffer, "
        "and minimum-p strategies are not supported for production site-level "
        "inference while the statistical model is being corrected.",
        "Minimum-p compatibility mode is intended only for historical "
        "reproducibility and can bias significance.",
    ]
    if compatibility_mode_warning:
        assumptions.append(
            "Compatibility warning: minimum peptide p-value selection treats "
            "peptides as competitors rather than combined evidence."
        )
    return ScientificPolicyRecord(
        id=ScientificPolicyId.PEPTIDE_TO_SITE_AGGREGATION,
        name=f"peptide_to_site_aggregation_{strategy}_experimental_internal_v1",
        version="1",
        description=(
            "Records an internal experimental post-hoc peptide-to-site "
            "differential aggregation run. This is not a supported site-level "
            "inferential lane."
        ),
        parameters={
            "support_status": "experimental_internal_compatibility_only",
            "strategy": str(strategy),
            "min_peptides_per_site": int(min_peptides_per_site),
            "missing_variance_policy": str(missing_variance_policy),
            "stouffer_weighting": str(stouffer_weighting),
            "random_effect_tau2_floor": float(random_effect_tau2_floor),
            "compatibility_mode_warning": bool(compatibility_mode_warning),
        },
        assumptions=tuple(assumptions),
        output_scale=(
            "Experimental post-hoc site summary generated from peptide-level "
            "statistics; not supported site-level uncertainty aggregation."
        ),
        quantitative_meaning="experimental_internal_posthoc_peptide_summary",
    )


__all__ = ["build_peptide_to_site_aggregation_policy"]
