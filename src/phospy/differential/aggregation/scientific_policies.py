"""Scientific policy records for peptide-to-site differential aggregation."""

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
        "Aggregation consumes peptide-level differential model outputs without "
        "refitting the upstream differential model.",
        "Site-level uncertainty is derived from peptide-level uncertainty statistics.",
        "Minimum-p compatibility mode is intended only for historical reproducibility "
        "and can bias significance.",
    ]
    if compatibility_mode_warning:
        assumptions.append(
            "Compatibility warning: minimum peptide p-value selection treats "
            "peptides as competitors rather than combined evidence."
        )
    return ScientificPolicyRecord(
        id=ScientificPolicyId.PEPTIDE_TO_SITE_AGGREGATION,
        name=f"peptide_to_site_aggregation_{strategy}_v1",
        version="1",
        description=(
            "Aggregates peptide-level differential statistics into site-level "
            "summaries with an explicit strategy."
        ),
        parameters={
            "strategy": str(strategy),
            "min_peptides_per_site": int(min_peptides_per_site),
            "missing_variance_policy": str(missing_variance_policy),
            "stouffer_weighting": str(stouffer_weighting),
            "random_effect_tau2_floor": float(random_effect_tau2_floor),
            "compatibility_mode_warning": bool(compatibility_mode_warning),
        },
        assumptions=tuple(assumptions),
        output_scale="Site-level log fold-change and uncertainty summaries.",
        quantitative_meaning="site_level_differential_summary",
    )


__all__ = ["build_peptide_to_site_aggregation_policy"]
