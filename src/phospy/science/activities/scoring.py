"""Scientific wrappers for activity-stage computation."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyRecord,
)
from phospy.science.activities.methods.ssgsea_substrate_enrichment import (
    SsgseaSubstrateEnrichmentActivityMethod,
)
from phospy.science.activities.methods.weighted_substrate_activity import (
    SimplifiedWeightedSubstrateActivityMethod,
)
from phospy.science.activities.models import KinaseActivityInputs, KinaseActivityResult
from phospy.science.activities.scientific_policies import (
    build_simplified_weighted_substrate_activity_policy,
    build_ssgsea_substrate_enrichment_activity_policy,
)


def compute_activity_from_inputs(inputs: KinaseActivityInputs) -> KinaseActivityResult:
    """Compute legacy weighted activity outputs from validated inputs."""

    policy = SimplifiedWeightedSubstrateActivityPolicy(
        threshold=float(inputs.threshold),
        min_substrates=int(inputs.min_substrates),
        top_n_substrates=int(inputs.top_n_substrates),
    )
    return policy.compute(inputs=inputs)


@dataclass(frozen=True, slots=True)
class SimplifiedWeightedSubstrateActivityPolicy:
    """Executable policy wrapper for the weighted heuristic activity method."""

    threshold: float
    min_substrates: int
    top_n_substrates: int

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_simplified_weighted_substrate_activity_policy(
            threshold=float(self.threshold),
            min_substrates=int(self.min_substrates),
            top_n_substrates=int(self.top_n_substrates),
        )

    def compute(self, *, inputs: KinaseActivityInputs) -> KinaseActivityResult:
        method = SimplifiedWeightedSubstrateActivityMethod(
            threshold=float(self.threshold),
            min_substrates=int(self.min_substrates),
            top_n_substrates=int(self.top_n_substrates),
        )
        return method.run(inputs)


@dataclass(frozen=True, slots=True)
class SsgseaSubstrateEnrichmentActivityPolicy:
    """Executable policy wrapper for ssGSEA-style substrate enrichment."""

    min_substrates: int
    ranking_direction: str
    permutation_count: int = 0
    random_seed: int | None = 0
    adjust_p_values: bool = True

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_ssgsea_substrate_enrichment_activity_policy(
            min_substrates=int(self.min_substrates),
            ranking_direction=str(self.ranking_direction),
            permutation_count=int(self.permutation_count),
            random_seed=self.random_seed,
            adjust_p_values=bool(self.adjust_p_values),
            q_value_method=(
                "benjamini_hochberg"
                if int(self.permutation_count) > 0 and bool(self.adjust_p_values)
                else None
            ),
        )

    def compute(
        self,
        *,
        effect_matrix,
        kinase_substrate_membership,
    ) -> KinaseActivityResult:
        method = SsgseaSubstrateEnrichmentActivityMethod(
            min_substrates=int(self.min_substrates),
            ranking_direction=str(self.ranking_direction),
            permutation_count=int(self.permutation_count),
            random_seed=self.random_seed,
            adjust_p_values=bool(self.adjust_p_values),
        )
        return method.run(
            effect_matrix=effect_matrix,
            kinase_substrate_membership=kinase_substrate_membership,
        )


__all__ = [
    "SimplifiedWeightedSubstrateActivityPolicy",
    "SsgseaSubstrateEnrichmentActivityPolicy",
    "compute_activity_from_inputs",
]
