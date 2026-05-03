"""Scientific wrappers for activity-stage computation."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.activities.methods.weighted_substrate_activity import (
    SimplifiedWeightedSubstrateActivityMethod,
)
from phospy.activities.models import KinaseActivityInputs, KinaseActivityResult
from phospy.scientific_policies import (
    ScientificPolicyRecord,
    build_simplified_weighted_substrate_activity_policy,
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


__all__ = [
    "SimplifiedWeightedSubstrateActivityPolicy",
    "compute_activity_from_inputs",
]
