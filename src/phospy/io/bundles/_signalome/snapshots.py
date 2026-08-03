"""Config snapshot contract for signalome bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from phospy.contracts.configs import (
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT,
    SignalomeConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import require_mapping
from phospy.io.bundles._signalome.config import (
    signalome_config_from_payload,
)
from phospy.io.bundles._signalome.primitives import _parse_optional_int

if TYPE_CHECKING:
    from phospy.contracts.requests import SignalomeWorkflowRequest


@dataclass(frozen=True, slots=True)
class SignalomeWorkflowConfigSnapshot:
    """Serializable snapshot of the signalome workflow configuration."""

    signalome_config: SignalomeConfig
    network_min_paired_finite_observations_effective: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.signalome_config, SignalomeConfig):
            raise PhosPyInputError(
                "config snapshot signalome_config must be a SignalomeConfig"
            )
        if self.network_min_paired_finite_observations_effective is not None:
            object.__setattr__(
                self,
                "network_min_paired_finite_observations_effective",
                int(self.network_min_paired_finite_observations_effective),
            )

    @classmethod
    def from_request(
        cls, request: SignalomeWorkflowRequest
    ) -> SignalomeWorkflowConfigSnapshot:
        """Create a config snapshot from a workflow request."""

        from phospy.contracts.requests import SignalomeWorkflowRequest

        if not isinstance(request, SignalomeWorkflowRequest):
            raise PhosPyInputError(
                "config snapshot request must be a SignalomeWorkflowRequest"
            )
        return cls(
            signalome_config=request.config,
            network_min_paired_finite_observations_effective=(
                SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT
                if request.config.output.network_min_paired_finite_observations is None
                else int(request.config.output.network_min_paired_finite_observations)
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a manifest-safe JSON payload for this config snapshot."""

        return {
            "signalome_config": {
                "mode": str(self.signalome_config.mode),
                "scientific": {
                    "substrate_support_cutoff": float(
                        self.signalome_config.scientific.substrate_support_cutoff
                    ),
                    "assignment_policy": str(
                        self.signalome_config.scientific.assignment_policy
                    ),
                },
                "clustering": {
                    "module_count": (
                        None
                        if self.signalome_config.clustering.module_count is None
                        else int(self.signalome_config.clustering.module_count)
                    ),
                    "module_selection_primary_correlation_threshold": float(
                        self.signalome_config.clustering.module_selection_primary_correlation_threshold
                    ),
                    "module_selection_fallback_correlation_threshold": float(
                        self.signalome_config.clustering.module_selection_fallback_correlation_threshold
                    ),
                    "module_selection_max_clusters": int(
                        self.signalome_config.clustering.module_selection_max_clusters
                    ),
                    "candidate_scoring_policy": str(
                        self.signalome_config.clustering.candidate_scoring_policy
                    ),
                    "clustering_engine": str(
                        self.signalome_config.clustering.clustering_engine
                    ),
                },
                "validation": {
                    "score_preconditioning_policy": str(
                        self.signalome_config.validation.score_preconditioning_policy
                    ),
                    "allow_mixed_total_protein_quantitative_meaning": bool(
                        self.signalome_config.validation.allow_mixed_total_protein_quantitative_meaning
                    ),
                    "reference_context_compatibility_policy": str(
                        self.signalome_config.validation.reference_context_compatibility_policy
                    ),
                    "localisation_requirement": {
                        "require_present": bool(
                            self.signalome_config.validation.localisation_requirement.require_present
                        ),
                        "minimum_probability": (
                            None
                            if (
                                self.signalome_config.validation.localisation_requirement.minimum_probability
                                is None
                            )
                            else float(
                                self.signalome_config.validation.localisation_requirement.minimum_probability
                            )
                        ),
                    },
                },
                "output": {
                    "network_correlation_threshold": float(
                        self.signalome_config.output.network_correlation_threshold
                    ),
                    "network_policy": str(self.signalome_config.output.network_policy),
                    "network_min_paired_finite_observations": (
                        None
                        if (
                            self.signalome_config.output.network_min_paired_finite_observations
                            is None
                        )
                        else int(
                            self.signalome_config.output.network_min_paired_finite_observations
                        )
                    ),
                    "network_min_paired_finite_observations_effective": (
                        None
                        if self.network_min_paired_finite_observations_effective is None
                        else int(self.network_min_paired_finite_observations_effective)
                    ),
                },
                "performance": {
                    "max_exact_tree_sites": int(
                        self.signalome_config.performance.max_exact_tree_sites
                    ),
                    "max_full_candidate_scoring_sites": int(
                        self.signalome_config.performance.max_full_candidate_scoring_sites
                    ),
                },
            }
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
        *,
        effective_network_min_paired_finite_observations: int | None = None,
    ) -> SignalomeWorkflowConfigSnapshot:
        """Create a config snapshot from a decoded JSON payload."""

        scope = "config snapshot"
        signalome_payload = require_mapping(
            payload.get("signalome_config"),
            field_name=f"{scope}.signalome_config",
        )
        signalome_config = signalome_config_from_payload(
            signalome_payload,
            scope=scope,
            allow_legacy_network_minimum=True,
        )
        output_payload = require_mapping(
            signalome_payload.get("output"),
            field_name=f"{scope}.signalome_config.output",
        )
        parsed_effective = _parse_optional_int(
            output_payload.get("network_min_paired_finite_observations_effective"),
            field_name=(
                f"{scope}.signalome_config.output."
                "network_min_paired_finite_observations_effective"
            ),
        )
        return cls(
            signalome_config=signalome_config,
            network_min_paired_finite_observations_effective=(
                effective_network_min_paired_finite_observations
                if effective_network_min_paired_finite_observations is not None
                else parsed_effective
            ),
        )
