"""Config snapshot contract for signalome bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from phospy.api.configs import SignalomeConfig
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import require_mapping
from phospy.io.bundles._signalome.config import (
    signalome_config_from_payload,
)

if TYPE_CHECKING:
    from phospy.api.requests import SignalomeWorkflowRequest


@dataclass(frozen=True, slots=True)
class SignalomeWorkflowConfigSnapshot:
    """Serializable snapshot of the signalome workflow configuration."""

    signalome_config: SignalomeConfig

    @classmethod
    def from_request(
        cls, request: SignalomeWorkflowRequest
    ) -> SignalomeWorkflowConfigSnapshot:
        """Create a config snapshot from a workflow request."""

        from phospy.api.requests import SignalomeWorkflowRequest

        if not isinstance(request, SignalomeWorkflowRequest):
            raise PhosPyInputError(
                "config snapshot request must be a SignalomeWorkflowRequest"
            )
        return cls(signalome_config=request.config)

    def to_payload(self) -> dict[str, object]:
        """Return a manifest-safe JSON payload for this config snapshot."""

        return {
            "signalome_config": {
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
        cls, payload: Mapping[str, object]
    ) -> SignalomeWorkflowConfigSnapshot:
        """Create a config snapshot from a decoded JSON payload."""

        scope = "config snapshot"
        signalome_payload = require_mapping(
            payload.get("signalome_config"),
            field_name=f"{scope}.signalome_config",
        )
        return cls(
            signalome_config=signalome_config_from_payload(
                signalome_payload,
                scope=scope,
            )
        )
