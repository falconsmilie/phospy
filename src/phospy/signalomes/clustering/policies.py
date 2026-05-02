"""Shared clustering policy constants and type aliases."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

from phospy.scientific_policies import ScientificPolicyId, ScientificPolicyRecord

# Performance contracts for module-count selection scoring:
# - At or below `MAX_FULL_CORRELATION_SITE_COUNT`, candidate scoring computes a
#   full site-by-site correlation matrix.
# - Above this threshold, candidate scoring uses sampled within-cluster
#   correlations with at most `MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER` sites
#   per cluster.
#
# These thresholds only control how module-selection scores are computed. They do
# not change the input scoring matrix, selected output table schema, or whether
# approximation use is surfaced in diagnostics (`diagnostics.reason`).
MAX_FULL_CORRELATION_SITE_COUNT = 2000
MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER = 256
NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE = 1e-12

SIGNALOME_CLUSTERING_SCORING_MODE_AUTO = "auto"
SIGNALOME_CLUSTERING_SCORING_MODE_EXACT = "exact"
SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE = "approximate"
SignalomeClusteringScoringMode = Literal["auto", "exact", "approximate"]

SIGNALOME_TREE_ENGINE_EXACT = "exact"
SignalomeTreeEngine = Literal["exact"]

SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS = "column_median_imputation_with_zero_for_all_missing_columns"
SignalomeClusteringMissingValuePolicy = Literal[
    "column_median_imputation_with_zero_for_all_missing_columns"
]
# Signalome clustering-matrix preparation policy:
# - non-finite values are normalised to missing before imputation
# - partially missing columns are imputed with the column median
# - fully missing columns are imputed with 0.0
# - imputation is used for clustering distance/tree construction inputs only
#   (output tables remain the original workflow outputs and do not expose
#   imputed clustering values)
SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO = (
    "clustering_distance_and_tree_construction_only"
)
SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES = (
    False
)

SIGNALOME_CANDIDATE_SCORING_POLICY_FULL = "full"
SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED = "sampled"
SignalomeCandidateScoringPolicy = Literal["full", "sampled"]


@dataclass(frozen=True, slots=True)
class SignalomeCandidateScoringPolicyDefinition:
    """Versioned scientific policy for candidate module-count scoring behavior."""

    name: str
    version: str
    parameters: Mapping[str, object]
    description: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "parameters",
            MappingProxyType(
                {str(key): value for key, value in self.parameters.items()}
            ),
        )

    @property
    def record(self) -> ScientificPolicyRecord:
        return ScientificPolicyRecord(
            id=ScientificPolicyId.SIGNALOME_CANDIDATE_SCORING,
            name=self.name,
            version=self.version,
            description=self.description,
            parameters=self.parameters,
            assumptions=(
                "Candidate-scoring policy only changes module-count evaluation, not "
                "final exact tree construction semantics.",
            ),
            output_scale="Candidate module-count support summaries.",
            quantitative_meaning="relative_module_candidate_support",
        )


SIGNALOME_CANDIDATE_SCORING_FULL_POLICY = SignalomeCandidateScoringPolicyDefinition(
    name="signalome_candidate_scoring_full_v1",
    version="1",
    parameters={
        "candidate_scoring_policy": SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        "correlation_mode": "full_site_by_site_matrix",
    },
    description=(
        "Evaluate candidate module counts using full within-cluster correlation "
        "calculations."
    ),
)

SIGNALOME_CANDIDATE_SCORING_SAMPLED_POLICY = SignalomeCandidateScoringPolicyDefinition(
    name="signalome_candidate_scoring_sampled_v1",
    version="1",
    parameters={
        "candidate_scoring_policy": SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        "correlation_mode": "sampled_within_cluster_correlations",
        "sampling_method": "deterministic_uniform_without_replacement",
    },
    description=(
        "Evaluate candidate module counts with deterministic sampled "
        "within-cluster correlations."
    ),
)


def resolve_candidate_scoring_policy_definition(
    *,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy,
) -> SignalomeCandidateScoringPolicyDefinition:
    if candidate_scoring_policy == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED:
        return SIGNALOME_CANDIDATE_SCORING_SAMPLED_POLICY
    return SIGNALOME_CANDIDATE_SCORING_FULL_POLICY


SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED = "not_evaluated"
_CandidateScoringMode = SignalomeCandidateScoringPolicy | Literal["not_evaluated"]

SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD = (
    "deterministic_uniform_without_replacement"
)
SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY = (
    "order_invariant_seed_from_row_hashes_and_sample_size"
)
SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT = "explicit_module_count"
SIGNALOME_CANDIDATE_SCORING_APPLIES_TO = "candidate_module_count_evaluation_only"

SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE = "exact_cluster_tree"
SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE = "single_module_assignment"


__all__ = [
    "MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER",
    "MAX_FULL_CORRELATION_SITE_COUNT",
    "NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE",
    "SIGNALOME_CANDIDATE_SCORING_APPLIES_TO",
    "SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_FULL",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED",
    "SIGNALOME_CANDIDATE_SCORING_FULL_POLICY",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLED_POLICY",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY",
    "SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES",
    "SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE",
    "SIGNALOME_CLUSTERING_SCORING_MODE_AUTO",
    "SIGNALOME_CLUSTERING_SCORING_MODE_EXACT",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE",
    "SIGNALOME_TREE_ENGINE_EXACT",
    "SignalomeCandidateScoringPolicyDefinition",
    "SignalomeCandidateScoringPolicy",
    "SignalomeClusteringMissingValuePolicy",
    "SignalomeClusteringScoringMode",
    "SignalomeTreeEngine",
    "resolve_candidate_scoring_policy_definition",
    "_CandidateScoringMode",
]
