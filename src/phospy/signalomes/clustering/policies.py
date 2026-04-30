"""Shared clustering policy constants and type aliases."""

from __future__ import annotations

from typing import Literal

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

SIGNALOME_CANDIDATE_SCORING_POLICY_FULL = "full"
SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED = "sampled"
SignalomeCandidateScoringPolicy = Literal["full", "sampled"]

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
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY",
    "SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT",
    "SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE",
    "SIGNALOME_CLUSTERING_SCORING_MODE_AUTO",
    "SIGNALOME_CLUSTERING_SCORING_MODE_EXACT",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE",
    "SIGNALOME_TREE_ENGINE_EXACT",
    "SignalomeCandidateScoringPolicy",
    "SignalomeClusteringScoringMode",
    "SignalomeTreeEngine",
    "_CandidateScoringMode",
]
