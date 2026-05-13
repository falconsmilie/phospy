"""Scientific policy records for signalome clustering behavior."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)


@dataclass(frozen=True, slots=True)
class SignalomeMissingValueClusteringPolicy:
    """Executable metadata policy for clustering missing-value behavior."""

    missing_value_policy: str
    applies_to: str
    imputed_values_exposed_in_output_tables: bool

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_signalome_missing_value_clustering_policy(
            missing_value_policy=self.missing_value_policy,
            applies_to=self.applies_to,
            imputed_values_exposed_in_output_tables=(
                self.imputed_values_exposed_in_output_tables
            ),
        )


def build_signalome_missing_value_clustering_policy(
    *,
    missing_value_policy: str,
    applies_to: str,
    imputed_values_exposed_in_output_tables: bool,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_MISSING_VALUE_CLUSTERING,
        name="Signalome Missing-Value Clustering Policy",
        version="1",
        description=(
            "Normalizes non-finite clustering inputs to missing values and imputes "
            "missing clustering cells for distance/tree construction."
        ),
        parameters={
            "missing_value_policy": str(missing_value_policy),
            "applies_to": str(applies_to),
            "imputed_values_exposed_in_output_tables": bool(
                imputed_values_exposed_in_output_tables
            ),
            "partial_missingness_handling": "column_median_imputation",
            "fully_missing_column_handling": "impute_zero",
        },
        assumptions=(
            "Imputation is used for clustering internals and may influence module "
            "selection and assignment outcomes.",
            "Output signalome tables do not expose the imputed clustering matrix.",
        ),
        output_scale=(
            "Prepared clustering values used for distance calculations and tree "
            "construction."
        ),
        quantitative_meaning="clustering_preconditioned_support_matrix",
    )


def build_signalome_module_candidate_score_policy(
    *,
    requested_policy: str,
    candidate_scoring_policy: str,
    candidate_scoring_mode: str,
    max_exact_tree_sites: int | None,
    max_full_candidate_scoring_sites: int,
    candidate_scoring_evaluated: bool,
    candidate_scoring_skip_reason: str | None,
    candidate_scoring_scope: str = "candidate_module_count_evaluation_only",
    tree_generation_mode: str = "full_exact_tree_construction",
    tree_generation_is_approximate: bool = False,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_MODULE_CANDIDATE_SCORE,
        name="Signalome Module Candidate Score",
        version="1",
        description=(
            "Ranks candidate module counts using within-cluster median "
            "correlation summaries over downstream kinase-score profiles. "
            "Candidate scoring policy does not alter tree-generation exactness "
            "in the current implementation."
        ),
        parameters={
            "requested_policy": str(requested_policy),
            "candidate_scoring_policy": str(candidate_scoring_policy),
            "candidate_scoring_mode": str(candidate_scoring_mode),
            "candidate_scoring_scope": str(candidate_scoring_scope),
            "max_exact_tree_sites": max_exact_tree_sites,
            "max_full_candidate_scoring_sites": int(max_full_candidate_scoring_sites),
            "tree_generation_mode": str(tree_generation_mode),
            "tree_generation_is_approximate": bool(tree_generation_is_approximate),
            "candidate_scoring_evaluated": bool(candidate_scoring_evaluated),
            "candidate_scoring_skip_reason": (
                None
                if candidate_scoring_skip_reason is None
                else str(candidate_scoring_skip_reason)
            ),
        },
        assumptions=(
            "Candidate quality is summarized by within-cluster correlation coherence.",
            "Degenerate or undefined profiles are excluded or tracked via diagnostics.",
            "Selected module count depends on thresholds and candidate-scoring policy.",
            "Tree generation remains exact and is guarded separately from candidate scoring.",
        ),
        output_scale=(
            "Candidate module-count support scores; higher values indicate stronger "
            "within-cluster profile coherence."
        ),
        quantitative_meaning="relative_module_candidate_support",
    )


PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.PROTEIN_MODULE_FROM_SITE_MEMBERSHIP,
    name="protein_module_from_site_membership_v1",
    version="1",
    description=(
        "Derives protein-level module IDs by grouping proteins with matching "
        "site-cluster membership patterns."
    ),
    parameters={
        "membership_representation": "binary site-cluster incidence vector",
        "module_id_assignment": "first-seen pattern order",
    },
    assumptions=(
        "Site-cluster membership captures protein-level signaling context.",
        "Proteins with identical site-cluster incidence vectors are grouped into "
        "the same module.",
    ),
    output_scale="Integer module IDs at the protein level.",
    quantitative_meaning="protein_module_membership_label",
)


__all__ = [
    "PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY",
    "SignalomeMissingValueClusteringPolicy",
    "build_signalome_missing_value_clustering_policy",
    "build_signalome_module_candidate_score_policy",
]
