"""Stable scientific-policy identifiers and serializable metadata records."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType

from phospy.policy_models import ThresholdMode

ScientificPolicyParameter = str | int | float | bool | None


class ScientificPolicyId(str, Enum):
    """Stable identifiers for scientific scoring and derivation behavior."""

    PROFILE_CORRELATION_SHIFTED_UNIT = "profile_correlation_shifted_unit_v1"
    KINASE_PROFILE_SCORING = "kinase_profile_scoring_v1"
    MOTIF_PROFILE_RANK_FUSION = "motif_profile_rank_fusion_v1"
    CANDIDATE_SUBSTRATE_SELECTION = "candidate_substrate_selection_v1"
    SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY = "simplified_weighted_substrate_activity_v1"
    KSEA_ZSCORE_ACTIVITY = "ksea_zscore_activity_v1"
    SIGNALOME_MISSING_VALUE_CLUSTERING = "signalome_missing_value_clustering_v1"
    SIGNALOME_SCORE_PRECONDITIONING = "signalome_score_preconditioning_v1"
    PREPROCESSING_STAGE_ORDER = "preprocessing_stage_order_v1"
    SIGNALOME_MODULE_CANDIDATE_SCORE = "signalome_module_candidate_score_v1"
    PROTEIN_MODULE_FROM_SITE_MEMBERSHIP = "protein_module_from_site_membership_v1"
    DUPLICATE_SITE_RESOLUTION = "duplicate_site_resolution_v1"
    ADAPTIVE_PREDICTION_SAMPLING = "adaptive_prediction_sampling_v1"
    SIGNALOME_DOWNSTREAM_SCORE_SELECTION = "signalome_downstream_score_selection_v1"
    SIGNALOME_CANDIDATE_SCORING = "signalome_candidate_scoring_v1"
    SIGNALOME_ASSIGNMENT_POLICY = "signalome_assignment_policy_v1"
    SIGNALOME_NETWORK_POLICY = "signalome_network_policy_v1"
    PEPTIDE_TO_SITE_AGGREGATION = "peptide_to_site_aggregation_v1"


@dataclass(frozen=True, slots=True)
class ScientificPolicyRecord:
    """Serializable metadata for one scientific scoring/derivation policy."""

    id: ScientificPolicyId
    name: str
    version: str
    description: str
    parameters: Mapping[str, object]
    assumptions: tuple[str, ...]
    output_scale: str | None = None
    quantitative_meaning: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "parameters",
            MappingProxyType(
                {str(key): value for key, value in self.parameters.items()}
            ),
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "id": self.id.value,
            "name": self.name,
            "version": self.version,
            "description": self.description,
            "parameters": dict(self.parameters),
            "assumptions": list(self.assumptions),
            "output_scale": self.output_scale,
            "quantitative_meaning": self.quantitative_meaning,
        }

    @classmethod
    def from_payload(cls, payload: dict[str, object]) -> ScientificPolicyRecord:
        assumptions = payload.get("assumptions", ())
        if isinstance(assumptions, (list, tuple)):
            resolved_assumptions = tuple(str(value) for value in assumptions)
        else:
            resolved_assumptions = ()
        parameters_raw = payload.get("parameters", {})
        parameters: dict[str, object]
        if isinstance(parameters_raw, dict):
            parameters = {}
            for key, value in parameters_raw.items():
                if value is None or isinstance(value, (str, int, float, bool)):
                    parameters[str(key)] = value
                else:
                    parameters[str(key)] = str(value)
        else:
            parameters = {}
        output_scale = payload.get("output_scale")
        resolved_output_scale = None if output_scale is None else str(output_scale)
        quantitative_meaning = payload.get("quantitative_meaning")
        resolved_quantitative_meaning = (
            None if quantitative_meaning is None else str(quantitative_meaning)
        )
        return cls(
            id=ScientificPolicyId(str(payload.get("id"))),
            name=str(payload.get("name")),
            version=str(payload.get("version")),
            description=str(payload.get("description")),
            parameters=parameters,
            assumptions=resolved_assumptions,
            output_scale=resolved_output_scale,
            quantitative_meaning=resolved_quantitative_meaning,
        )


@dataclass(frozen=True, slots=True)
class KinaseProfileScoringPolicy:
    """Executable metadata policy for kinase profile scoring behavior."""

    profile_missing_value_strategy: str
    min_substrates_floor: int
    requested_min_substrates: int
    self_inclusion_behavior: str = "self_inclusion"
    leave_one_out_enabled: bool = False

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_kinase_profile_scoring_policy(
            profile_missing_value_strategy=self.profile_missing_value_strategy,
            min_substrates_floor=self.min_substrates_floor,
            requested_min_substrates=self.requested_min_substrates,
            self_inclusion_behavior=self.self_inclusion_behavior,
            leave_one_out_enabled=self.leave_one_out_enabled,
        )


@dataclass(frozen=True, slots=True)
class CandidateSubstrateSelectionPolicy:
    """Executable metadata policy for candidate substrate selection behavior."""

    top_k: int
    score_threshold: float
    inclusion: int
    threshold_operator: ThresholdMode = ThresholdMode.GREATER_THAN
    ranking_rule: str = "top_n_scores_per_kinase_then_threshold"
    site_restriction: str = "none"

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "threshold_operator",
            ThresholdMode.parse(
                self.threshold_operator,
                field_name=("candidate substrate selection policy threshold_operator"),
            ),
        )

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_candidate_substrate_selection_policy(
            top_k=self.top_k,
            score_threshold=self.score_threshold,
            inclusion=self.inclusion,
            threshold_operator=self.threshold_operator,
            ranking_rule=self.ranking_rule,
            site_restriction=self.site_restriction,
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


@dataclass(frozen=True, slots=True)
class ScorePreconditioningPolicy:
    """Executable metadata policy for score preconditioning behavior."""

    policy: str
    input_row_count: int
    dropped_all_missing_row_count: int
    retained_row_count: int
    row_retention_rule: str = "drop_rows_with_all_scores_missing"
    retained_partial_missing_rows: bool = True

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_score_preconditioning_policy(
            policy=self.policy,
            input_row_count=self.input_row_count,
            dropped_all_missing_row_count=self.dropped_all_missing_row_count,
            retained_row_count=self.retained_row_count,
            row_retention_rule=self.row_retention_rule,
            retained_partial_missing_rows=self.retained_partial_missing_rows,
        )


@dataclass(frozen=True, slots=True)
class PreprocessingStageOrderPolicy:
    """Executable metadata policy for preprocessing stage-order behavior."""

    configured_stage_order: tuple[str, ...]
    default_stage_order: tuple[str, ...]
    supported_stage_order: tuple[str, ...]

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_preprocessing_stage_order_policy(
            configured_stage_order=self.configured_stage_order,
            default_stage_order=self.default_stage_order,
            supported_stage_order=self.supported_stage_order,
        )


PROFILE_CORRELATION_SHIFTED_UNIT_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.PROFILE_CORRELATION_SHIFTED_UNIT,
    name="profile_correlation_v1",
    version="1",
    description=(
        "Transforms Pearson correlation from [-1, 1] to [0, 1] using (r + 1) / 2."
    ),
    parameters={
        "transform": "(r + 1) / 2",
        "clip_to_unit_interval": True,
        "preserve_undefined_as_nan": True,
    },
    assumptions=(
        "Higher positive correlation indicates stronger support.",
        "Negative correlation is treated as lower support, not explicit "
        "inhibitory evidence.",
        "Undefined correlations remain missing (NaN).",
    ),
    output_scale=(
        "Relative support score in [0, 1] where larger means stronger positive "
        "profile agreement."
    ),
    quantitative_meaning="relative_support_score",
)


def build_motif_profile_rank_fusion_policy(
    *,
    allow_profile_only_fallback: bool,
    emit_weights: bool,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.MOTIF_PROFILE_RANK_FUSION,
        name="rank_weighted_motif_profile_fusion_v1",
        version="1",
        description=(
            "Combines motif-frequency and profile-correlation scores using "
            "rank-derived logarithmic weights."
        ),
        parameters={
            "motif_weight_formula": "log(rank(motif_size)+1) / total_weight",
            "profile_weight_formula": "log(rank(profile_size)+1) / total_weight",
            "allow_profile_only_fallback": bool(allow_profile_only_fallback),
            "emit_weights": bool(emit_weights),
        },
        assumptions=(
            "Motif-library size and quantified-substrate count proxy evidence "
            "strength.",
            "When motif evidence is missing for a kinase/site, profile evidence can "
            "be propagated.",
            "Outputs are relative support scores and are not calibrated probabilities.",
        ),
        output_scale="Relative downstream support score for kinase-site ranking.",
        quantitative_meaning="relative_support_score",
    )


def build_kinase_profile_scoring_policy(
    *,
    profile_missing_value_strategy: str,
    min_substrates_floor: int,
    requested_min_substrates: int,
    self_inclusion_behavior: str = "self_inclusion",
    leave_one_out_enabled: bool = False,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.KINASE_PROFILE_SCORING,
        name="Kinase Profile Scoring Policy",
        version="1",
        description=(
            "Builds kinase reference profiles from quantified substrates and scores "
            "sites against those profiles with shifted Pearson correlation support."
        ),
        parameters={
            "profile_missing_value_strategy": str(profile_missing_value_strategy),
            "self_inclusion_behavior": str(self_inclusion_behavior),
            "leave_one_out_enabled": bool(leave_one_out_enabled),
            "min_substrates_floor": int(min_substrates_floor),
            "requested_min_substrates": int(requested_min_substrates),
        },
        assumptions=(
            "Profiles can include the same substrate site that is later scored when "
            "that site is present in the kinase profile definition.",
            "Leave-one-out profile recomputation is not applied in this policy.",
            "Profile missing-value strategy affects profile medians and can change "
            "site-level downstream support.",
        ),
        output_scale=(
            "Relative downstream support score in [0, 1] after shifted-correlation "
            "transformation."
        ),
        quantitative_meaning="relative_support_score",
    )


def build_candidate_substrate_selection_policy(
    *,
    top_k: int,
    score_threshold: float,
    inclusion: int,
    threshold_operator: ThresholdMode | str = ThresholdMode.GREATER_THAN,
    ranking_rule: str = "top_n_scores_per_kinase_then_threshold",
    site_restriction: str = "none",
) -> ScientificPolicyRecord:
    resolved_threshold_mode = ThresholdMode.parse(
        threshold_operator,
        field_name="candidate substrate selection policy threshold_operator",
    )
    return ScientificPolicyRecord(
        id=ScientificPolicyId.CANDIDATE_SUBSTRATE_SELECTION,
        name="Candidate Substrate Selection Policy",
        version="1",
        description=(
            "Selects per-kinase candidate substrate sites from downstream support "
            "scores using top-k ranking, threshold filtering, and minimum inclusion."
        ),
        parameters={
            "top_k": int(top_k),
            "score_threshold": float(score_threshold),
            "inclusion": int(inclusion),
            "threshold_operator": resolved_threshold_mode.value,
            "ranking_rule": str(ranking_rule),
            "site_restriction": str(site_restriction),
        },
        assumptions=(
            "Only finite scores are eligible for candidate selection.",
            "Threshold and inclusion rules jointly determine which kinases are "
            "considered to have usable candidate substrate support.",
            "Changing selection thresholds changes downstream kinase ranking and "
            "prediction outputs.",
        ),
        output_scale=(
            "Per-kinase candidate substrate sets for downstream ranking/prediction."
        ),
        quantitative_meaning="candidate_support_set",
    )


def build_simplified_weighted_substrate_activity_policy(
    *,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
        name="Simplified Weighted Substrate Activity",
        version="1",
        description=(
            "Computes weighted activity and thresholded substrate-mean activity "
            "from predicted substrate support."
        ),
        parameters={
            "threshold": float(threshold),
            "min_substrates": int(min_substrates),
            "top_n_substrates": int(top_n_substrates),
            "weighted_activity_rule": (
                "prediction-weighted mean over top-N predicted substrates"
            ),
            "thresholded_activity_rule": (
                "mean phospho over predicted substrates with "
                f"{ThresholdMode.GREATER_THAN_OR_EQUAL.value}"
            ),
        },
        assumptions=(
            "Predicted substrate support approximates kinase-substrate relevance.",
            "Higher weighted/thresholded values indicate stronger relative activity "
            "support in-run.",
            "This is not full KSEA-style statistical enrichment.",
        ),
        output_scale=(
            "Sample-by-kinase relative activity summaries (weighted mean and "
            "thresholded mean)."
        ),
        quantitative_meaning="relative_activity_support",
    )


def build_ksea_zscore_activity_policy(
    *,
    evidence_threshold: float,
    min_substrates: int,
    p_value_method: str,
    adjust_p_values: bool,
    q_value_method: str | None,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.KSEA_ZSCORE_ACTIVITY,
        name="ksea_zscore_activity_v1",
        version="1",
        description=(
            "Computes KSEA-style substrate-set enrichment activity z-scores using "
            "unweighted substrate membership after evidence thresholding."
        ),
        parameters={
            "evidence_threshold": float(evidence_threshold),
            "min_substrates": int(min_substrates),
            "membership_rule": (
                "finite_evidence "
                f"{_threshold_operator_token(ThresholdMode.GREATER_THAN_OR_EQUAL)} "
                "evidence_threshold"
            ),
            "weighting_rule": "unweighted_membership",
            "z_score_formula": "(mean_S - mean_U) * sqrt(n) / sd_U",
            "background_sd_ddof": 1,
            "p_value_method": str(p_value_method),
            "adjust_p_values": bool(adjust_p_values),
            "q_value_method": None if q_value_method is None else str(q_value_method),
        },
        assumptions=(
            "Substrate evidence contributes as binary membership after thresholding.",
            "Background phosphosite values define per-condition mean and sample variance.",
            "Scores with insufficient substrates or invalid background variance are not computable.",
            "KSEA z-scores are statistical enrichment summaries and are not PhosR-equivalent activity inference.",
        ),
        output_scale=(
            "Condition-by-kinase z-score substrate-set enrichment activity matrix "
            "with normal-approximation p-values."
        ),
        quantitative_meaning="substrate_set_enrichment_z_score",
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


def build_score_preconditioning_policy(
    *,
    policy: str,
    input_row_count: int,
    dropped_all_missing_row_count: int,
    retained_row_count: int,
    row_retention_rule: str = "drop_rows_with_all_scores_missing",
    retained_partial_missing_rows: bool = True,
) -> ScientificPolicyRecord:
    base_policy = resolve_score_preconditioning_policy(policy=policy)
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
        name=base_policy.name,
        version="1",
        description=(
            "Preconditions aligned downstream score rows before signalome "
            "construction by handling unsupported all-missing rows explicitly."
        ),
        parameters={
            "policy": str(policy),
            "row_retention_rule": str(row_retention_rule),
            "retained_partial_missing_rows": bool(retained_partial_missing_rows),
            "input_row_count": int(input_row_count),
            "dropped_all_missing_row_count": int(dropped_all_missing_row_count),
            "retained_row_count": int(retained_row_count),
        },
        assumptions=(
            "All-missing score rows are scientifically unsupported for score-driven "
            "signalome construction.",
            "Preconditioning policy determines whether row dropping is allowed or "
            "treated as a boundary error.",
            "Row retention changes site coverage and therefore can change final "
            "signalome assignments and module summaries.",
        ),
        output_scale="Retained downstream score matrix rows for signalome execution.",
        quantitative_meaning="retained_signalome_score_rows",
    )


SCORE_PRECONDITIONING_ERROR_ON_DROP_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
    name="score_preconditioning_error_on_drop_v1",
    version="1",
    description=(
        "Treat dropped all-missing downstream score rows as a hard workflow boundary "
        "error."
    ),
    parameters={
        "policy": "error_on_drop",
        "row_retention_rule": "drop_rows_with_all_scores_missing",
        "retained_partial_missing_rows": True,
    },
    assumptions=(
        "All interpreted sites must retain at least one finite downstream score.",
        "Any unsupported all-missing row invalidates signalome construction.",
    ),
    output_scale="Validation-only policy for preconditioning gate behavior.",
    quantitative_meaning="preconditioning_validation_rule",
)


SCORE_PRECONDITIONING_ALLOW_AND_REPORT_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
    name="score_preconditioning_allow_and_report_v1",
    version="1",
    description=(
        "Allow dropping all-missing downstream score rows and report the drop "
        "diagnostics in provenance."
    ),
    parameters={
        "policy": "allow_and_report",
        "row_retention_rule": "drop_rows_with_all_scores_missing",
        "retained_partial_missing_rows": True,
    },
    assumptions=(
        "All-missing rows provide no usable downstream score evidence.",
        "Retained-site coverage can change as unsupported rows are removed.",
    ),
    output_scale="Validation-and-row-retention policy for score preconditioning.",
    quantitative_meaning="preconditioning_row_retention_rule",
)


def resolve_score_preconditioning_policy(*, policy: str) -> ScientificPolicyRecord:
    if policy == "error_on_drop":
        return SCORE_PRECONDITIONING_ERROR_ON_DROP_POLICY
    if policy == "allow_and_report":
        return SCORE_PRECONDITIONING_ALLOW_AND_REPORT_POLICY
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING,
        name=f"score_preconditioning_{policy}_v1",
        version="1",
        description="Configured score preconditioning policy.",
        parameters={"policy": str(policy)},
        assumptions=(),
        output_scale="Validation-and-row-retention policy for score preconditioning.",
        quantitative_meaning="preconditioning_row_retention_rule",
    )


def build_preprocessing_stage_order_policy(
    *,
    configured_stage_order: tuple[str, ...],
    default_stage_order: tuple[str, ...],
    supported_stage_order: tuple[str, ...],
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.PREPROCESSING_STAGE_ORDER,
        name="Preprocessing Stage Order Policy",
        version="1",
        description=(
            "Defines the explicit stage execution order used to transform dataset "
            "inputs into analysis-ready workflow matrices."
        ),
        parameters={
            "configured_stage_order": " -> ".join(configured_stage_order),
            "configured_stage_count": int(len(configured_stage_order)),
            "default_stage_order": " -> ".join(default_stage_order),
            "default_stage_count": int(len(default_stage_order)),
            "supported_stage_order": " -> ".join(supported_stage_order),
            "supported_stage_count": int(len(supported_stage_order)),
        },
        assumptions=(
            "Stage order can change output row retention, transformed values, and "
            "derived comparison tables.",
            "Configured order must be interpreted as part of the scientific method "
            "for reproducibility.",
        ),
        output_scale="Ordered preprocessing execution plan for dataset construction.",
        quantitative_meaning="preprocessing_execution_order",
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


DUPLICATE_SITE_RESOLUTION_AGGREGATE_MEAN_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.DUPLICATE_SITE_RESOLUTION,
    name="duplicate_site_resolution_aggregate_mean_v1",
    version="1",
    description=(
        "Resolves duplicate constructed site IDs by aggregating duplicate rows "
        "with a column-wise arithmetic mean."
    ),
    parameters={
        "duplicate_site_policy": "aggregate_mean",
        "aggregation_function": "column_mean",
    },
    assumptions=(
        "All duplicate rows contribute numerically to one retained site row.",
        "Aggregation can blend peptide-context-specific measurements.",
    ),
    output_scale="Deduplicated site-matrix rows for downstream scoring.",
    quantitative_meaning="duplicate_site_resolved_matrix",
)


def build_duplicate_site_resolution_policy(
    *,
    duplicate_site_policy: str,
) -> ScientificPolicyRecord:
    if duplicate_site_policy == "aggregate_mean":
        return DUPLICATE_SITE_RESOLUTION_AGGREGATE_MEAN_POLICY
    return ScientificPolicyRecord(
        id=ScientificPolicyId.DUPLICATE_SITE_RESOLUTION,
        name=f"duplicate_site_resolution_{duplicate_site_policy}_v1",
        version="1",
        description=(
            "Resolves duplicate constructed site IDs according to the configured "
            "site-matrix duplicate-site policy."
        ),
        parameters={"duplicate_site_policy": str(duplicate_site_policy)},
        assumptions=(
            "Duplicate-site resolution policy changes retained rows and/or values.",
            "Resolved rows become the authoritative site matrix for downstream use.",
        ),
        output_scale="Deduplicated site-matrix rows for downstream scoring.",
        quantitative_meaning="duplicate_site_resolved_matrix",
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


def _threshold_operator_token(mode: ThresholdMode) -> str:
    if mode is ThresholdMode.GREATER_THAN:
        return ">"
    if mode is ThresholdMode.GREATER_THAN_OR_EQUAL:
        return ">="
    return ">"


__all__ = [
    "CandidateSubstrateSelectionPolicy",
    "KinaseProfileScoringPolicy",
    "PreprocessingStageOrderPolicy",
    "PROFILE_CORRELATION_SHIFTED_UNIT_POLICY",
    "PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY",
    "DUPLICATE_SITE_RESOLUTION_AGGREGATE_MEAN_POLICY",
    "ScorePreconditioningPolicy",
    "SCORE_PRECONDITIONING_ALLOW_AND_REPORT_POLICY",
    "SCORE_PRECONDITIONING_ERROR_ON_DROP_POLICY",
    "SignalomeMissingValueClusteringPolicy",
    "ScientificPolicyId",
    "ScientificPolicyRecord",
    "build_candidate_substrate_selection_policy",
    "build_duplicate_site_resolution_policy",
    "build_ksea_zscore_activity_policy",
    "build_kinase_profile_scoring_policy",
    "build_motif_profile_rank_fusion_policy",
    "build_peptide_to_site_aggregation_policy",
    "build_preprocessing_stage_order_policy",
    "build_score_preconditioning_policy",
    "build_signalome_missing_value_clustering_policy",
    "build_signalome_module_candidate_score_policy",
    "build_simplified_weighted_substrate_activity_policy",
    "resolve_score_preconditioning_policy",
]
