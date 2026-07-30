"""Kinase workflow provenance assembly."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence

import pandas as pd

from phospy.contracts.configs import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
    KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY,
    LocalisationRequirement,
)
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table_normalized_axes
from phospy.provenance.models import (
    EnvironmentProvenance,
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyRecord,
)
from phospy.science.activities.methods import (
    KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG,
    SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG,
    SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE,
    SSGSEA_SIGNIFICANCE_STATUS_P_VALUE_AVAILABLE_Q_VALUE_DISABLED,
    SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS,
)
from phospy.science.activities.models import KinaseActivityResult
from phospy.science.activities.scientific_policies import (
    SSGSEA_PERMUTATION_RNG_SEED_POLICY,
    SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION,
    build_ksea_zscore_activity_policy,
    build_simplified_weighted_substrate_activity_policy,
    build_ssgsea_substrate_enrichment_activity_policy,
)
from phospy.science.datasets.preprocessing.scientific_policies import (
    build_duplicate_site_resolution_policy,
)
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.prediction.scientific_policies import (
    PROFILE_CORRELATION_SHIFTED_UNIT_POLICY,
    CandidateSubstrateSelectionPolicy,
    KinaseProfileScoringPolicy,
    build_kinase_library_motif_scoring_policy,
    build_motif_profile_rank_fusion_policy,
)
from phospy.science.tables.kinase import (
    KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
)
from phospy.workflows.intensity_scale_evidence import (
    input_intensity_scale_evidence_payload,
)
from phospy.workflows.kinase.attrition_metrics import (
    build_kinase_attrition_provenance_payload,
    kinase_attrition_policy_to_payload,
)
from phospy.workflows.kinase.component_models import (
    CANDIDATE_MIN_INCLUSION,
    CANDIDATE_SCORE_THRESHOLD,
)
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.row_attrition import (
    build_kinase_row_attrition_provenance,
)


class KinaseProvenanceBuilder:
    """Build workflow-level provenance for kinase execution."""

    def __init__(
        self,
        *,
        collect_environment: Callable[
            [], EnvironmentProvenance
        ] = collect_environment_provenance,
    ) -> None:
        self._collect_environment = collect_environment

    def run(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult,
        activity_result: KinaseActivityResult | None,
        substrate_contributions: pd.DataFrame | None = None,
    ) -> RunProvenance:
        input_tables = _build_input_table_fingerprints(request)
        output_tables = _build_output_table_fingerprints(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
            substrate_contributions=substrate_contributions,
        )
        workflow_parameters = _build_workflow_parameters(
            request=request,
            config=config,
            scoring_result=scoring_result,
            activity_result=activity_result,
        )
        duplicate_site_policy = self._resolve_duplicate_site_resolution_policy(
            request=request
        )
        scientific_policies = _build_scientific_policy_records(
            config=config,
            scoring_result=scoring_result,
            duplicate_site_policy=duplicate_site_policy,
        )

        return RunProvenance(
            environment=self._collect_environment(),
            input_tables=input_tables,
            preprocessing_stages=self._dataset_preprocessing_stages(request),
            reference=request.references.provenance,
            workflow_name="kinase_workflow",
            workflow_parameters=workflow_parameters,
            random_state=config.prediction_random_state,
            random_seed_policy=self._resolve_seed_policy(config),
            output_tables=output_tables,
            scientific_policies=scientific_policies,
            reference_context=request.dataset.reference_context,
        )

    @staticmethod
    def _resolve_seed_policy(config: ResolvedKinaseExecutionConfig) -> str | None:
        if config.prediction_mode == KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING:
            return None
        if config.prediction_mode != KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE:
            return None
        return config.prediction_sampling_policy.seed_strategy

    @staticmethod
    def _dataset_preprocessing_stages(
        request: ResolvedKinaseWorkflowRequest,
    ) -> tuple[PreprocessingStageProvenance, ...]:
        provenance = request.dataset.provenance
        if provenance is None:
            return ()
        return tuple(provenance.preprocessing_stages)

    @staticmethod
    def _resolve_duplicate_site_resolution_policy(
        *,
        request: ResolvedKinaseWorkflowRequest,
    ) -> ScientificPolicyRecord | None:
        dataset_provenance = request.dataset.provenance
        if dataset_provenance is None:
            return None
        workflow_parameters = dataset_provenance.workflow_parameters
        if not isinstance(workflow_parameters, Mapping):
            return None
        preprocessing_plan = workflow_parameters.get("preprocessing_plan")
        if not isinstance(preprocessing_plan, Mapping):
            return None
        stage_order = preprocessing_plan.get("stage_order")
        if (
            not isinstance(stage_order, Sequence)
            or isinstance(stage_order, (str, bytes, bytearray))
            or "site_matrix" not in stage_order
        ):
            return None
        duplicate_site_policy = preprocessing_plan.get(
            "site_matrix_duplicate_site_policy"
        )
        if not isinstance(duplicate_site_policy, str) or not duplicate_site_policy:
            return None
        return build_duplicate_site_resolution_policy(
            duplicate_site_policy=duplicate_site_policy
        )


def _build_input_table_fingerprints(
    request: ResolvedKinaseWorkflowRequest,
) -> tuple[TableFingerprint, ...]:
    return _collect_fingerprints(
        (
            ("dataset.phospho", request.dataset.phospho),
            ("dataset.site_metadata", request.dataset.site_metadata),
            ("dataset.sample_metadata", request.dataset.sample_metadata),
            ("dataset.total", request.dataset.total),
            ("dataset.comparisons", request.dataset.comparisons),
            (
                "references.kinase_substrate_map",
                request.references.kinase_substrate_map,
            ),
            ("references.site_sequences", request.references.site_sequences),
        )
    )


def _build_output_table_fingerprints(
    *,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult,
    activity_result: KinaseActivityResult | None,
    substrate_contributions: pd.DataFrame | None = None,
) -> tuple[TableFingerprint, ...]:
    return _collect_fingerprints(
        (
            (
                "outputs.scoring.profile_scores",
                _profile_scores_for_output_fingerprint(scoring_result),
            ),
            (
                "outputs.scoring.profile_score_diagnostics",
                scoring_result.profile_score_diagnostics,
            ),
            ("outputs.scoring.motif_scores", scoring_result.motif_scores),
            (
                "outputs.scoring.rank_weighted_fusion_scores",
                scoring_result.rank_weighted_fusion_scores,
            ),
            (
                "outputs.scoring.kinase_library_motif_scores",
                scoring_result.kinase_library_motif_scores,
            ),
            (
                "outputs.scoring.combined_profile_motif_scores",
                scoring_result.combined_profile_motif_scores,
            ),
            (
                "outputs.scoring.score_fusion_weights",
                scoring_result.score_fusion_weights,
            ),
            (
                "outputs.scoring.kinase_library_site_diagnostics",
                scoring_result.kinase_library_site_diagnostics,
            ),
            (
                "outputs.scoring.kinase_library_kinase_diagnostics",
                scoring_result.kinase_library_kinase_diagnostics,
            ),
            ("outputs.prediction.pred_mat", prediction_result.pred_mat),
            ("outputs.prediction.substrate_list", prediction_result.substrate_list),
            ("outputs.substrate_contributions", substrate_contributions),
            (
                "outputs.activity.weighted_activity",
                None if activity_result is None else activity_result.activity_matrix,
            ),
            (
                "outputs.activity.thresholded_substrate_mean_activity",
                (
                    None
                    if activity_result is None
                    else activity_result.thresholded_substrate_mean_activity
                ),
            ),
            (
                "outputs.activity.thresholded_substrate_counts",
                (
                    None
                    if activity_result is None
                    else activity_result.thresholded_substrate_counts.to_frame(
                        name="n_substrates"
                    )
                ),
            ),
            (
                "outputs.activity.activity_substrate_counts",
                (
                    None
                    if activity_result is None
                    else activity_result.activity_substrate_counts
                ),
            ),
            (
                "outputs.activity.target_counts",
                (
                    None
                    if activity_result is None
                    else activity_result.target_counts.to_frame(name="n_targets")
                ),
            ),
            (
                "outputs.activity.target_table",
                None if activity_result is None else activity_result.target_table,
            ),
            (
                "outputs.activity.statistics_table",
                None if activity_result is None else activity_result.statistics_table,
            ),
        )
    )


def _profile_scores_for_output_fingerprint(
    scoring_result: KinaseScoringResult,
) -> pd.DataFrame | None:
    metadata = (
        {}
        if scoring_result.score_scale_metadata is None
        else dict(scoring_result.score_scale_metadata)
    )
    if metadata.get("uses_profile_correlation") is False:
        return None
    return scoring_result.profile_scores


def _collect_fingerprints(
    entries: tuple[tuple[str, pd.DataFrame | None], ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table_normalized_axes(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _build_workflow_parameters(
    *,
    request: ResolvedKinaseWorkflowRequest,
    config: ResolvedKinaseExecutionConfig,
    scoring_result: KinaseScoringResult,
    activity_result: KinaseActivityResult | None,
) -> dict[str, object]:
    payload = input_intensity_scale_evidence_payload(request.dataset)
    row_attrition = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=scoring_result,
    )
    payload.update(
        {
            "site_token_validation": _build_site_token_validation_payload(request),
            "scoring_config": _build_scoring_config_payload(
                request=request,
                scoring_result=scoring_result,
                config=config,
            ),
            "attrition_provenance": _build_attrition_provenance_payload(
                request=request,
                config=config,
            ),
            "scoring_diagnostics": _build_scoring_diagnostics_payload(
                request=request,
                scoring_result=scoring_result,
            ),
            "prediction_config": _build_prediction_config_payload(config),
            "activity_config": _build_activity_config_payload(
                config=config,
                activity_result=activity_result,
            ),
            **row_attrition.to_workflow_parameters(),
        }
    )
    return payload


def _build_site_token_validation_payload(
    request: ResolvedKinaseWorkflowRequest,
) -> dict[str, object]:
    return {
        "mode": (
            "opaque_opt_in"
            if request.dataset.opaque_site_values_allowed
            else "strict_sty_residue_position"
        )
    }


def _build_prediction_config_payload(
    config: ResolvedKinaseExecutionConfig,
) -> dict[str, object]:
    return {
        "top_k": int(config.prediction_top_k),
        "deterministic_max_selected_kinases": int(
            config.prediction_deterministic_max_selected_kinases
        ),
        "adaptive_ensemble_runs": int(config.prediction_adaptive_ensemble_runs),
        "mode": str(config.prediction_mode),
        "adaptive_policy": str(config.prediction_adaptive_policy),
        "n_iterations": int(config.prediction_n_iterations),
        "random_state": config.prediction_random_state,
    }


def _build_activity_config_payload(
    *,
    config: ResolvedKinaseExecutionConfig,
    activity_result: KinaseActivityResult | None,
) -> dict[str, object] | None:
    if config.activity is None:
        return None
    ssgsea_significance = _build_ssgsea_significance_payload(
        config=config,
        activity_result=activity_result,
    )
    return {
        "method": str(config.activity.method),
        "threshold": float(config.activity.threshold),
        "min_substrates": int(config.activity.min_substrates),
        "top_n_substrates": int(config.activity.top_n_substrates),
        "ksea_min_substrates": int(config.activity.ksea_min_substrates),
        "ksea_evidence_threshold": float(config.activity.ksea_evidence_threshold),
        "ksea_p_value_method": str(config.activity.ksea_p_value_method),
        "ksea_adjust_p_values": bool(config.activity.ksea_adjust_p_values),
        "ksea_formula_version": "v1",
        "ksea_q_value_method": (
            KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG
            if config.activity.ksea_adjust_p_values
            else None
        ),
        "ssgsea_min_substrates": int(config.activity.ssgsea_min_substrates),
        "ssgsea_ranking_direction": str(config.activity.ssgsea_ranking_direction),
        "ssgsea_permutations": int(config.activity.ssgsea_permutations),
        "ssgsea_random_seed": (
            None
            if config.activity.ssgsea_random_seed is None
            else int(config.activity.ssgsea_random_seed)
        ),
        "ssgsea_permutation_rng_seed_policy": (
            SSGSEA_PERMUTATION_RNG_SEED_POLICY
            if int(config.activity.ssgsea_permutations) > 0
            else None
        ),
        "ssgsea_permutation_rng_seed_policy_version": (
            SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION
            if int(config.activity.ssgsea_permutations) > 0
            else None
        ),
        "ssgsea_adjust_p_values": bool(config.activity.ssgsea_adjust_p_values),
        "ssgsea_q_value_method": (
            SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG
            if config.activity.ssgsea_permutations > 0
            and config.activity.ssgsea_adjust_p_values
            else None
        ),
        "ssgsea_significance_status": ssgsea_significance["status"],
        "ssgsea_significance_status_counts": ssgsea_significance["status_counts"],
        "activity_method": (
            None
            if activity_result is None
            else activity_result.activity_method.to_payload()
        ),
        "activity_method_summary": (
            None
            if activity_result is None or activity_result.method_summary is None
            else activity_result.method_summary.to_payload()
        ),
        "threshold_membership_diagnostics": (
            None
            if activity_result is None
            or activity_result.threshold_membership_diagnostics is None
            else activity_result.threshold_membership_diagnostics.to_payload()
        ),
    }


def _build_scoring_config_payload(
    *,
    request: ResolvedKinaseWorkflowRequest,
    scoring_result: KinaseScoringResult,
    config: ResolvedKinaseExecutionConfig,
) -> dict[str, object]:
    profile_correlation_enabled = (
        str(config.scoring_mode) != KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY
    )
    payload: dict[str, object] = {
        "requested_reliability_profile": (
            None
            if config.requested_reliability_profile is None
            else str(config.requested_reliability_profile)
        ),
        "effective_reliability_profile": str(config.effective_reliability_profile),
        "include_diagnostic_scoring_tables": bool(
            config.include_diagnostic_scoring_tables
        ),
        "profile_self_inclusion_policy": str(config.profile_self_inclusion_policy),
        "localisation_requirement": _localisation_requirement_to_payload(
            config.localisation_requirement
        ),
        "reference_context_compatibility_policy": str(
            config.reference_context_compatibility_policy
        ),
        "attrition_policy": kinase_attrition_policy_to_payload(config.attrition_policy),
    }
    if profile_correlation_enabled:
        payload.update(
            {
                "min_substrates": int(config.scoring_min_substrates),
                "profile_missing_value_strategy": str(
                    config.profile_missing_value_strategy
                ),
            }
        )
    else:
        payload.update(
            {
                "uses_profile_correlation": False,
                "uses_reference_substrate_profiles": False,
                "uses_sequence_motif_resource": True,
            }
        )
    if config.include_substrate_contributions:
        payload["include_substrate_contributions"] = True
    if config.scoring_mode == KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED:
        return payload
    payload.update(
        {
            "scoring_mode": str(config.scoring_mode),
            "score_source": scoring_result.score_source,
            "score_scale": scoring_result.score_scale,
            "kinase_library_resource": (
                None
                if request.kinase_library_resource is None
                else dict(request.kinase_library_resource.provenance.manifest or {})
            ),
        }
    )
    return payload


def _localisation_requirement_to_payload(
    requirement: LocalisationRequirement,
) -> dict[str, object]:
    return {
        "policy": str(requirement.policy),
        "require_present": bool(requirement.require_present),
        "minimum_probability": (
            None
            if requirement.minimum_probability is None
            else float(requirement.minimum_probability)
        ),
    }


def _build_ssgsea_significance_payload(
    *,
    config: ResolvedKinaseExecutionConfig,
    activity_result: KinaseActivityResult | None,
) -> dict[str, object]:
    activity_config = config.activity
    if activity_config is None:
        return {"status": None, "status_counts": None}
    if activity_config.method != KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT:
        return {"status": None, "status_counts": None}
    status: str
    if int(activity_config.ssgsea_permutations) <= 0:
        status = SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS
    elif not bool(activity_config.ssgsea_adjust_p_values):
        status = SSGSEA_SIGNIFICANCE_STATUS_P_VALUE_AVAILABLE_Q_VALUE_DISABLED
    else:
        status = SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE
    status_counts = _ssgsea_significance_status_counts(activity_result)
    return {"status": status, "status_counts": status_counts}


def _ssgsea_significance_status_counts(
    activity_result: KinaseActivityResult | None,
) -> dict[str, int] | None:
    if activity_result is None:
        return None
    statistics_table = activity_result.statistics_table
    if statistics_table is None or "significance_status" not in statistics_table:
        return None
    return {
        str(key): int(value)
        for key, value in statistics_table.loc[:, "significance_status"]
        .astype(str)
        .value_counts()
        .sort_index()
        .items()
    }


def _build_attrition_provenance_payload(
    *,
    request: ResolvedKinaseWorkflowRequest,
    config: ResolvedKinaseExecutionConfig,
) -> dict[str, object]:
    metrics = request.attrition_metrics
    if metrics is None:
        raise RuntimeError("kinase provenance requires resolved attrition metrics")
    return build_kinase_attrition_provenance_payload(
        metrics=metrics,
        policy=config.attrition_policy,
        violations=request.attrition_policy_violations,
    )


def _build_scoring_diagnostics_payload(
    *,
    request: ResolvedKinaseWorkflowRequest,
    scoring_result: KinaseScoringResult,
) -> dict[str, object]:
    scoring_diagnostics: dict[str, object] = (
        {}
        if scoring_result.motif_sequence_validation is None
        else dict(scoring_result.motif_sequence_validation.summary())
    )
    if scoring_result.motif_sequence_validation is not None:
        scoring_diagnostics["motif_site_sequence_coverage"] = (
            scoring_result.motif_sequence_validation.site_sequence_coverage_summary()
        )
    if scoring_result.motif_library_validation is not None:
        scoring_diagnostics["motif_library_validation"] = (
            scoring_result.motif_library_validation.summary()
        )
    if scoring_result.kinase_library_site_diagnostics is not None:
        scoring_diagnostics["kinase_library_site_status_counts"] = {
            str(key): int(value)
            for key, value in scoring_result.kinase_library_site_diagnostics.loc[
                :, "status"
            ]
            .astype(str)
            .value_counts()
            .sort_index()
            .items()
        }
    if scoring_result.kinase_library_kinase_diagnostics is not None:
        scoring_diagnostics["kinase_library_matrix_status_counts"] = {
            str(key): int(value)
            for key, value in scoring_result.kinase_library_kinase_diagnostics.loc[
                :, "status"
            ]
            .astype(str)
            .value_counts()
            .sort_index()
            .items()
        }
    if scoring_result.score_scale_metadata is not None:
        score_scale_metadata = dict(scoring_result.score_scale_metadata)
        scoring_diagnostics["score_scale_metadata"] = score_scale_metadata
        for evidence_key in (
            "uses_profile_correlation",
            "uses_reference_substrate_profiles",
            "uses_sequence_motif_resource",
        ):
            if evidence_key in score_scale_metadata:
                scoring_diagnostics[evidence_key] = bool(
                    score_scale_metadata[evidence_key]
                )
    if scoring_result.score_source_summary is not None:
        score_source_summary = scoring_result.score_source_summary
        by_kinase: dict[str, dict[str, int]] = {}
        for kinase, row in score_source_summary.iterrows():
            by_kinase[str(kinase)] = {
                str(column): int(value) for column, value in row.items()
            }
        scoring_diagnostics["kinase_score_source_counts_by_kinase"] = by_kinase
        scoring_diagnostics["kinase_score_source_counts_total"] = {
            str(column): int(value)
            for column, value in score_source_summary.sum(axis=0).items()
        }
    profile_score_diagnostics = scoring_result.profile_score_diagnostics
    if profile_score_diagnostics is not None:
        scoring_diagnostics["profile_score_diagnostics"] = (
            _profile_score_diagnostics_payload(profile_score_diagnostics)
        )
    if request.site_sequence_merge_diagnostics:
        scoring_diagnostics["site_sequence_merge"] = dict(
            request.site_sequence_merge_diagnostics
        )
    if request.attrition_metrics is not None:
        scoring_diagnostics["attrition_metrics"] = (
            request.attrition_metrics.to_payload()
        )
    if request.attrition_policy_violations:
        scoring_diagnostics["attrition_policy_violations"] = [
            violation.to_payload() for violation in request.attrition_policy_violations
        ]
    return scoring_diagnostics


def _profile_score_diagnostics_payload(
    diagnostics: pd.DataFrame,
) -> dict[str, object]:
    status_counts = {
        str(key): int(value)
        for key, value in diagnostics.loc[:, "status"]
        .astype(str)
        .value_counts()
        .items()
    }
    reason_series = diagnostics.loc[
        diagnostics.loc[:, "status"].astype(str)
        == KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
        "reason",
    ]
    reason_counts = {
        str(key): int(value)
        for key, value in reason_series.astype(str).value_counts().items()
    }
    return {
        "row_count": int(diagnostics.shape[0]),
        "status_counts": status_counts,
        "unscored_reason_counts": reason_counts,
    }


def _build_scientific_policy_records(
    *,
    config: ResolvedKinaseExecutionConfig,
    scoring_result: KinaseScoringResult,
    duplicate_site_policy: ScientificPolicyRecord | None,
) -> tuple[ScientificPolicyRecord, ...]:
    scientific_policies: list[ScientificPolicyRecord] = []
    if str(config.scoring_mode) != KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY:
        scientific_policies.extend(
            [
                PROFILE_CORRELATION_SHIFTED_UNIT_POLICY,
                KinaseProfileScoringPolicy(
                    profile_missing_value_strategy=str(
                        config.profile_missing_value_strategy
                    ),
                    min_substrates_floor=KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
                    requested_min_substrates=int(config.scoring_min_substrates),
                    profile_self_inclusion_policy=(
                        config.profile_self_inclusion_policy
                    ),
                ).record,
            ]
        )
    if config.scoring_mode in {
        KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    }:
        scientific_policies.append(
            build_motif_profile_rank_fusion_policy(
                allow_profile_only_fallback=True,
                emit_weights=bool(config.include_diagnostic_scoring_tables),
            )
        )
    kinase_library_policy = _build_kinase_library_scoring_policy(
        config=config,
        scoring_result=scoring_result,
    )
    if kinase_library_policy is not None:
        scientific_policies.append(kinase_library_policy)
    scientific_policies.append(
        CandidateSubstrateSelectionPolicy(
            top_k=int(config.prediction_top_k),
            score_threshold=CANDIDATE_SCORE_THRESHOLD,
            inclusion=CANDIDATE_MIN_INCLUSION,
        ).record
    )
    if config.prediction_mode == KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE:
        scientific_policies.append(config.prediction_sampling_policy.record)
    if duplicate_site_policy is not None:
        scientific_policies.append(duplicate_site_policy)
    activity_policy = _build_activity_policy_record(config)
    if activity_policy is not None:
        scientific_policies.append(activity_policy)
    return tuple(scientific_policies)


def _build_kinase_library_scoring_policy(
    *,
    config: ResolvedKinaseExecutionConfig,
    scoring_result: KinaseScoringResult,
) -> ScientificPolicyRecord | None:
    if config.scoring_mode not in KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY:
        return None
    score_scale_metadata = (
        {}
        if scoring_result.score_scale_metadata is None
        else dict(scoring_result.score_scale_metadata)
    )
    raw_sequence_window = score_scale_metadata.get("sequence_window")
    sequence_window: Mapping[str, object] | None = (
        {str(key): value for key, value in raw_sequence_window.items()}
        if isinstance(raw_sequence_window, Mapping)
        else None
    )
    return build_kinase_library_motif_scoring_policy(
        scoring_mode=str(config.scoring_mode),
        resource_source_name=_optional_metadata_text(
            score_scale_metadata.get("resource_source_name")
        ),
        resource_source_version=_optional_metadata_text(
            score_scale_metadata.get("resource_source_version")
        ),
        resource_score_scale=_optional_metadata_text(
            score_scale_metadata.get("resource_score_scale")
        ),
        workflow_score_scale=str(scoring_result.score_scale),
        sequence_window=sequence_window,
    )


def _build_activity_policy_record(
    config: ResolvedKinaseExecutionConfig,
) -> ScientificPolicyRecord | None:
    if config.activity is None:
        return None
    if (
        config.activity.method
        == KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
    ):
        return build_simplified_weighted_substrate_activity_policy(
            threshold=float(config.activity.threshold),
            min_substrates=int(config.activity.min_substrates),
            top_n_substrates=int(config.activity.top_n_substrates),
        )
    if config.activity.method == KINASE_ACTIVITY_METHOD_KSEA_ZSCORE:
        return build_ksea_zscore_activity_policy(
            evidence_threshold=float(config.activity.ksea_evidence_threshold),
            min_substrates=int(config.activity.ksea_min_substrates),
            p_value_method=str(config.activity.ksea_p_value_method),
            adjust_p_values=bool(config.activity.ksea_adjust_p_values),
            q_value_method=(
                KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG
                if config.activity.ksea_adjust_p_values
                else None
            ),
        )
    if config.activity.method == KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT:
        return build_ssgsea_substrate_enrichment_activity_policy(
            min_substrates=int(config.activity.ssgsea_min_substrates),
            ranking_direction=str(config.activity.ssgsea_ranking_direction),
            permutation_count=int(config.activity.ssgsea_permutations),
            random_seed=config.activity.ssgsea_random_seed,
            adjust_p_values=bool(config.activity.ssgsea_adjust_p_values),
            q_value_method=(
                SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG
                if config.activity.ssgsea_permutations > 0
                and config.activity.ssgsea_adjust_p_values
                else None
            ),
        )
    return None


def _optional_metadata_text(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


__all__ = ["KinaseProvenanceBuilder"]
