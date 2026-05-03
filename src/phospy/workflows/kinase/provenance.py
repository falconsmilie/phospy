"""Kinase workflow provenance assembly."""

from __future__ import annotations

from collections.abc import Callable, Mapping

import pandas as pd

from phospy.activities.methods import KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG
from phospy.activities.models import KinaseActivityResult
from phospy.api.configs import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
)
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    EnvironmentProvenance,
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.scientific_policies import (
    PROFILE_CORRELATION_SHIFTED_UNIT_POLICY,
    CandidateSubstrateSelectionPolicy,
    KinaseProfileScoringPolicy,
    ScientificPolicyRecord,
    build_duplicate_site_resolution_policy,
    build_ksea_zscore_activity_policy,
    build_motif_profile_rank_fusion_policy,
    build_simplified_weighted_substrate_activity_policy,
)
from phospy.workflows.kinase.component_models import (
    CANDIDATE_MIN_INCLUSION,
    CANDIDATE_SCORE_THRESHOLD,
)
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
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
    ) -> RunProvenance:
        input_tables = self._collect_fingerprints(
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
        output_tables = self._collect_fingerprints(
            (
                ("outputs.scoring.profile_scores", scoring_result.profile_scores),
                ("outputs.scoring.motif_scores", scoring_result.motif_scores),
                (
                    "outputs.scoring.rank_weighted_fusion_scores",
                    scoring_result.rank_weighted_fusion_scores,
                ),
                (
                    "outputs.scoring.score_fusion_weights",
                    scoring_result.score_fusion_weights,
                ),
                ("outputs.prediction.pred_mat", prediction_result.pred_mat),
                ("outputs.prediction.substrate_list", prediction_result.substrate_list),
                (
                    "outputs.activity.weighted_activity",
                    None
                    if activity_result is None
                    else activity_result.activity_scores,
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
                    None
                    if activity_result is None
                    else activity_result.statistics_table,
                ),
            )
        )
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
        scientific_policies = [
            PROFILE_CORRELATION_SHIFTED_UNIT_POLICY,
            KinaseProfileScoringPolicy(
                profile_missing_value_strategy=str(
                    config.profile_missing_value_strategy
                ),
                min_substrates_floor=KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
                requested_min_substrates=int(config.scoring_min_substrates),
            ).record,
            build_motif_profile_rank_fusion_policy(
                allow_profile_only_fallback=True,
                emit_weights=bool(config.include_diagnostic_scoring_tables),
            ),
            CandidateSubstrateSelectionPolicy(
                top_k=int(config.prediction_top_k),
                score_threshold=CANDIDATE_SCORE_THRESHOLD,
                inclusion=CANDIDATE_MIN_INCLUSION,
            ).record,
        ]
        if config.prediction_mode == KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE:
            scientific_policies.append(config.prediction_sampling_policy.record)
        duplicate_site_policy = self._resolve_duplicate_site_resolution_policy(
            request=request
        )
        if duplicate_site_policy is not None:
            scientific_policies.append(duplicate_site_policy)
        if config.activity is not None:
            if (
                config.activity.method
                == KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
            ):
                scientific_policies.append(
                    build_simplified_weighted_substrate_activity_policy(
                        threshold=float(config.activity.threshold),
                        min_substrates=int(config.activity.min_substrates),
                        top_n_substrates=int(config.activity.top_n_substrates),
                    )
                )
            elif config.activity.method == KINASE_ACTIVITY_METHOD_KSEA_ZSCORE:
                scientific_policies.append(
                    build_ksea_zscore_activity_policy(
                        evidence_threshold=float(
                            config.activity.ksea_evidence_threshold
                        ),
                        min_substrates=int(config.activity.ksea_min_substrates),
                        p_value_method=str(config.activity.ksea_p_value_method),
                        adjust_p_values=bool(config.activity.ksea_adjust_p_values),
                        q_value_method=(
                            KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG
                            if config.activity.ksea_adjust_p_values
                            else None
                        ),
                    )
                )

        return RunProvenance(
            environment=self._collect_environment(),
            input_tables=input_tables,
            preprocessing_stages=self._dataset_preprocessing_stages(request),
            reference=request.references.provenance,
            workflow_name="kinase_workflow",
            workflow_parameters={
                "scoring_config": {
                    "min_substrates": int(config.scoring_min_substrates),
                    "include_diagnostic_scoring_tables": bool(
                        config.include_diagnostic_scoring_tables
                    ),
                    "profile_missing_value_strategy": str(
                        config.profile_missing_value_strategy
                    ),
                },
                "scoring_diagnostics": scoring_diagnostics,
                "prediction_config": {
                    "top_k": int(config.prediction_top_k),
                    "deterministic_max_selected_kinases": int(
                        config.prediction_deterministic_max_selected_kinases
                    ),
                    "adaptive_ensemble_runs": int(
                        config.prediction_adaptive_ensemble_runs
                    ),
                    "mode": str(config.prediction_mode),
                    "adaptive_policy": str(config.prediction_adaptive_policy),
                    "n_iterations": int(config.prediction_n_iterations),
                    "random_state": config.prediction_random_state,
                },
                "activity_config": (
                    None
                    if config.activity is None
                    else {
                        "method": str(config.activity.method),
                        "threshold": float(config.activity.threshold),
                        "min_substrates": int(config.activity.min_substrates),
                        "top_n_substrates": int(config.activity.top_n_substrates),
                        "ksea_min_substrates": int(config.activity.ksea_min_substrates),
                        "ksea_evidence_threshold": float(
                            config.activity.ksea_evidence_threshold
                        ),
                        "ksea_p_value_method": str(config.activity.ksea_p_value_method),
                        "ksea_adjust_p_values": bool(
                            config.activity.ksea_adjust_p_values
                        ),
                        "ksea_formula_version": "v1",
                        "ksea_q_value_method": (
                            KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG
                            if config.activity.ksea_adjust_p_values
                            else None
                        ),
                        "activity_method": (
                            None
                            if activity_result is None
                            else activity_result.activity_method.to_payload()
                        ),
                        "activity_method_summary": (
                            None
                            if activity_result is None
                            or activity_result.method_summary is None
                            else activity_result.method_summary.to_payload()
                        ),
                    }
                ),
            },
            random_state=config.prediction_random_state,
            random_seed_policy=self._resolve_seed_policy(config),
            output_tables=output_tables,
            scientific_policies=tuple(scientific_policies),
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
        if not isinstance(stage_order, list) or "site_matrix" not in stage_order:
            return None
        duplicate_site_policy = preprocessing_plan.get(
            "site_matrix_duplicate_site_policy"
        )
        if not isinstance(duplicate_site_policy, str) or not duplicate_site_policy:
            return None
        return build_duplicate_site_resolution_policy(
            duplicate_site_policy=duplicate_site_policy
        )

    @staticmethod
    def _collect_fingerprints(
        entries: tuple[tuple[str, pd.DataFrame | None], ...],
    ) -> tuple[TableFingerprint, ...]:
        fingerprints: list[TableFingerprint] = []
        for name, table in entries:
            fingerprint = fingerprint_optional_table(table, name=name)
            if fingerprint is None:
                continue
            fingerprints.append(fingerprint)
        return tuple(fingerprints)


__all__ = ["KinaseProvenanceBuilder"]
