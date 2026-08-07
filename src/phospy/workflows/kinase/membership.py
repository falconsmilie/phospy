"""Workflow construction of activity membership-selection provenance."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.provenance.hashing import fingerprint_optional_table_normalized_axes
from phospy.provenance.models import TableFingerprint
from phospy.science.activities.membership import (
    ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
    ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF,
    ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED,
    ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF,
    KSEA_SELECTION_QUANTITATIVE_MATRIX_FINGERPRINT_NAME,
    KSEA_TESTED_QUANTITATIVE_MATRIX_FINGERPRINT_NAME,
    ActivityMembershipSelection,
    fingerprint_ksea_selection_quantitative_matrix,
    fingerprint_ksea_tested_quantitative_matrix,
    selected_substrate_universe_from_prediction_matrix,
)
from phospy.science.prediction.models import KinasePredictionResult
from phospy.science.scoring.policy_models import DownstreamScoreSource
from phospy.workflows.kinase.component_models import (
    CANDIDATE_MIN_INCLUSION,
    CANDIDATE_SCORE_THRESHOLD,
    KinaseScoringRunResult,
)
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)


def build_ksea_membership_selection(
    *,
    request: ResolvedKinaseWorkflowRequest,
    config: ResolvedKinaseExecutionConfig,
    scoring_execution: KinaseScoringRunResult,
    prediction_result: KinasePredictionResult,
    evidence_threshold: float,
    membership_matrix: pd.DataFrame | None = None,
) -> ActivityMembershipSelection:
    """Build KSEA membership provenance from resolved workflow stage outputs."""

    downstream_score_source = scoring_execution.downstream_score_source
    consumed_tested_matrix = _score_source_consumed_tested_matrix(
        downstream_score_source
    )
    source_category = _source_category(downstream_score_source)
    selection_quantitative_fingerprint = (
        _selection_quantitative_matrix_fingerprint(request)
        if consumed_tested_matrix
        else None
    )
    tested_quantitative_fingerprint = _tested_quantitative_matrix_fingerprint(request)
    pred_mat = (
        prediction_result.pred_mat if membership_matrix is None else membership_matrix
    )
    selected_kinases = tuple(str(value) for value in pred_mat.columns.tolist())
    selected_substrates = selected_substrate_universe_from_prediction_matrix(
        pred_mat,
        threshold=float(evidence_threshold),
    )
    threshold_top_k_policy = {
        "prediction_mode": str(config.prediction_mode),
        "prediction_top_k": int(config.prediction_top_k),
        "prediction_deterministic_max_selected_kinases": int(
            config.prediction_deterministic_max_selected_kinases
        ),
        "prediction_score_threshold": float(CANDIDATE_SCORE_THRESHOLD),
        "prediction_candidate_min_inclusion": int(CANDIDATE_MIN_INCLUSION),
        "ksea_evidence_threshold": float(evidence_threshold),
        "ksea_evidence_threshold_operator": ">=",
        "selection_quantitative_input": (
            KSEA_SELECTION_QUANTITATIVE_MATRIX_FINGERPRINT_NAME
            if consumed_tested_matrix
            else None
        ),
        "tested_quantitative_input": KSEA_TESTED_QUANTITATIVE_MATRIX_FINGERPRINT_NAME,
    }
    source_reference_fingerprints = _source_reference_fingerprints(request)
    if source_category == ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF:
        return ActivityMembershipSelection.sequence_only_motif(
            provider_method_identifier="prediction_matrix_thresholded_membership",
            provider_method_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
            threshold_top_k_policy=threshold_top_k_policy,
            source_reference_fingerprints=source_reference_fingerprints,
            tested_quantitative_matrix_fingerprint=tested_quantitative_fingerprint,
            selected_kinase_universe=selected_kinases,
            selected_substrate_universe=selected_substrates,
        )
    if source_category == ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF:
        return ActivityMembershipSelection.fused_profile_motif(
            selection_method="prediction_matrix_thresholded_membership",
            score_source=downstream_score_source.value,
            threshold_top_k_policy=threshold_top_k_policy,
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=(
                selection_quantitative_fingerprint
            ),
            tested_quantitative_matrix_fingerprint=tested_quantitative_fingerprint,
            consumed_tested_matrix=consumed_tested_matrix,
            selected_kinase_universe=selected_kinases,
            selected_substrate_universe=selected_substrates,
        )
    return ActivityMembershipSelection.profile_derived(
        selection_method="prediction_matrix_thresholded_membership",
        score_source=downstream_score_source.value,
        threshold_top_k_policy=threshold_top_k_policy,
        source_reference_fingerprints=source_reference_fingerprints,
        selection_quantitative_matrix_fingerprint=selection_quantitative_fingerprint,
        tested_quantitative_matrix_fingerprint=tested_quantitative_fingerprint,
        consumed_tested_matrix=consumed_tested_matrix,
        selected_kinase_universe=selected_kinases,
        selected_substrate_universe=selected_substrates,
    )


def _score_source_consumed_tested_matrix(
    score_source: DownstreamScoreSource,
) -> bool:
    if score_source is DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES:
        return False
    return True


def _source_category(score_source: DownstreamScoreSource) -> str:
    if score_source is DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES:
        return ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF
    if score_source is DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES:
        return ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF
    return ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED


def _source_reference_fingerprints(
    request: ResolvedKinaseWorkflowRequest,
) -> tuple[TableFingerprint, ...]:
    kinase_library_resource = request.kinase_library_resource
    if kinase_library_resource is not None:
        fingerprints = tuple(kinase_library_resource.provenance.table_fingerprints)
        if fingerprints:
            return fingerprints
    reference_provenance = request.references.provenance
    if reference_provenance is not None:
        fingerprints = tuple(reference_provenance.table_fingerprints)
        if fingerprints:
            return fingerprints
    return tuple(
        fingerprint
        for fingerprint in (
            fingerprint_optional_table_normalized_axes(
                request.references.kinase_substrate_map,
                name="references.kinase_substrate_map",
            ),
            fingerprint_optional_table_normalized_axes(
                request.references.site_sequences,
                name="references.site_sequences",
            ),
        )
        if fingerprint is not None
    )


def _selection_quantitative_matrix_fingerprint(
    request: ResolvedKinaseWorkflowRequest,
) -> TableFingerprint:
    return fingerprint_ksea_selection_quantitative_matrix(
        request.scoring_phospho_matrix
    )


def _tested_quantitative_matrix_fingerprint(
    request: ResolvedKinaseWorkflowRequest,
) -> TableFingerprint:
    return fingerprint_ksea_tested_quantitative_matrix(
        request.ksea_background_phospho_matrix
    )


def membership_selection_payload(
    selection: ActivityMembershipSelection | None,
) -> Mapping[str, object] | None:
    """Return a JSON-safe payload for optional workflow provenance."""

    return None if selection is None else selection.to_payload()


__all__ = ["build_ksea_membership_selection", "membership_selection_payload"]
