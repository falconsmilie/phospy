"""Structured caveats for kinase workflow results."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs import (
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    ReferenceContextCompatibilityPolicy,
)
from phospy.contracts.result_caveats import ResultCaveat
from phospy.science.prediction.models import KinaseScoringResult
from phospy.science.prediction.scoring import (
    KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT,
    KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP,
)
from phospy.science.scoring.policy_models import ProfileSelfInclusionPolicy
from phospy.science.tables.kinase import (
    KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
)
from phospy.validation.identity_contracts import (
    validate_reference_context_compatibility,
)
from phospy.workflows.intensity_scale_evidence import (
    build_declared_input_intensity_scale_caveat,
)
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.workflows.result_caveat_helpers import (
    build_direct_trusted_dataset_construction_caveat,
    build_localisation_policy_details,
    build_reference_context_compatibility_caveat,
    deduplicate_caveats,
    is_permissive_localisation_requirement,
)

KINASE_ATTRITION_POLICY_CAVEAT_CODE = "kinase_attrition_policy_violation"
KINASE_PERMISSIVE_LOCALISATION_POLICY_CAVEAT_CODE = (
    "kinase_permissive_localisation_policy"
)
KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE = "kinase_non_default_reference_source"
KINASE_REFERENCE_AUTO_RESOLUTION_CAVEAT_CODE = "kinase_reference_auto_resolution"
KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE = "kinase_reference_score_fallback"
KINASE_SCORING_LIMITATION_CAVEAT_CODE = "kinase_non_phosr_equivalent_scoring"
KINASE_LIBRARY_MOTIF_ONLY_CAVEAT_CODE = "kinase_library_motif_only_sequence_evidence"
KINASE_PROFILE_SELF_INCLUSION_CAVEAT_CODE = "kinase_profile_self_inclusion_allowed"
KINASE_PROFILE_LEAVE_ONE_OUT_CAVEAT_CODE = "kinase_profile_leave_one_out"
KINASE_DIRECT_TRUSTED_DATASET_CAVEAT_CODE = "kinase_direct_trusted_dataset_construction"


def build_kinase_result_caveats(
    *,
    request: ResolvedKinaseWorkflowRequest,
    scoring_result: KinaseScoringResult,
) -> tuple[ResultCaveat, ...]:
    """Build compact machine-readable caveats for kinase workflow results."""

    caveats: list[ResultCaveat] = []
    declared_input_scale = build_declared_input_intensity_scale_caveat(
        dataset=request.dataset,
        workflow_scope="kinase_scoring",
    )
    if declared_input_scale is not None:
        caveats.append(declared_input_scale)
    direct_construction = build_direct_trusted_dataset_construction_caveat(
        dataset=request.dataset,
        code=KINASE_DIRECT_TRUSTED_DATASET_CAVEAT_CODE,
        workflow_scope="kinase_scoring",
        workflow_label="kinase workflow",
    )
    if direct_construction is not None:
        caveats.append(direct_construction)
    caveats.extend(_attrition_policy_caveats(request))
    localisation = _permissive_localisation_caveat(request)
    if localisation is not None:
        caveats.append(localisation)
    reference_source = _non_default_reference_source_caveat(request)
    if reference_source is not None:
        caveats.append(reference_source)
    reference_context = _reference_context_unknown_caveat(request)
    if reference_context is not None:
        caveats.append(reference_context)
    reference_resolution = _reference_auto_resolution_caveat(request)
    if reference_resolution is not None:
        caveats.append(reference_resolution)
    score_fallback = _reference_score_fallback_caveat(scoring_result)
    if score_fallback is not None:
        caveats.append(score_fallback)
    motif_only = _kinase_library_motif_only_caveat(request, scoring_result)
    if motif_only is not None:
        caveats.append(motif_only)
    self_inclusion = _profile_self_inclusion_caveat(request, scoring_result)
    if self_inclusion is not None:
        caveats.append(self_inclusion)
    caveats.append(_scoring_limitation_caveat(request, scoring_result))
    return deduplicate_caveats(caveats)


def _attrition_policy_caveats(
    request: ResolvedKinaseWorkflowRequest,
) -> tuple[ResultCaveat, ...]:
    caveats: list[ResultCaveat] = []
    for violation in request.attrition_policy_violations:
        payload = violation.to_payload()
        code = payload.get("code", KINASE_ATTRITION_POLICY_CAVEAT_CODE)
        caveats.append(
            ResultCaveat(
                code=str(code),
                severity="warning",
                message=violation.message,
                details=payload,
            )
        )
    return tuple(caveats)


def _permissive_localisation_caveat(
    request: ResolvedKinaseWorkflowRequest,
) -> ResultCaveat | None:
    requirement = request.execution_config.localisation_requirement
    if not is_permissive_localisation_requirement(requirement):
        return None
    details = build_localisation_policy_details(
        site_metadata=request.dataset.site_metadata,
        requirement=requirement,
        workflow_scope="kinase_scoring",
    )
    policy = str(requirement.policy)
    return ResultCaveat(
        code=KINASE_PERMISSIVE_LOCALISATION_POLICY_CAVEAT_CODE,
        severity="warning",
        message=(
            "Kinase workflow localisation policy does not enforce a minimum "
            "localisation probability; unknown or low-confidence phosphosite "
            "localisation can remain in scoring inputs."
        ),
        details=details | {"policy_is_permissive": True, "policy": policy},
    )


def _non_default_reference_source_caveat(
    request: ResolvedKinaseWorkflowRequest,
) -> ResultCaveat | None:
    provenance = request.references.provenance
    source_type = None if provenance is None else provenance.source_type
    if source_type == "bundled":
        return None
    details = _reference_details(request)
    details["expected_default_source_type"] = "bundled"
    return ResultCaveat(
        code=KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE,
        severity="warning",
        message=(
            "Kinase workflow used a non-default reference source; scoring and "
            "prediction depend on the supplied reference bundle content."
        ),
        details=details,
    )


def _reference_context_unknown_caveat(
    request: ResolvedKinaseWorkflowRequest,
) -> ResultCaveat | None:
    policy = request.execution_config.reference_context_compatibility_policy
    allow_unknown = (
        policy is ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
    )
    if not allow_unknown:
        return None
    warning = validate_reference_context_compatibility(
        request.dataset.reference_context,
        None
        if request.references.provenance is None
        else request.references.provenance.reference_context,
        operation="kinase workflow result dataset/reference bundle",
        allow_unknown=allow_unknown,
    )
    if warning is None:
        return None
    return build_reference_context_compatibility_caveat(
        warning,
        policy=policy,
        workflow_scope="kinase_scoring",
    )


def _reference_auto_resolution_caveat(
    request: ResolvedKinaseWorkflowRequest,
) -> ResultCaveat | None:
    resolution = request.reference_resolution_details
    if str(resolution.get("reference_input_value", "")) != "auto":
        return None
    details = _reference_details(request)
    details.update({str(key): value for key, value in resolution.items()})
    return ResultCaveat(
        code=KINASE_REFERENCE_AUTO_RESOLUTION_CAVEAT_CODE,
        severity="info",
        message=(
            "Kinase workflow resolved ReferencePreset.AUTO from dataset.organism; "
            "results depend on the bundled reference lane selected at runtime."
        ),
        details=details,
    )


def _reference_score_fallback_caveat(
    scoring_result: KinaseScoringResult,
) -> ResultCaveat | None:
    summary = scoring_result.score_source_summary
    if summary is None or summary.empty:
        return None
    details = _score_source_summary_details(summary)
    fallback_count = _sum_int(
        summary, "profile_only_motif_missing_or_constant_count"
    ) + _sum_int(summary, "profile_only_no_motif_overlap_count")
    if fallback_count <= 0:
        return None
    return ResultCaveat(
        code=KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE,
        severity="warning",
        message=(
            "Kinase rank-weighted scoring used profile-only fallback for one or "
            "more kinase-site scores because motif/reference evidence was "
            "unavailable or unusable."
        ),
        details=details,
    )


def _kinase_library_motif_only_caveat(
    request: ResolvedKinaseWorkflowRequest,
    scoring_result: KinaseScoringResult,
) -> ResultCaveat | None:
    scoring_mode = str(request.execution_config.scoring_mode)
    if scoring_mode != KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY:
        return None
    authoritative = scoring_result.authoritative_scores
    return ResultCaveat(
        code=KINASE_LIBRARY_MOTIF_ONLY_CAVEAT_CODE,
        severity="warning",
        message=(
            "Kinase Library motif-only scoring uses sequence motif evidence only. "
            "It does not use quantified known-substrate profile correlation and "
            "should not be interpreted as kinase activity or causal kinase "
            "assignment."
        ),
        details={
            "scoring_mode": scoring_mode,
            "score_source": scoring_result.score_source,
            "score_scale": scoring_result.score_scale,
            "site_count": int(authoritative.shape[0]),
            "kinase_count": int(authoritative.shape[1]),
            "uses_profile_correlation": False,
            "uses_reference_substrate_profiles": False,
            "uses_sequence_motif_resource": True,
        },
    )


def _profile_self_inclusion_caveat(
    request: ResolvedKinaseWorkflowRequest,
    scoring_result: KinaseScoringResult,
) -> ResultCaveat | None:
    policy = request.execution_config.profile_self_inclusion_policy
    authoritative = scoring_result.authoritative_scores
    if policy is ProfileSelfInclusionPolicy.ALLOW:
        return ResultCaveat(
            code=KINASE_PROFILE_SELF_INCLUSION_CAVEAT_CODE,
            severity="warning",
            message=(
                "Kinase profile scoring allowed self-inclusion: a known substrate site "
                "may contribute to the kinase profile used to score that same site. "
                "Scores are exploratory and may be inflated for known substrates."
            ),
            details={
                "profile_self_inclusion_policy": policy.value,
                "self_inclusion_allowed": True,
                "leave_one_out_enabled": False,
                "site_count": int(authoritative.shape[0]),
                "kinase_count": int(authoritative.shape[1]),
            },
        )
    if policy is ProfileSelfInclusionPolicy.LEAVE_ONE_OUT:
        details = {
            "profile_self_inclusion_policy": policy.value,
            "self_inclusion_allowed": False,
            "leave_one_out_enabled": True,
            "site_count": int(authoritative.shape[0]),
            "kinase_count": int(authoritative.shape[1]),
        }
        details.update(_leave_one_out_diagnostic_details(scoring_result))
        return ResultCaveat(
            code=KINASE_PROFILE_LEAVE_ONE_OUT_CAVEAT_CODE,
            severity="info",
            message=(
                "Leave-one-out profile scoring was used. Known substrate sites "
                "were excluded from their own kinase profiles before scoring "
                "where applicable."
            ),
            details=details,
        )
    return None


def _scoring_limitation_caveat(
    request: ResolvedKinaseWorkflowRequest,
    scoring_result: KinaseScoringResult,
) -> ResultCaveat:
    scoring_mode = str(request.execution_config.scoring_mode)
    authoritative = scoring_result.authoritative_scores
    details = {
        "scoring_mode": scoring_mode,
        "score_source": scoring_result.score_source,
        "score_scale": scoring_result.score_scale,
        "default_scoring_mode": KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        "score_interpretation": "relative_support_within_run",
        "not_calibrated_probability": True,
        "not_causal_activity_proof": True,
        "site_count": int(authoritative.shape[0]),
        "kinase_count": int(authoritative.shape[1]),
        "kinase_library_resource_provided": (
            request.kinase_library_resource is not None
        ),
    }
    return ResultCaveat(
        code=KINASE_SCORING_LIMITATION_CAVEAT_CODE,
        severity="info",
        message=(
            "Kinase scoring outputs are relative support values under PhosPy "
            "workflow scoring semantics; they are not exact PhosR numerical "
            "equivalents and are not causal kinase activity proof."
        ),
        details=details,
    )


def _reference_details(request: ResolvedKinaseWorkflowRequest) -> dict[str, object]:
    provenance = request.references.provenance
    details: dict[str, object] = {
        "reference_organism": request.references.organism.value,
        "kinase_substrate_record_count": int(
            request.references.kinase_substrate_map.shape[0]
        ),
        "site_sequence_record_count": int(request.references.site_sequences.shape[0]),
    }
    if provenance is None:
        details["source_type"] = "unknown"
        return details
    details.update(
        {
            "source_type": provenance.source_type,
            "bundle_id": provenance.bundle_id,
            "source_name": provenance.source_name,
            "source_version": provenance.source_version,
            "identifier_namespace": provenance.identifier_namespace,
        }
    )
    return details


def _score_source_summary_details(summary: pd.DataFrame) -> dict[str, object]:
    fallback_columns = (
        "profile_only_motif_missing_or_constant_count",
        "profile_only_no_motif_overlap_count",
    )
    fused_column = "fused_motif_profile_evidence_count"
    unavailable_column = "unavailable_no_score_count"
    affected = summary.loc[:, list(fallback_columns)].sum(axis=1) > 0
    details: dict[str, object] = {
        "affected_kinase_count": int(affected.sum()),
        "kinase_count": int(summary.shape[0]),
        "score_source_summary_columns": [str(column) for column in summary],
        "profile_only_source_labels": [
            KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT,
            KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP,
        ],
    }
    for column in (*fallback_columns, fused_column, unavailable_column):
        details[column] = _sum_int(summary, column)
    details["total_score_cells"] = _sum_int(summary, "total_sites_count")
    return details


def _leave_one_out_diagnostic_details(
    scoring_result: KinaseScoringResult,
) -> dict[str, object]:
    diagnostics = scoring_result.profile_score_diagnostics
    if diagnostics is None:
        return {
            "profile_score_diagnostics_present": False,
            "leave_one_out_diagnostic_row_count": 0,
            "leave_one_out_unscored_cell_count": 0,
            "leave_one_out_unscored_reasons": {},
        }
    unscored = diagnostics.loc[
        diagnostics.loc[:, "status"].astype(str)
        == KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
        :,
    ]
    return {
        "profile_score_diagnostics_present": True,
        "leave_one_out_diagnostic_row_count": int(diagnostics.shape[0]),
        "leave_one_out_unscored_cell_count": int(unscored.shape[0]),
        "leave_one_out_unscored_reasons": {
            str(key): int(value)
            for key, value in unscored.loc[:, "reason"]
            .astype(str)
            .value_counts()
            .items()
        },
    }


def _sum_int(frame: pd.DataFrame, column: str) -> int:
    if column not in frame.columns:
        return 0
    return int(frame.loc[:, column].sum())


__all__ = [
    "KINASE_ATTRITION_POLICY_CAVEAT_CODE",
    "KINASE_DIRECT_TRUSTED_DATASET_CAVEAT_CODE",
    "KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE",
    "KINASE_PERMISSIVE_LOCALISATION_POLICY_CAVEAT_CODE",
    "KINASE_REFERENCE_AUTO_RESOLUTION_CAVEAT_CODE",
    "KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE",
    "KINASE_LIBRARY_MOTIF_ONLY_CAVEAT_CODE",
    "KINASE_PROFILE_SELF_INCLUSION_CAVEAT_CODE",
    "KINASE_PROFILE_LEAVE_ONE_OUT_CAVEAT_CODE",
    "KINASE_SCORING_LIMITATION_CAVEAT_CODE",
    "build_kinase_result_caveats",
]
